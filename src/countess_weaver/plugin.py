"""The HGVS Weaver CountESS plugin: each value of an HGVS column is parsed by weaver and passed
through a chain of steps (normalise, project between coordinate systems, predict the protein
consequence), with the result of any step written to a column.

Every distinct input is converted once, independently of every other, and the results are joined
back onto the rows, so a count table that repeats a variant across many rows costs one conversion
per variant.
"""

from __future__ import annotations

import dataclasses
import importlib.util
import logging
import os
import secrets
from collections.abc import Callable, Sequence
from typing import Any, Optional

import duckdb
import duckdb.func as duckdb_func
import weaver
from countess.core.parameters import (
    ArrayParam,
    BooleanParam,
    ColumnChoiceParam,
    ColumnOrStringParam,
    DictChoiceParam,
    FileParam,
    FramedMultiParam,
    StringParam,
    TabularMultiParam,
)
from countess.core.plugins import DuckdbSimplePlugin
from countess.utils.duckdb import duckdb_escape_identifier, duckdb_source_to_table

from countess_weaver import VERSION
from countess_weaver.fasta import FastaProvider

logger = logging.getLogger(__name__)

# Individual failures logged per run in the "log" error modes; the rest are counted in one line.
MAX_LOGGED_FAILURES = 20

Apply = Callable[[weaver.VariantMapper, weaver.Variant, Optional[str]], weaver.Variant]


@dataclasses.dataclass(frozen=True)
class Operation:
    label: str
    apply: Apply
    # What the step's target accession names, or None if the operation takes no target.
    target: Optional[str] = None
    target_required: bool = False
    # A p. result is already written as weaver predicts it; normalisation applies to nucleotide results.
    normalisable: bool = True


OPERATIONS: dict[str, Operation] = {
    "normalize": Operation("Normalise", lambda m, v, t: m.normalize_variant(v), normalisable=False),
    "g_to_c": Operation("g. → c.", lambda m, v, t: m.g_to_c(v, t), target="transcript", target_required=True),
    "c_to_g": Operation("c. → g.", lambda m, v, t: m.c_to_g(v, t), target="genomic reference"),
    "n_to_g": Operation("n. → g.", lambda m, v, t: m.n_to_g(v, t), target="genomic reference"),
    "r_to_g": Operation("r. → g.", lambda m, v, t: m.r_to_g(v, t), target="genomic reference"),
    "c_to_r": Operation("c. → r.", lambda m, v, t: m.c_to_r(v)),
    "n_to_r": Operation("n. → r.", lambda m, v, t: m.n_to_r(v)),
    "r_to_c": Operation("r. → c.", lambda m, v, t: m.r_to_c(v)),
    "r_to_n": Operation("r. → n.", lambda m, v, t: m.r_to_n(v)),
    "c_to_p": Operation("c./r. → p.", lambda m, v, t: m.c_to_p(v, t), target="protein", normalisable=False),
}

REFERENCE_SOURCES = {
    "fasta": "FASTA (sequences only)",
    "refseq": "RefSeq GFF3 + FASTA",
}

ERROR_MODES = {
    "column": "Error column",
    "log": "Log warning",
    "both": "Error column and log warning",
    "silent": "Neither (failures are only NULL outputs)",
}

FASTA_FILE_TYPES = [("FASTA", ["*.fa", "*.fasta", "*.fna", "*.fa.gz", "*.fasta.gz", "*.fna.gz"]), ("Any", "*")]
GFF3_FILE_TYPES = [("GFF3", ["*.gff", "*.gff3", "*.gff.gz", "*.gff3.gz"]), ("Any", "*")]


class ReferenceParam(FramedMultiParam):
    source = DictChoiceParam("Source", choices=REFERENCE_SOURCES)
    fasta = FileParam("FASTA File", file_types=FASTA_FILE_TYPES)
    gff3 = FileParam("GFF3 File (RefSeq only)", file_types=GFF3_FILE_TYPES)


class StepParam(TabularMultiParam):
    operation = DictChoiceParam("Operation", choices={key: op.label for key, op in OPERATIONS.items()})
    target = ColumnOrStringParam("Target Accession")
    normalise = BooleanParam("Normalise Result", True)
    output = StringParam("Output Column")


@dataclasses.dataclass(frozen=True)
class Step:
    number: int
    operation: Operation
    target: Optional[str]
    target_column: Optional[str]
    normalise: bool
    output: Optional[str]

    @property
    def name(self) -> str:
        return f"step {self.number} ({self.operation.label})"


@dataclasses.dataclass(frozen=True)
class Conversion:
    """One input's outputs, by step (None where the step was not reached or failed), and its failure."""

    outputs: list[Optional[str]]
    error: Optional[str]


class WeaverPlugin(DuckdbSimplePlugin):
    """Applies weaver's HGVS operations to a column of variant descriptions."""

    name = "HGVS Weaver"
    description = "Normalise HGVS variants and project them between g., c., n., r. and p. with weaver"
    link = "https://github.com/populationgenomics/hgvs-weaver"
    version = VERSION

    column = ColumnChoiceParam("Input Column")
    reference = ReferenceParam("Reference")
    steps = ArrayParam("Steps", StepParam("Step"))
    on_error = DictChoiceParam("On Error", choices=ERROR_MODES)
    error_column = StringParam("Error Column", "weaver_error")

    _mapper: Optional[weaver.VariantMapper] = None
    _mapper_key: Optional[tuple] = None
    _tables: tuple[str, ...] = ()

    def prepare(self, ddbc: duckdb.DuckDBPyConnection, source: Optional[duckdb.DuckDBPyRelation]) -> None:
        # CountESS swaps a configured column the input lacks for another (or blanks it); restore the
        # configured names so that execute refuses them and a saved pipeline keeps them.
        column = self.column.value
        targets = [param.target.value for param in self.steps]
        super().prepare(ddbc, source)
        if column:
            self.column.value = column
        for param, target in zip(self.steps, targets):
            if isinstance(target, str) and target.startswith(ColumnOrStringParam.PREFIX):
                param.target.value = target

    def execute(
        self, ddbc: duckdb.DuckDBPyConnection, source: duckdb.DuckDBPyRelation, row_limit: Optional[int] = None
    ) -> Optional[duckdb.DuckDBPyRelation]:
        steps = self.configured_steps(source.columns)
        error_column = self.error_column.value if self.on_error.choice in ("column", "both") else None
        check_output_columns(steps, error_column)
        mapper = self.mapper()

        # The previous run's tables: its result has been copied out by now, or the run was abandoned.
        for name in self._tables:
            ddbc.sql(f"DROP TABLE IF EXISTS {name}")
        self._tables = ()

        # Materialised once: the distinct inputs and the join below would otherwise each run the upstream query.
        table = duckdb_source_to_table(ddbc, source)
        self._tables = (table.alias,)
        key_columns = [self.column.value] + sorted({s.target_column for s in steps if s.target_column})
        results = results_table(ddbc, table, key_columns, mapper, steps)
        self._tables = (table.alias, results.alias)

        if self.on_error.choice in ("log", "both"):
            log_failures(results)
        return join_results(ddbc, table, results, key_columns, steps, error_column)

    def configured_steps(self, columns: Sequence[str]) -> list[Step]:
        """The steps as configured, checked against the source's columns.

        Raises:
            ValueError: If the configuration cannot run: no input column or steps, a target given
                to an operation that takes none or missing from one that needs it, or no output
                from the last step.
        """
        if not self.column.value:
            raise ValueError("Choose an input column")
        if self.column.value not in columns:
            raise ValueError(f"Input column {self.column.value!r} is not a column of the input")
        if not self.steps:
            raise ValueError("Add at least one step")

        steps = []
        for number, param in enumerate(self.steps, start=1):
            operation = OPERATIONS[param.operation.choice]
            target_column = param.target.get_column_name()
            if target_column is not None and target_column not in columns:
                raise ValueError(f"Step {number}: target column {target_column!r} is not a column of the input")
            target = None if target_column else (param.target.get_string_value() or None)
            if operation.target is None and (target or target_column):
                raise ValueError(f"Step {number} ({operation.label}) takes no target accession")
            if operation.target_required and not (target or target_column):
                raise ValueError(f"Step {number} ({operation.label}) needs a {operation.target} accession")
            steps.append(
                Step(
                    number=number,
                    operation=operation,
                    target=target,
                    target_column=target_column,
                    normalise=bool(param.normalise.value) and operation.normalisable,
                    output=param.output.value or None,
                )
            )
        if steps[-1].output is None:
            raise ValueError(f"The last step (step {steps[-1].number}) needs an output column")
        return steps

    def mapper(self) -> weaver.VariantMapper:
        """The mapper over the configured reference, kept while the reference files are unchanged."""
        key = self.reference_key()
        if self._mapper is None or key != self._mapper_key:
            self._mapper = weaver.VariantMapper(self.provider())
            self._mapper_key = key
        return self._mapper

    def reference_key(self) -> tuple:
        source = self.reference.source.choice
        files = [self.reference.fasta] + ([self.reference.gff3] if source == "refseq" else [])
        for param in files:
            if not param.value:
                raise ValueError(f"Reference: choose a {param.label}")
        return (source, *(file_identity(param.get_file_path()) for param in files))

    def provider(self) -> Any:  # a weaver DataProvider
        fasta_path = self.reference.fasta.get_file_path()
        if self.reference.source.choice == "fasta":
            return FastaProvider(fasta_path)
        return refseq_provider(self.reference.gff3.get_file_path(), fasta_path)


def refseq_provider(gff3_path: str, fasta_path: str) -> Any:
    # weaver.cli.provider exits the interpreter when hgvs or pysam is missing; check first.
    missing = [name for name in ("hgvs", "pysam") if importlib.util.find_spec(name) is None]
    if missing:
        raise ImportError(
            f"A RefSeq GFF3 + FASTA reference needs {' and '.join(missing)}: " "install countess-weaver[refseq]"
        )
    from weaver.cli import provider  # pylint: disable=import-outside-toplevel

    return provider.RefSeqDataProvider(gff3_path, fasta_path)


def file_identity(path: str) -> tuple[str, int, int]:
    stat = os.stat(path)
    return (path, stat.st_mtime_ns, stat.st_size)


def check_output_columns(steps: Sequence[Step], error_column: Optional[str]) -> None:
    if error_column == "":
        raise ValueError("Choose a name for the error column")
    names = [s.output for s in steps if s.output] + ([error_column] if error_column is not None else [])
    # DuckDB column names are case-insensitive.
    folded = [n.casefold() for n in names]
    duplicates = sorted({n for n in names if folded.count(n.casefold()) > 1})
    if duplicates:
        raise ValueError(f"Output column names used more than once: {', '.join(duplicates)}")


def results_table(
    ddbc: duckdb.DuckDBPyConnection,
    table: duckdb.DuckDBPyRelation,
    key_columns: Sequence[str],
    mapper: weaver.VariantMapper,
    steps: Sequence[Step],
) -> duckdb.DuckDBPyRelation:
    """A table of each distinct key (input and target columns, as __k<n>) with each step's output
    (__o<n>) and the error (__error).

    The conversion runs inside DuckDB, as a function over the distinct keys, so that a query
    interrupt (CountESS's GUI interrupts a preview it is replacing) stops it.
    """
    outputs = [f"__o{n}" for n in range(len(steps))] + ["__error"]
    return_type = duckdb.struct_type({name: duckdb.sqltype("VARCHAR") for name in outputs})

    def convert_key(key: list[Optional[str]]) -> dict[str, Optional[str]]:
        conversion = convert(mapper, steps, dict(zip(key_columns, key)))
        return dict(zip(outputs, conversion.outputs + [conversion.error]))

    function_name = "weaver_" + secrets.token_hex(16)
    ddbc.create_function(
        function_name,
        convert_key,
        [duckdb.list_type(duckdb.sqltype("VARCHAR"))],
        return_type,
        null_handling=duckdb_func.FunctionNullHandling.SPECIAL,
    )
    try:
        keys = ", ".join(f"CAST({duckdb_escape_identifier(c)} AS VARCHAR) AS __k{n}" for n, c in enumerate(key_columns))
        key_names = ", ".join(f"__k{n}" for n in range(len(key_columns)))
        input_id = duckdb_escape_identifier(key_columns[0])
        relation = ddbc.sql(
            f"SELECT {key_names}, UNNEST({function_name}([{key_names}]))"
            f" FROM (SELECT DISTINCT {keys} FROM {table.alias} WHERE {input_id} IS NOT NULL)"
        )
        return duckdb_source_to_table(ddbc, relation)
    finally:
        ddbc.remove_function(function_name)


def convert(mapper: weaver.VariantMapper, steps: Sequence[Step], key: dict[str, Optional[str]]) -> Conversion:
    """Runs one input through the steps; the first step to fail ends the chain."""
    outputs: list[Optional[str]] = [None] * len(steps)
    hgvs = next(iter(key.values()))
    try:
        variant = weaver.parse(hgvs)
    except weaver.HGVSError as exc:
        return Conversion(outputs, f"parse: {describe(exc)}")

    for index, step in enumerate(steps):
        target = key[step.target_column] if step.target_column else step.target
        if step.operation.target_required and not target:
            return Conversion(outputs, f"{step.name}: no {step.operation.target} accession in this row")
        try:
            variant = step.operation.apply(mapper, variant, target or None)
            if step.normalise:
                variant = mapper.normalize_variant(variant)
        # ValueError: weaver's answer to an operation given the wrong kind of variant (a c. step given a g.).
        except (weaver.HGVSError, ValueError) as exc:
            return Conversion(outputs, f"{step.name}: {describe(exc)}")
        outputs[index] = variant.format()
    return Conversion(outputs, None)


def describe(exc: Exception) -> str:
    # One line for a table cell: weaver's parse errors draw a caret diagram over several lines.
    return f"{type(exc).__name__}: {' '.join(str(exc).split())}"


def log_failures(results: duckdb.DuckDBPyRelation) -> None:
    failures = results.filter("__error IS NOT NULL")
    for hgvs, error in failures.project("__k0, __error").order("__k0").limit(MAX_LOGGED_FAILURES).fetchall():
        logger.warning("%s: %s", hgvs, error)
    [(count,)] = failures.count("*").fetchall()
    if count > MAX_LOGGED_FAILURES:
        logger.warning("%d more variants failed", count - MAX_LOGGED_FAILURES)


def join_results(
    ddbc: duckdb.DuckDBPyConnection,
    table: duckdb.DuckDBPyRelation,
    results: duckdb.DuckDBPyRelation,
    key_columns: Sequence[str],
    steps: Sequence[Step],
    error_column: Optional[str],
) -> duckdb.DuckDBPyRelation:
    outputs = {f"__o{n}": s.output for n, s in enumerate(steps) if s.output}
    if error_column is not None:
        outputs["__error"] = error_column

    folded = {column.casefold() for column in outputs.values()}
    replaced = [c for c in table.columns if c.casefold() in folded]
    exclude = f" EXCLUDE ({', '.join(duckdb_escape_identifier(c) for c in replaced)})" if replaced else ""
    selects = ", ".join(f"r.{name} AS {duckdb_escape_identifier(column)}" for name, column in outputs.items())
    conditions = [f"CAST(s.{duckdb_escape_identifier(key_columns[0])} AS VARCHAR) = r.__k0"] + [
        f"CAST(s.{duckdb_escape_identifier(c)} AS VARCHAR) IS NOT DISTINCT FROM r.__k{n}"
        for n, c in enumerate(key_columns[1:], start=1)
    ]
    sql = (
        f"SELECT s.*{exclude}, {selects} FROM {table.alias} AS s"
        f" LEFT JOIN {results.alias} AS r ON {' AND '.join(conditions)}"
    )
    logger.debug("WeaverPlugin sql %s", sql)
    return ddbc.sql(sql)
