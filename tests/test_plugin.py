"""The HGVS Weaver plugin, driven through its configuration keys as a CountESS .ini sets them.

Expected values are worked out from the synthetic reference's layout (tests/data/make_synthetic.py):
c.1 is g.121, the CDS codes MAEKK..., c.8-c.14 (g.128-g.134) is a run of seven As, and g.401-g.405
is a run of five Gs.
"""

from __future__ import annotations

import logging
import pathlib
import shutil
from typing import Any

import duckdb
import pandas as pd
import pytest

from countess_weaver import plugin as weaver_plugin

DATA = pathlib.Path(__file__).parent / "data"


@pytest.fixture
def reference_dir(tmp_path: pathlib.Path) -> pathlib.Path:
    # A copy: the RefSeq provider indexes the FASTA, writing a .fai beside it.
    for name in ("synthetic.fa", "synthetic.gff3"):
        shutil.copy(DATA / name, tmp_path / name)
    return tmp_path


def configure(reference_dir: pathlib.Path, config: dict[str, Any]) -> weaver_plugin.WeaverPlugin:
    plugin = weaver_plugin.WeaverPlugin()
    for key, value in config.items():
        plugin.set_parameter(key, value, str(reference_dir))
    return plugin


def run(plugin: weaver_plugin.WeaverPlugin, rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Runs the plugin over rows (each with a distinct "id"), returning its output rows in id order."""
    ddbc = duckdb.connect()
    ddbc.from_df(pd.DataFrame(rows)).create("n_0")
    source = ddbc.table("n_0")
    plugin.prepare(ddbc, source)
    result = plugin.execute(ddbc, source)
    return sorted((dict(zip(result.columns, row)) for row in result.fetchall()), key=lambda row: row["id"])


FASTA_REFERENCE = {"reference.source": "fasta", "reference.fasta": "synthetic.fa", "column": "hgvs"}
REFSEQ_REFERENCE = {
    "reference.source": "refseq",
    "reference.fasta": "synthetic.fa",
    "reference.gff3": "synthetic.gff3",
    "column": "hgvs",
}
NORMALISE = {"steps.0.operation": "normalize", "steps.0.output": "normalised"}


@pytest.mark.parametrize(
    ("hgvs", "normalised"),
    [
        ("NC_000099.1:g.401del", "NC_000099.1:g.405del"),
        ("NC_000099.1:g.400_401insG", "NC_000099.1:g.405dup"),
        ("NC_000099.1:g.401_402del", "NC_000099.1:g.404_405del"),
        ("NC_000099.1:g.121A>G", "NC_000099.1:g.121A>G"),
    ],
)
def test_a_fasta_reference_normalises_a_genomic_variant_by_the_3_prime_rule(
    reference_dir: pathlib.Path, hgvs: str, normalised: str
) -> None:
    plugin = configure(reference_dir, FASTA_REFERENCE | NORMALISE)

    [row] = run(plugin, [{"id": 1, "hgvs": hgvs}])

    assert row["normalised"] == normalised
    assert row["weaver_error"] is None


def test_every_row_carrying_a_variant_gets_its_result(reference_dir: pathlib.Path) -> None:
    plugin = configure(reference_dir, FASTA_REFERENCE | NORMALISE)
    rows = [
        {"id": 1, "hgvs": "NC_000099.1:g.401del", "count": 10},
        {"id": 2, "hgvs": "NC_000099.1:g.402del", "count": 20},
        {"id": 3, "hgvs": "NC_000099.1:g.401del", "count": 30},
    ]

    result = run(plugin, rows)

    assert [(r["id"], r["count"], r["normalised"]) for r in result] == [
        (1, 10, "NC_000099.1:g.405del"),
        (2, 20, "NC_000099.1:g.405del"),
        (3, 30, "NC_000099.1:g.405del"),
    ]


def test_a_null_input_gives_null_outputs_and_no_error(reference_dir: pathlib.Path) -> None:
    plugin = configure(reference_dir, FASTA_REFERENCE | NORMALISE)

    result = run(plugin, [{"id": 1, "hgvs": None}, {"id": 2, "hgvs": "NC_000099.1:g.401del"}])

    assert (result[0]["normalised"], result[0]["weaver_error"]) == (None, None)
    assert result[1]["normalised"] == "NC_000099.1:g.405del"


def test_a_step_can_write_over_its_input_column(reference_dir: pathlib.Path) -> None:
    plugin = configure(reference_dir, FASTA_REFERENCE | {"steps.0.operation": "normalize", "steps.0.output": "hgvs"})

    [row] = run(plugin, [{"id": 1, "hgvs": "NC_000099.1:g.401del"}])

    assert row["hgvs"] == "NC_000099.1:g.405del"


def test_a_base_the_reference_does_not_carry_is_reported_in_the_error_column(reference_dir: pathlib.Path) -> None:
    plugin = configure(reference_dir, FASTA_REFERENCE | NORMALISE)

    [row] = run(plugin, [{"id": 1, "hgvs": "NC_000099.1:g.121C>G"}])

    assert row["normalised"] is None
    assert "ReferenceMismatchError" in row["weaver_error"]


def test_an_unparseable_description_is_reported_as_a_parse_failure(reference_dir: pathlib.Path) -> None:
    plugin = configure(reference_dir, FASTA_REFERENCE | NORMALISE)

    [row] = run(plugin, [{"id": 1, "hgvs": "not a variant"}])

    assert row["normalised"] is None
    assert row["weaver_error"].startswith("parse: ")


def test_an_error_is_written_on_one_line(reference_dir: pathlib.Path) -> None:
    plugin = configure(reference_dir, FASTA_REFERENCE | NORMALISE)

    [row] = run(plugin, [{"id": 1, "hgvs": "NC_000099.1:g.124G>"}])

    assert "\n" not in row["weaver_error"]
    assert "expected dna" in row["weaver_error"]


def test_a_variant_on_a_sequence_the_fasta_lacks_is_reported(reference_dir: pathlib.Path) -> None:
    plugin = configure(reference_dir, FASTA_REFERENCE | NORMALISE)

    [row] = run(plugin, [{"id": 1, "hgvs": "NC_000098.1:g.401del"}])

    assert row["normalised"] is None
    assert "NC_000098.1" in row["weaver_error"]


def test_a_projection_from_a_fasta_reference_is_refused_as_having_no_transcripts(reference_dir: pathlib.Path) -> None:
    plugin = configure(
        reference_dir,
        FASTA_REFERENCE | {"steps.0.operation": "g_to_c", "steps.0.target": "NM_000099.1", "steps.0.output": "c"},
    )

    [row] = run(plugin, [{"id": 1, "hgvs": "NC_000099.1:g.124G>C"}])

    assert row["c"] is None
    assert "no transcript models" in row["weaver_error"]


def test_a_chain_writes_each_steps_result(reference_dir: pathlib.Path) -> None:
    plugin = configure(
        reference_dir,
        REFSEQ_REFERENCE
        | {
            "steps.0.operation": "g_to_c",
            "steps.0.target": "NM_000099.1",
            "steps.0.output": "c",
            "steps.1.operation": "c_to_p",
            "steps.1.output": "p",
        },
    )

    [row] = run(plugin, [{"id": 1, "hgvs": "NC_000099.1:g.124G>C"}])

    assert (row["c"], row["p"]) == ("NM_000099.1:c.4G>C", "NP_000099.1:p.(Ala2Pro)")


def test_a_projection_normalises_its_result_by_default(reference_dir: pathlib.Path) -> None:
    plugin = configure(
        reference_dir,
        REFSEQ_REFERENCE | {"steps.0.operation": "g_to_c", "steps.0.target": "NM_000099.1", "steps.0.output": "c"},
    )

    [row] = run(plugin, [{"id": 1, "hgvs": "NC_000099.1:g.128del"}])

    assert row["c"] == "NM_000099.1:c.14del"


def test_a_projection_without_normalise_result_keeps_the_position_projected(reference_dir: pathlib.Path) -> None:
    plugin = configure(
        reference_dir,
        REFSEQ_REFERENCE
        | {
            "steps.0.operation": "g_to_c",
            "steps.0.target": "NM_000099.1",
            "steps.0.normalise": False,
            "steps.0.output": "c",
        },
    )

    [row] = run(plugin, [{"id": 1, "hgvs": "NC_000099.1:g.128del"}])

    assert row["c"] == "NM_000099.1:c.8del"


def test_a_target_accession_can_come_from_a_column(reference_dir: pathlib.Path) -> None:
    plugin = configure(
        reference_dir,
        REFSEQ_REFERENCE | {"steps.0.operation": "g_to_c", "steps.0.target": "— transcript", "steps.0.output": "c"},
    )
    rows = [
        {"id": 1, "hgvs": "NC_000099.1:g.124G>C", "transcript": "NM_000099.1"},
        {"id": 2, "hgvs": "NC_000099.1:g.124G>C", "transcript": "NM_000098.1"},
    ]

    result = run(plugin, rows)

    assert [r["id"] for r in result] == [1, 2]
    assert result[0]["c"] == "NM_000099.1:c.4G>C"
    assert result[1]["c"] is None
    assert "NM_000098.1" in result[1]["weaver_error"]


def test_a_failing_step_keeps_the_results_of_the_steps_before_it(reference_dir: pathlib.Path) -> None:
    plugin = configure(
        reference_dir,
        REFSEQ_REFERENCE | NORMALISE | {"steps.1.operation": "c_to_p", "steps.1.output": "p"},
    )

    [row] = run(plugin, [{"id": 1, "hgvs": "NC_000099.1:g.128del"}])

    assert row["normalised"] == "NC_000099.1:g.134del"
    assert row["p"] is None
    assert row["weaver_error"].startswith("step 2 ")


@pytest.mark.parametrize(
    ("mode", "has_error_column", "logs"),
    [
        ("column", True, False),
        ("log", False, True),
        ("both", True, True),
        ("silent", False, False),
    ],
)
def test_the_error_mode_chooses_where_a_failure_is_reported(
    reference_dir: pathlib.Path,
    caplog: pytest.LogCaptureFixture,
    mode: str,
    has_error_column: bool,
    logs: bool,
) -> None:
    plugin = configure(reference_dir, FASTA_REFERENCE | NORMALISE | {"on_error": mode})

    with caplog.at_level(logging.WARNING, logger=weaver_plugin.__name__):
        [row] = run(plugin, [{"id": 1, "hgvs": "NC_000099.1:g.121C>G"}])

    assert row["normalised"] is None
    assert ("weaver_error" in row) == has_error_column
    assert any("NC_000099.1:g.121C>G" in r.getMessage() for r in caplog.records) == logs


def test_the_error_column_is_named_by_its_parameter(reference_dir: pathlib.Path) -> None:
    plugin = configure(reference_dir, FASTA_REFERENCE | NORMALISE | {"error_column": "why"})

    [row] = run(plugin, [{"id": 1, "hgvs": "NC_000099.1:g.121C>G"}])

    assert "ReferenceMismatchError" in row["why"]


@pytest.mark.parametrize(
    ("config", "message"),
    [
        ({"steps.0.operation": "g_to_c", "steps.0.output": "c"}, "needs a transcript accession"),
        (
            {"steps.0.operation": "normalize", "steps.0.target": "NM_000099.1", "steps.0.output": "n"},
            "takes no target accession",
        ),
        ({"steps.0.operation": "normalize"}, "needs an output column"),
        (
            {
                "steps.0.operation": "normalize",
                "steps.0.output": "x",
                "steps.1.operation": "normalize",
                "steps.1.output": "x",
            },
            "used more than once",
        ),
        ({}, "at least one step"),
    ],
)
def test_a_configuration_that_cannot_run_is_refused(
    reference_dir: pathlib.Path, config: dict[str, Any], message: str
) -> None:
    plugin = configure(reference_dir, REFSEQ_REFERENCE | config)

    with pytest.raises(ValueError, match=message):
        run(plugin, [{"id": 1, "hgvs": "NC_000099.1:g.401del"}])


def test_a_reference_with_no_fasta_is_refused(reference_dir: pathlib.Path) -> None:
    plugin = configure(reference_dir, NORMALISE | {"reference.source": "fasta", "column": "hgvs"})

    with pytest.raises(ValueError, match="FASTA"):
        run(plugin, [{"id": 1, "hgvs": "NC_000099.1:g.401del"}])


def test_a_changed_reference_file_is_read_again(reference_dir: pathlib.Path) -> None:
    plugin = configure(reference_dir, FASTA_REFERENCE | NORMALISE)
    run(plugin, [{"id": 1, "hgvs": "NC_000099.1:g.401del"}])
    # g.2 is the second of two Ts in the synthetic reference, the second of four As here.
    (reference_dir / "synthetic.fa").write_text(">NC_000099.1\nAAAACCCC\n")

    [row] = run(plugin, [{"id": 1, "hgvs": "NC_000099.1:g.2del"}])

    assert row["normalised"] == "NC_000099.1:g.4del"


def test_each_distinct_input_is_converted_once(reference_dir: pathlib.Path, monkeypatch: pytest.MonkeyPatch) -> None:
    converted = []
    convert = weaver_plugin.convert

    def counting_convert(mapper: Any, steps: Any, key: dict[str, Any]) -> weaver_plugin.Conversion:
        converted.append(key["hgvs"])
        return convert(mapper, steps, key)

    monkeypatch.setattr(weaver_plugin, "convert", counting_convert)
    plugin = configure(reference_dir, FASTA_REFERENCE | NORMALISE)
    rows = [{"id": n, "hgvs": f"NC_000099.1:g.{401 + n % 2}del"} for n in range(10)]

    run(plugin, rows)

    assert sorted(converted) == ["NC_000099.1:g.401del", "NC_000099.1:g.402del"]


def test_an_input_column_the_input_lacks_is_refused_and_kept(reference_dir: pathlib.Path) -> None:
    plugin = configure(reference_dir, FASTA_REFERENCE | NORMALISE)

    with pytest.raises(ValueError, match="'hgvs' is not a column"):
        run(plugin, [{"id": 1, "variant": "NC_000099.1:g.401del"}])
    assert ("column", "hgvs") in plugin.get_parameters("")


def test_a_target_column_the_input_lacks_is_refused_and_kept(reference_dir: pathlib.Path) -> None:
    plugin = configure(
        reference_dir,
        REFSEQ_REFERENCE | {"steps.0.operation": "c_to_g", "steps.0.target": "— genome", "steps.0.output": "g"},
    )

    with pytest.raises(ValueError, match="'genome' is not a column"):
        run(plugin, [{"id": 1, "hgvs": "NM_000099.1:c.4G>C"}])
    assert ("steps.0.target", "— genome") in plugin.get_parameters("")


def test_an_output_column_differing_only_in_case_replaces_the_existing_column(reference_dir: pathlib.Path) -> None:
    plugin = configure(reference_dir, FASTA_REFERENCE | {"steps.0.operation": "normalize", "steps.0.output": "HGVS"})

    [row] = run(plugin, [{"id": 1, "hgvs": "NC_000099.1:g.401del"}])

    assert row == {"id": 1, "HGVS": "NC_000099.1:g.405del", "weaver_error": None}


def test_output_columns_differing_only_in_case_are_refused(reference_dir: pathlib.Path) -> None:
    plugin = configure(
        reference_dir,
        FASTA_REFERENCE
        | {
            "steps.0.operation": "normalize",
            "steps.0.output": "Norm",
            "steps.1.operation": "normalize",
            "steps.1.output": "norm",
        },
    )

    with pytest.raises(ValueError, match="used more than once"):
        run(plugin, [{"id": 1, "hgvs": "NC_000099.1:g.401del"}])


def test_a_rerun_leaves_no_more_tables_than_the_first_run(reference_dir: pathlib.Path) -> None:
    plugin = configure(reference_dir, FASTA_REFERENCE | NORMALISE)
    ddbc = duckdb.connect()
    ddbc.from_df(pd.DataFrame([{"id": 1, "hgvs": "NC_000099.1:g.401del"}])).create("n_0")
    source = ddbc.table("n_0")
    plugin.prepare(ddbc, source)
    plugin.execute(ddbc, source).to_table("n_1")
    [(after_first,)] = ddbc.sql("SELECT count(*) FROM duckdb_tables()").fetchall()

    ddbc.sql("DROP TABLE n_1")
    plugin.execute(ddbc, source).to_table("n_1")

    assert ddbc.sql("SELECT count(*) FROM duckdb_tables()").fetchall() == [(after_first,)]
