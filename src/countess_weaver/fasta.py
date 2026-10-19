"""A weaver DataProvider over the sequences in a FASTA file, for references weaver's other
providers do not cover: an amplicon, a plasmid, a construct, a GenBank record.

It answers sequence lookups only. With no transcript models, the variants it supports are those
described on the records themselves (g. and m.); c., n., r. and p. need a reference with
transcript models.
"""

from __future__ import annotations

import gzip

import weaver


class FastaProvider:
    """The sequences of one FASTA file, held in memory and named by the first word of each header.
    Whitespace within a line is not a base, and lines starting ";" are comments.

    The whole file is read when the provider is made, so the FASTA needs no index and may have
    lines of any length; a reference large enough for that to matter belongs in one of weaver's
    indexed providers instead.
    """

    def __init__(self, path: str) -> None:
        """Reads every record of the FASTA at ``path`` (gzipped if it ends in ``.gz``).

        Raises:
            ValueError: If the file holds no records, a record has no name or no bases, or two
                records share a name.
        """
        self.path = path
        self.sequences = _read_fasta(path)

    def get_seq(self, ac: str, start: int, end: int | None, kind: str | weaver.IdentifierType) -> str:
        try:
            sequence = self.sequences[ac]
        except KeyError:
            raise weaver.DataProviderError(f"{ac}: not a record in {self.path}") from None
        # A window running past the end is answered with the bases there are; one starting past it is a
        # position the record does not have.
        if start > len(sequence):
            raise weaver.DataProviderError(
                f"{ac}: position {start + 1} is past the end of the record ({len(sequence)} bp)"
            )
        return sequence[start:end]

    def get_transcript(self, transcript_ac: str, reference_ac: str | None) -> weaver.TranscriptData:
        raise weaver.DataProviderError(
            f"{transcript_ac}: a FASTA reference has no transcript models; use a RefSeq GFF3 + FASTA reference"
        )

    def get_symbol_accessions(self, symbol: str, source_kind: str, target_kind: str) -> list[tuple[str, str]]:
        raise weaver.DataProviderError(f"{symbol}: a FASTA reference maps no symbols or accessions")

    def get_identifier_type(self, identifier: str) -> weaver.IdentifierType:
        if identifier in self.sequences:
            return weaver.IdentifierType.GenomicAccession
        return weaver.IdentifierType.Unknown


def _read_fasta(path: str) -> dict[str, str]:
    opener = gzip.open if path.endswith(".gz") else open
    sequences: dict[str, list[str]] = {}
    name = None
    with opener(path, "rt") as handle:
        for number, line in enumerate(handle, start=1):
            line = line.strip()
            if line.startswith(";"):
                continue
            if line.startswith(">"):
                words = line[1:].split()
                if not words:
                    raise ValueError(f"{path}:{number}: a FASTA header with no name")
                name = words[0]
                if name in sequences:
                    raise ValueError(f"{path}:{number}: a second record named {name}")
                sequences[name] = []
            elif line:
                if name is None:
                    raise ValueError(f"{path}:{number}: bases before the first FASTA header")
                sequences[name].append("".join(line.split()).upper())
    if not sequences:
        raise ValueError(f"{path}: no FASTA records")
    joined = {name: "".join(parts) for name, parts in sequences.items()}
    if empty := [name for name, sequence in joined.items() if not sequence]:
        raise ValueError(f"{path}: records with no bases: {', '.join(empty)}")
    return joined
