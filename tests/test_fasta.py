"""The FASTA provider, through weaver's mapper: what a FASTA reference can and cannot be."""

from __future__ import annotations

import gzip
import pathlib

import pytest
import weaver

from countess_weaver import fasta


def normalise(provider: fasta.FastaProvider, hgvs: str) -> str:
    return weaver.VariantMapper(provider).normalize_variant(weaver.parse(hgvs)).format()


def test_a_fasta_with_uneven_line_lengths_is_read(tmp_path: pathlib.Path) -> None:
    # Lines of different lengths cannot be faidx-indexed; pebbles' test reference is one such file.
    path = tmp_path / "amplicon.fa"
    path.write_text(">AMP1 an amplicon\nTTG\nGGGCTT\nC\n")

    assert normalise(fasta.FastaProvider(str(path)), "AMP1:g.3del") == "AMP1:g.6del"


def test_a_gzipped_fasta_is_read(tmp_path: pathlib.Path) -> None:
    path = tmp_path / "amplicon.fa.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(">AMP1\nTTGGGGCTT\n")

    assert normalise(fasta.FastaProvider(str(path)), "AMP1:g.3del") == "AMP1:g.6del"


def test_lowercase_bases_are_read_as_bases(tmp_path: pathlib.Path) -> None:
    path = tmp_path / "soft_masked.fa"
    path.write_text(">AMP1\nttggggctt\n")

    assert normalise(fasta.FastaProvider(str(path)), "AMP1:g.3G>A") == "AMP1:g.3G>A"


def test_each_record_is_its_own_reference(tmp_path: pathlib.Path) -> None:
    path = tmp_path / "two.fa"
    path.write_text(">AMP1\nTTGGGGCTT\n>AMP2\nTTAAAACTT\n")

    assert normalise(fasta.FastaProvider(str(path)), "AMP2:g.3del") == "AMP2:g.6del"


@pytest.mark.parametrize(
    ("content", "message"),
    [
        ("", "no FASTA records"),
        (">AMP1\nTTG\n>AMP1\nCCA\n", "a second record named AMP1"),
        ("TTG\n>AMP1\nCCA\n", "bases before the first FASTA header"),
        (">AMP1\n>AMP2\nCCA\n", "records with no bases: AMP1"),
        (">\nCCA\n", "a FASTA header with no name"),
    ],
)
def test_a_malformed_fasta_is_refused(tmp_path: pathlib.Path, content: str, message: str) -> None:
    path = tmp_path / "bad.fa"
    path.write_text(content)

    with pytest.raises(ValueError, match=message):
        fasta.FastaProvider(str(path))


def test_whitespace_within_a_line_is_not_a_base(tmp_path: pathlib.Path) -> None:
    path = tmp_path / "spaced.fa"
    path.write_text(">AMP1\nTTG GGG CTT\n")

    assert normalise(fasta.FastaProvider(str(path)), "AMP1:g.7C>A") == "AMP1:g.7C>A"


def test_a_comment_line_is_not_a_base(tmp_path: pathlib.Path) -> None:
    path = tmp_path / "commented.fa"
    path.write_text(">AMP1\n; a comment\nTTGGGGCTT\n")

    assert normalise(fasta.FastaProvider(str(path)), "AMP1:g.7C>A") == "AMP1:g.7C>A"


def test_a_sequence_window_starting_past_the_end_of_the_record_is_refused(tmp_path: pathlib.Path) -> None:
    # weaver reads in blocks of 4096 bases, so on a 9-base record g.5000 asks for a window starting at 4096.
    # A position past the end inside the record's first block is weaver's to refuse, not the provider's.
    path = tmp_path / "amplicon.fa"
    path.write_text(">AMP1\nTTGGGGCTT\n")

    with pytest.raises(weaver.DataProviderError, match="past the end"):
        normalise(fasta.FastaProvider(str(path)), "AMP1:g.5000del")
