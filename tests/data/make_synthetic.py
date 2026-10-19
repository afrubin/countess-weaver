"""Writes the synthetic reference the tests use: one chromosome, NC_000099.1, carrying one
two-exon plus-strand gene, FWD (transcript NM_000099.1, protein NP_000099.1).

Layout, 1-based genomic coordinates:

    1-100    filler
    101-160  exon 1: 5' UTR 101-120, CDS 121-160 (c.1-c.40)
    161-260  intron (GT ... AG)
    261-340  exon 2: CDS 261-310 (c.41-c.90, stop codon c.88-c.90), 3' UTR 311-340
    341-399  filler
    400-406  T GGGGG C  (an intergenic homopolymer, g.401-g.405)
    407-600  filler

The CDS codes MAEKKLSGTLDRIHWPEFNVIAGRLPYDQ*; c.8-c.14 (g.128-g.134) is AAAAAAA (GAA AAA AAG).

Run from this directory to regenerate synthetic.fa and synthetic.gff3.
"""

from __future__ import annotations

import pathlib
import random

CHROM = "NC_000099.1"
LENGTH = 600
CODONS = "ATG GCT GAA AAA AAG CTG TCC GGT ACC CTG GAT CGT ATT CAC TGG CCA GAG TTC AAC GTC ATC GCA GGC CGA TTG CCT TAC GAC CAG TAA"


def filler(rng: random.Random, length: int) -> str:
    # C and T only, so filler never makes a start codon, a GGGGG run or a splice dinucleotide by chance.
    return "".join(rng.choice("CT") for _ in range(length))


def build() -> str:
    rng = random.Random(99)
    cds = CODONS.replace(" ", "")
    assert len(cds) == 90
    seq = (
        filler(rng, 100)  # 1-100
        + filler(rng, 20)  # 101-120, 5' UTR
        + cds[:40]  # 121-160
        + "GT"  # 161-260, intron
        + filler(rng, 96)
        + "AG"
        + cds[40:]  # 261-310
        + filler(rng, 30)  # 311-340, 3' UTR
        + filler(rng, 59)  # 341-399
        + "TGGGGGC"  # 400-406
        + filler(rng, 194)  # 407-600
    )
    assert len(seq) == LENGTH
    return seq


GFF3 = f"""##gff-version 3
{CHROM}\tRefSeq\tgene\t101\t340\t.\t+\t.\tID=gene-FWD;gene=FWD;Name=FWD
{CHROM}\tRefSeq\tmRNA\t101\t340\t.\t+\t.\tID=rna-NM_000099.1;Parent=gene-FWD;gene=FWD;transcript_id=NM_000099.1
{CHROM}\tRefSeq\texon\t101\t160\t.\t+\t.\tID=exon-NM_000099.1-1;Parent=rna-NM_000099.1;gene=FWD;transcript_id=NM_000099.1
{CHROM}\tRefSeq\texon\t261\t340\t.\t+\t.\tID=exon-NM_000099.1-2;Parent=rna-NM_000099.1;gene=FWD;transcript_id=NM_000099.1
{CHROM}\tRefSeq\tCDS\t121\t160\t.\t+\t0\tID=cds-NP_000099.1;Parent=rna-NM_000099.1;gene=FWD;protein_id=NP_000099.1
{CHROM}\tRefSeq\tCDS\t261\t310\t.\t+\t2\tID=cds-NP_000099.1;Parent=rna-NM_000099.1;gene=FWD;protein_id=NP_000099.1
"""


def main() -> None:
    here = pathlib.Path(__file__).parent
    seq = build()
    lines = [seq[i : i + 60] for i in range(0, len(seq), 60)]
    (here / "synthetic.fa").write_text(
        f">{CHROM} synthetic chromosome for countess-weaver tests\n" + "\n".join(lines) + "\n"
    )
    (here / "synthetic.gff3").write_text(GFF3)


if __name__ == "__main__":
    main()
