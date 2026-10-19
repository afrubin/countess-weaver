# Example: annotating and collapsing saturation-editing counts

This example takes variant counts from a small saturation genome editing experiment on a synthetic gene, FWD. It
normalises each variant, describes it on the transcript and the protein, and sums the counts of variants that turn out
to be the same change.

Open it in the GUI from this directory:

```shell
countess_gui fwd_saturation.ini
```

or run it without the GUI:

```shell
countess_cmd fwd_saturation.ini
```

The results are written to `output/annotated.csv` and `output/collapsed.csv`.

## The data

- `synthetic.fa` and `synthetic.gff3` are the reference: a 600-base synthetic chromosome, `NC_000099.1`, carrying one
  plus-strand gene, FWD. Its transcript `NM_000099.1` has two exons (genomic 101–160 and 261–340). Its 30-codon CDS
  starts at g.121 and codes `MAEKKLSGTLDRIHWPEFNVIAGRLPYDQ*`. The sequence is made up for testing and documentation;
  it is not from any organism.
- `counts.csv` has one row per variant, written in genomic HGVS the way a read-level variant caller reports it, with
  counts before (`pre`) and after (`post`) selection.

Three of the rows, `g.128del`, `g.131del` and `g.134del`, are one change written three ways. The CDS has a run of seven
As (g.128–g.134), and deleting any one of them gives the same sequence. Aligners place such a deletion at different
positions in the run, often the 5'-most, but HGVS describes it at the 3'-most. Two more rows are mistakes, there to show
how weaver refuses a variant: `g.124C>G` states a reference base the sequence does not have (g.124 is a G), and
`g.124G>` is not valid HGVS.

## The pipeline

**Load Counts** reads `counts.csv`.

**HGVS Weaver** reads the `hgvs` column against the RefSeq GFF3 + FASTA reference and runs three steps on each variant:

1. **Normalise** writes the 3'-most genomic description to `hgvs_g`. All three spellings of the deletion become
   `NC_000099.1:g.134del`.
2. **g. → c.** projects it onto `NM_000099.1` and writes `hgvs_c`, normalised on the transcript (`c.14del`).
   Intronic positions are numbered from the nearest exon (`g.161G>A` is `c.40+1G>A`, a splice donor).
3. **c./r. → p.** predicts the protein consequence and writes `hgvs_p`: missense `p.(Ala2Pro)`, synonymous
   `p.(Ala2=)`, nonsense `p.(Glu3Ter)`, start loss `p.(Met1?)`, stop loss `p.(Ter30GlnextTer?)`, an in-frame deletion
   `p.(Lys5del)`, and `p.?` for the splice-region and intronic variants, whose effect on the protein cannot be read from
   the DNA.

The two mistaken rows get no descriptions. Their `weaver_error` column names the step that refused them and why.

**Save Annotated** writes every input row with its three descriptions.

**Collapse by Variant** groups the rows by `hgvs_c` and `hgvs_p` and sums the counts. The three deletion rows collapse
into one `c.14del` row with 296 reads before selection and 23 after. The refused rows group together under empty
descriptions.

**Save Collapsed** writes the grouped table.

## Things to try in the GUI

- Change **On Error** to *Log warning* and watch the `weaver_error` column disappear and the reasons move to the log.
- Delete the **Normalise** step: the transcript descriptions are unchanged, because **Normalise Result** on the
  g. → c. step normalises the projection.
- Then also turn off **Normalise Result**: each deletion keeps the position the aligner gave it (`c.8del`, `c.11del`,
  `c.14del`), and the collapsed table keeps them as three rows, though all three predict the same protein change.
