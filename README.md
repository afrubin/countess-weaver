# countess-weaver

A [CountESS](https://github.com/CountESS-Project/CountESS) plugin that applies
[hgvs-weaver](https://github.com/populationgenomics/hgvs-weaver) to a column of HGVS variant descriptions:
normalising them by the 3' rule and projecting them between genomic (`g.`), coding (`c.`), non-coding (`n.`), RNA
(`r.`) and protein (`p.`) descriptions.

## The HGVS Weaver plugin

Each value of the **input column** is parsed by weaver and passed through a chain of **steps**. A step is one weaver
operation; any step can write its result to an output column, and the last step must. Each variant is converted on its
own, and each distinct value only once, however many rows carry it.

| Operation    | Target accession                          |
|--------------|-------------------------------------------|
| Normalise    | none                                      |
| g. → c.      | the transcript (required)                 |
| c. → g.      | the genomic reference (optional)          |
| n. → g.      | the genomic reference (optional)          |
| r. → g.      | the genomic reference (optional)          |
| c. → r., n. → r., r. → c., r. → n. | none                |
| c./r. → p.   | the protein (optional)                    |

A target is either fixed for every row or read from a column. Projections do not normalise their result on their own;
**Normalise Result**, on by default, normalises each nucleotide result.

The **reference** is one of:

- **FASTA (sequences only)**: any FASTA file — an amplicon, a construct, a GenBank record. It has no transcript models,
  so it supports the variants described on its own records (`g.`, `m.`): normalisation, not projection.
- **RefSeq GFF3 + FASTA**: weaver's RefSeq provider, with transcript models for projection and protein prediction.
  Needs the `refseq` extra (`pip install countess-weaver[refseq]`). That pulls in the biocommons `hgvs` package, which
  depends on `psycopg2`; weaver's provider never uses it, and where it will not build, leave it out of the install
  (with uv, an override such as `psycopg2; sys_platform == 'never'`).

When weaver refuses a variant — it does not parse, a base it states is not the reference's, it is the wrong kind of
variant for a step — that row's outputs from the failing step on are NULL, and **On Error** says what else happens:
an error column holding the step and weaver's message (the default), a logged warning, both, or neither.

## Development

```shell
pip install -e ".[refseq,dev]"
pytest
```

The tests run without network access against a synthetic reference, `tests/data/synthetic.fa` and
`tests/data/synthetic.gff3`, written by `tests/data/make_synthetic.py`.
