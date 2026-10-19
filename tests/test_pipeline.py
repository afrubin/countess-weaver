"""The plugin inside a CountESS pipeline run from its .ini file, as countess_cmd runs one."""

from __future__ import annotations

import pathlib
import shutil

import countess.core.cmd

DATA = pathlib.Path(__file__).parent / "data"
EXAMPLE = pathlib.Path(__file__).parent.parent / "examples" / "fwd-saturation"

PIPELINE = """
[Load Counts]
_module = countess.plugins.csv
_class = LoadCsvPlugin
files.0.filename = 'counts.csv'
header = True
columns.0.name = 'hgvs'
columns.0.type = 'STRING'
columns.1.name = 'count'
columns.1.type = 'INTEGER'

[Weaver]
_module = countess_weaver.plugin
_class = WeaverPlugin
_parent.0 = Load Counts
column = 'hgvs'
reference.source = 'refseq'
reference.fasta = 'synthetic.fa'
reference.gff3 = 'synthetic.gff3'
steps.0.operation = 'g_to_c'
steps.0.target = 'NM_000099.1'
steps.0.output = 'hgvs_c'
steps.1.operation = 'c_to_p'
steps.1.output = 'hgvs_p'

[Save]
_module = countess.plugins.csv
_class = SaveCsvPlugin
_parent.0 = Weaver
filename = 'scored.csv'
quoting = False
sorting.0.order_by = 'count'
"""

COUNTS = """hgvs,count
NC_000099.1:g.124G>C,7
NC_000099.1:g.128del,11
NC_000099.1:g.121C>G,13
"""


def test_a_pipeline_run_from_its_ini_file_writes_the_converted_variants(tmp_path: pathlib.Path) -> None:
    for name in ("synthetic.fa", "synthetic.gff3"):
        shutil.copy(DATA / name, tmp_path / name)
    (tmp_path / "counts.csv").write_text(COUNTS)
    ini = tmp_path / "pipeline.ini"
    ini.write_text(PIPELINE)

    countess.core.cmd.run([str(ini)])

    lines = (tmp_path / "scored.csv").read_text().splitlines()
    assert lines[0].split(",") == ["hgvs", "count", "hgvs_c", "hgvs_p", "weaver_error"]
    assert lines[1].split(",")[:4] == ["NC_000099.1:g.124G>C", "7", "NM_000099.1:c.4G>C", "NP_000099.1:p.(Ala2Pro)"]
    assert lines[2].split(",")[:3] == ["NC_000099.1:g.128del", "11", "NM_000099.1:c.14del"]
    assert lines[3].split(",")[:4] == ["NC_000099.1:g.121C>G", "13", "", ""]
    assert "ReferenceMismatchError" in lines[3]


def test_the_documented_example_collapses_the_three_spellings_of_its_deletion(tmp_path: pathlib.Path) -> None:
    example = tmp_path / "example"
    shutil.copytree(EXAMPLE, example, ignore=shutil.ignore_patterns("output"))

    countess.core.cmd.run([str(example / "fwd_saturation.ini")])

    collapsed = (example / "output" / "collapsed.csv").read_text().splitlines()
    assert '"NM_000099.1:c.14del","NP_000099.1:p.(Lys5SerfsTer?)",296.0,23.0' in collapsed
