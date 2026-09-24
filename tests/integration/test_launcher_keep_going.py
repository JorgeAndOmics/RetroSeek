"""The launcher's failure summary against a real Snakemake run (ADR-020, ADR-021).

A toy workflow in a temporary folder has one job per "genome"; one of them fails.
With --keep-going the other finishes, and the tally read from Snakemake's own
output names the failed rule, its genome, its log and the script's last error:
exactly what the closing summary shows. This pins the Snakemake 9 text the
launcher relies on, which unit fixtures alone could let drift.
"""

import io
import shutil
import sys
from pathlib import Path

import pytest

import console

SNAKEFILE = """
rule all:
    input: expand("out/{genome}.txt", genome=["Good_genome", "Bad_genome"])

rule make:
    output: "out/{genome}.txt"
    log: "logs/make/{genome}.log"
    shell:
        '''
        {python} -c "import sys; g = '{wildcards.genome}'; \\
print('10:00:00 OK make ' + g + ' | fine') if g.startswith('Good') else \\
(print('10:00:00 ERROR make ' + g + ' | the input is empty. Fix: rerun the stage before', file=sys.stderr), sys.exit(1))" \\
            && touch {output}
        '''
"""


@pytest.mark.skipif(shutil.which("snakemake") is None, reason="needs snakemake")
def test_keep_going_and_the_failure_summary(tmp_path: Path) -> None:
    (tmp_path / "Snakefile").write_text(SNAKEFILE.replace("{python}", sys.executable))
    out = io.StringIO()
    screen = console.Screen(file=out, verbosity="normal")
    tally = console.Tally()
    cmd = [
        "snakemake",
        "all",
        "--cores",
        "1",
        "--keep-going",
        "--directory",
        str(tmp_path),
        "--snakefile",
        str(tmp_path / "Snakefile"),
    ]

    code = console.stream(cmd, screen, tally, tmp_path / "runs" / "r.log")

    assert code != 0
    assert (tmp_path / "out" / "Good_genome.txt").exists()  # the other job finished
    assert [(f.rule, f.genome) for f in tally.failures] == [("make", "Bad_genome")]
    assert tally.failures[0].log.endswith("logs/make/Bad_genome.log")
    summary = "\n".join(console.summary_lines(tally, "failed", 1.0, "r.log"))
    assert "make Bad_genome: the input is empty. Fix: rerun the stage before" in summary
    assert "fine" in out.getvalue()  # the OK line of the good job reached the screen
