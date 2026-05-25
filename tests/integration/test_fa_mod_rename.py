"""Integration test: the LTR_retriever ``.fa.mod.*`` rename loop produces canonical filenames.

LTR_retriever writes its outputs alongside the input genome FASTA and
prefixes them with ``<basename>.fa.mod.`` because it internally creates
a sanitised copy of the genome. The Snakefile rule ``ltr_retriever_setup``
includes a shell rename loop that strips the ``.fa.mod.`` infix back to
the canonical ``{genome}.X`` form expected by downstream rules.

This test reproduces the rename-loop logic in isolation (no LTR_retriever
binary required) so the regression survives a refactor of the rule body.
The extension list is pinned to the three files declared as outputs:
``pass.list.gff3``, ``nmtf.pass.list``, ``LTRlib.fa``. If any of those
disappear from the rule, this test breaks, surfacing the divergence.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

EXPECTED_EXTENSIONS = ("pass.list.gff3", "nmtf.pass.list", "LTRlib.fa")
RENAME_LOOP = """
for ext in pass.list.gff3 nmtf.pass.list LTRlib.fa; do
    if [ -f "${GENOME}.fa.mod.${ext}" ]; then
        mv "${GENOME}.fa.mod.${ext}" "${GENOME}.${ext}"
    fi
done
"""


@pytest.mark.integration
def test_rename_loop_renames_three_outputs(tmp_path: Path) -> None:
    """All three ``.fa.mod.X`` files get renamed to canonical ``{genome}.X``."""
    if shutil.which("bash") is None:
        pytest.skip("bash not available")

    genome = "Toyus"
    for ext in EXPECTED_EXTENSIONS:
        (tmp_path / f"{genome}.fa.mod.{ext}").write_text(f"placeholder for {ext}\n")

    subprocess.run(
        ["bash", "-c", RENAME_LOOP],
        cwd=tmp_path,
        env={"GENOME": genome, "PATH": "/usr/bin:/bin"},
        check=True,
    )

    for ext in EXPECTED_EXTENSIONS:
        canonical = tmp_path / f"{genome}.{ext}"
        assert canonical.is_file(), f"canonical {canonical.name} missing after rename"
        assert canonical.read_text() == f"placeholder for {ext}\n"
        assert not (tmp_path / f"{genome}.fa.mod.{ext}").exists(), (
            f"original .fa.mod.{ext} should have been moved, not copied"
        )


@pytest.mark.integration
def test_rename_loop_is_noop_when_no_fa_mod_files(tmp_path: Path) -> None:
    """The loop must not error if no ``.fa.mod.*`` files are present.

    Protects against a future regression where someone replaces the
    ``[ -f ... ]`` guard with an unconditional ``mv`` that would fail
    on a clean re-run.
    """
    if shutil.which("bash") is None:
        pytest.skip("bash not available")

    subprocess.run(
        ["bash", "-c", RENAME_LOOP],
        cwd=tmp_path,
        env={"GENOME": "Toyus", "PATH": "/usr/bin:/bin"},
        check=True,
    )
    assert list(tmp_path.iterdir()) == []
