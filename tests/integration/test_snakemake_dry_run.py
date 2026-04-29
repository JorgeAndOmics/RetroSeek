"""Integration test: Snakemake dry-run resolves the DAG against the default config.

Skipped automatically when ``snakemake`` is not on PATH so the unit suite
still runs in minimal environments. Real fixture-based end-to-end runs
will live alongside this file once fixture genomes are available.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest


@pytest.mark.integration
def test_snakemake_dry_run_resolves_dag(project_root: Path) -> None:
    """Snakemake should build the rule graph without attempting to execute.

    Snakemake 8 requires an explicit target (no implicit default rule);
    we pass ``genome_downloader``, an aggregate rule that resolves over
    all configured ``SPECIES`` and exercises the early DAG without
    attempting any heavy work (since ``-n`` short-circuits execution).
    """
    if shutil.which("snakemake") is None:
        pytest.skip("snakemake not installed; activate the conda env first")

    result = subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(project_root / "workflow" / "Snakefile"),
            "--configfile",
            str(project_root / "data" / "config" / "config.yaml"),
            "-n",
            "--cores",
            "1",
            "genome_downloader",
        ],
        cwd=project_root,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, (
        f"snakemake dry-run failed:\n--- stdout ---\n{result.stdout}\n"
        f"--- stderr ---\n{result.stderr}"
    )
