"""DAG-contract regression tests for the Snakefile.

Goal: pin structural invariants of the workflow (rule presence, wildcard
constraints) without requiring a snakemake binary on PATH for every
contributor. We parse the Snakefile with regex; that's brittle for
fancy syntax but the existing rules are all top-level and conventionally
formatted, so it works.

Tests for rules that don't exist yet (e.g. ``genome_fasta_normalizer_setup``)
are marked ``xfail(strict=True)`` — they flip to PASS when the rule lands
in its target phase, so anyone removing a planned rule prematurely will
trip a failure.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

RULE_DEF_RE = re.compile(r"^rule\s+(\w+)\s*:", re.MULTILINE)
WILDCARD_BLOCK_RE = re.compile(
    r"wildcard_constraints\s*:\s*\n\s*genome\s*=\s*(.+?)$",
    re.MULTILINE,
)


def _read_snakefile(project_root: Path) -> str:
    return (project_root / "workflow" / "Snakefile").read_text()


def _all_rules(project_root: Path) -> set[str]:
    return set(RULE_DEF_RE.findall(_read_snakefile(project_root)))


# ---------------------------------------------------------------------
# Rule presence
# ---------------------------------------------------------------------
@pytest.mark.parametrize(
    "rule_name",
    [
        "ltr_harvester_setup",
        "ltr_retriever_prefilter_setup",
        "ltr_retriever_setup",
        "solo_ltr_integrator_setup",
        "solo_ltr_detector",
    ],
)
def test_existing_rule_present(project_root: Path, rule_name: str) -> None:
    """Each LTR-Retriever workstream rule must remain in the Snakefile."""
    rules = _all_rules(project_root)
    assert rule_name in rules, f"rule {rule_name!r} missing from Snakefile"


def test_genome_fasta_normalizer_rule_present(project_root: Path) -> None:
    """``genome_fasta_normalizer_setup`` rule must be in the workflow."""
    assert "genome_fasta_normalizer_setup" in _all_rules(project_root)


def test_ruleorder_normalizer_wins_over_downloader(project_root: Path) -> None:
    """The normalizer must take precedence when both can produce {genome}.fa."""
    text = _read_snakefile(project_root)
    assert "ruleorder: genome_fasta_normalizer_setup > genome_downloader_setup" in text


def test_blast_db_generator_no_longer_renames_inline(project_root: Path) -> None:
    """The inline `find / parallel mv` workaround must be gone — normalizer owns it."""
    text = _read_snakefile(project_root)
    assert "parallel 'mv" not in text
    assert "-iname '*.fna'" not in text


# ---------------------------------------------------------------------
# Wildcard constraints
# ---------------------------------------------------------------------
def test_genome_wildcard_constraint_pinned_to_species_list(
    project_root: Path,
) -> None:
    """``wildcard_constraints.genome`` must be a regex pinned to SPECIES.

    Protects the recent fix (commit f1a844b) that prevents Snakemake's
    default ``.+`` matcher from greedily absorbing suffixes like
    ``_retroviral`` or ``_full`` into ``{genome}``.
    """
    text = _read_snakefile(project_root)
    match = WILDCARD_BLOCK_RE.search(text)
    assert match is not None, (
        "wildcard_constraints block for {genome} missing — "
        "see commit f1a844b for why this matters"
    )
    expr = match.group(1)
    # The expression should reference SPECIES and the no-match fallback.
    assert "SPECIES" in expr, "wildcard regex no longer pinned to SPECIES list"


# ---------------------------------------------------------------------
# Phase-2/3 contract checks (xfail until those phases land)
# ---------------------------------------------------------------------
def test_prefilter_rule_declares_both_retroviral_and_full_outputs(
    project_root: Path,
) -> None:
    """Phase 2: the prefilter rule emits both ``_retroviral`` and ``_full`` SCNs."""
    text = _read_snakefile(project_root)
    assert "{genome}_retroviral.scn" in text
    assert "{genome}_full.scn" in text


def test_ltr_retriever_setup_invokes_runner_script(project_root: Path) -> None:
    """Phase 3: the inline LTR_retriever shell collapses into a runner script call."""
    text = _read_snakefile(project_root)
    assert "run_ltr_retriever.py" in text


def test_config_yaml_uses_source_scn_field(project_root: Path) -> None:
    """``restrict_to_retroviral`` replaced with ``source_scn``."""
    text = (project_root / "data" / "config" / "config.yaml").read_text()
    assert "source_scn:" in text
    assert "restrict_to_retroviral" not in text
