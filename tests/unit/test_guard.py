"""Unit tests for workflow/scripts/guard.py, the heavy-rule guard.

The guard reads a Snakemake dry run and refuses to start when a heavy rule (a
day-per-genome search) would run without its own stage being asked for. The
fixture is a real Snakemake 9.24 dry run from 2026-09-23 (paths replaced): pinning
the Pfam URL made the downloader rerun and LTRdigest follow on every genome.
"""

from pathlib import Path

import pytest

import guard

FIXTURE = (
    Path(__file__).parents[1] / "fixtures" / "snakemake" / "dry_run_pfam_cascade.txt"
)


@pytest.fixture
def dry_run() -> str:
    return FIXTURE.read_text(encoding="utf-8")


def test_job_counts_reads_the_job_stats_table(dry_run: str) -> None:
    counts = guard.job_counts(dry_run)
    assert counts["pfam_hmm_downloader"] == 1
    assert counts["ltr_digester_setup"] == 5
    assert counts["taxonomy_classify_setup"] == 5
    assert "total" not in counts
    assert len(counts) == 15


def test_job_counts_takes_the_first_table_only() -> None:
    """Snakemake prints the table again at the end of a dry run."""
    text = (
        "Job stats:\njob  count\n---  -----\nrule_a  2\ntotal  2\n\n"
        "Job stats:\njob  count\n---  -----\nrule_b  9\ntotal  9\n"
    )
    assert guard.job_counts(text) == {"rule_a": 2}


def test_job_counts_is_empty_when_nothing_would_run() -> None:
    assert guard.job_counts("Building DAG of jobs...\nNothing to be done.\n") == {}


def test_first_reasons_keeps_one_reason_per_rule(dry_run: str) -> None:
    reasons = guard.first_reasons(dry_run)
    assert reasons["pfam_hmm_downloader"].startswith("Params have changed")
    assert reasons["ltr_digester_setup"].startswith(
        "Input files updated by another job"
    )


def test_first_reasons_reads_local_rules_too() -> None:
    text = "localrule solo_blaster_setup:\n    jobid: 3\n    reason: Missing output files: x\n"
    assert guard.first_reasons(text) == {
        "solo_blaster_setup": "Missing output files: x"
    }


def test_blocked_lists_heavy_rules_nobody_asked_for(dry_run: str) -> None:
    counts = guard.job_counts(dry_run)
    assert guard.blocked(counts, allowed=set()) == [
        "pfam_hmm_downloader",
        "ltr_digester_setup",
    ]


def test_blocked_lets_a_requested_heavy_stage_through(dry_run: str) -> None:
    counts = guard.job_counts(dry_run)
    allowed = {"pfam_hmm_downloader", "ltr_digester_setup", "ltr_digester"}
    assert guard.blocked(counts, allowed=allowed) == []


@pytest.mark.parametrize(
    ("reason", "expected"),
    [
        ("Params have changed since last execution: ...", "--cleanup-metadata"),
        ("Code has changed since last execution", "--rerun-triggers mtime"),
        ("Input files updated by another job: /a/Pfam-A.hmm", "upstream"),
        ("Updated input files: /a/Pfam-A.hmm", "backdate"),
        ("Missing output files: /a/b.gff3", "paths"),
    ],
)
def test_remedy_matches_the_reason(reason: str, expected: str) -> None:
    assert expected in guard.remedy(reason)


def test_report_names_rules_counts_reasons_and_the_override(dry_run: str) -> None:
    counts = guard.job_counts(dry_run)
    reasons = guard.first_reasons(dry_run)
    text = guard.report(["pfam_hmm_downloader", "ltr_digester_setup"], counts, reasons)
    assert "ltr_digester_setup (5 jobs)" in text
    assert "pfam_hmm_downloader (1 job)" in text
    assert "Params have changed" in text
    assert "--allow-heavy" in text
