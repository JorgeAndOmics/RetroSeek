# =============================================================================
# guard.py
# =============================================================================
# The heavy-rule guard (ADR-020): read a Snakemake dry run and say whether a heavy
# rule would run without being asked for.
#
# Why this exists: LTRdigest costs about a day per genome, and several harmless-
# looking changes wake it for every genome at once (a redownloaded Pfam file, a
# pinned download URL, a reworded rule). The run protocol always said "dry-run
# first and stop if a heavy rule appears"; this module makes the launcher do it.
#
# It reads exactly two things from Snakemake's text output, both stable across
# Snakemake 7 to 9: the "Job stats" table and each job's "reason:" line. If a
# future Snakemake changes them, the guard finds no jobs and lets the run
# through; the unit tests, built on a real 9.24 dry run, flag that change.
# =============================================================================

from __future__ import annotations

import re

import stages

_RULE_LINE = re.compile(r"^(?:local)?rule (\w+):\s*$")
_REASON_LINE = re.compile(r"^\s+reason: (.+)$")


def _job_stats_sections(dry_run: str) -> list[list[str]]:
    """The lines after each "Job stats:" header, one list per header."""
    sections: list[list[str]] = []
    for line in dry_run.splitlines():
        if line.startswith("Job stats:"):
            sections.append([])
        elif sections:
            sections[-1].append(line)
    return sections


def _rule_rows(lines: list[str]) -> tuple[dict[str, int], bool]:
    """(jobs per rule, whether the table's "total" line was reached)."""
    counts: dict[str, int] = {}
    for line in lines:
        fields = line.split()
        if fields[:1] == ["total"]:
            return counts, True
        if len(fields) == 2 and fields[1].isdigit():
            counts[fields[0]] = int(fields[1])
    return counts, False


def job_counts(dry_run: str) -> dict[str, int]:
    """Jobs per rule from the first "Job stats" table of a dry run.

    Snakemake prints the table again at the end, so reading stops at the first
    table's "total" line, or at the next header once rows were counted. Empty
    when nothing would run.
    """
    counts: dict[str, int] = {}
    for section in _job_stats_sections(dry_run):
        rows, reached_total = _rule_rows(section)
        counts.update(rows)
        if reached_total or counts:
            break
    return counts


def first_reasons(dry_run: str) -> dict[str, str]:
    """The first "reason:" Snakemake gives for each rule in a dry run."""
    reasons: dict[str, str] = {}
    rule: str | None = None
    for line in dry_run.splitlines():
        if match := _RULE_LINE.match(line):
            rule = match.group(1)
        elif (match := _REASON_LINE.match(line)) and rule and rule not in reasons:
            reasons[rule] = match.group(1).strip()
    return reasons


def blocked(counts: dict[str, int], allowed: set[str]) -> list[str]:
    """Heavy rules that would run although no requested stage owns them."""
    return [r for r in counts if r in stages.HEAVY_RULES and r not in allowed]


def remedy(reason: str) -> str:
    """What to do about one reason Snakemake gave for rerunning a heavy rule."""
    lowered = reason.lower()
    if "params have changed" in lowered or "code has changed" in lowered:
        return (
            "The rule's command or settings changed, not its inputs. If its outputs "
            "are still right, clear the record once with --cleanup-metadata <its "
            "output files>, or add --rerun-triggers mtime for this run."
        )
    if "updated by another job" in lowered:
        return (
            "An upstream rule would run first and refresh this rule's inputs. Deal "
            "with that upstream rule; this one follows it."
        )
    if "updated input files" in lowered:
        return (
            "An input file is newer than the outputs (for example a redownloaded "
            "Pfam-A.hmm). If its contents did not change, backdate it: "
            "touch -d 2000-01-01 <file>."
        )
    if "missing output" in lowered:
        return (
            "Its output files are not where the config says. Check the paths in the "
            "config before running; do not let it rebuild them."
        )
    return "Read the reason above; if unsure, do not run."


def report(rules: list[str], counts: dict[str, int], reasons: dict[str, str]) -> str:
    """The message shown when the guard stops a run."""
    lines = ["These heavy rules would run although their stage was not requested:"]
    for rule in rules:
        n = counts.get(rule, 0)
        reason = reasons.get(rule, "no reason given")
        lines.append(f"  {rule} ({n} job{'s' if n != 1 else ''})")
        lines.append(f"    reason: {reason[:300]}")
        lines.append(f"    fix: {remedy(reason)}")
    lines.append(
        "Nothing was run. If this is intended, add the rule's own stage flag or "
        "--allow-heavy."
    )
    return "\n".join(lines)
