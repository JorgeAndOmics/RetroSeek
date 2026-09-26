# =============================================================================
# RetroSeek.py: the pipeline side of the launcher (ADR-020, ADR-021)
# =============================================================================
# Called by the root `./RetroSeek` shim. One run goes:
#
#   banner -> checks (preflight always, validation unless -skp) -> the heavy-rule
#   guard (a captured dry run) -> ONE Snakemake call streamed through console.py
#   -> summary. The exit code is Snakemake's, or 1 when a check or the guard
#   stopped the run.
#
# Every line of the run, ours and Snakemake's, lands in LOG_DIR/runs/<time>.log;
# the screen shows what `display.verbosity` (or --verbosity) asks for.
# =============================================================================

"""The pipeline side of the `./RetroSeek` launcher (ADR-020, ADR-021)."""

import argparse
import logging
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import console
import defaults
import guard
import stages
from log import OK
from validator import preflight, uses_pfam, validation_run

logger = logging.getLogger(__name__)

# -----------------------------
# FASTA EXTENSION STANDARDIZATION
# -----------------------------


def standardize_fasta_extensions(fasta_dir_path: str | Path) -> None:
    """Rename every .fasta, .fna or .fas file in a directory to .fa.

    The match ignores case. Each file is renamed in place, so an existing
    file with the new name is replaced.

    Args:
        fasta_dir_path: The directory holding the FASTA files.
    """
    pattern = re.compile(r"\.(fasta|fna|fas)$", re.IGNORECASE)

    for file in Path(fasta_dir_path).iterdir():
        if file.is_file() and pattern.search(file.name):
            new_name: Path = file.with_name(f"{file.stem}.fa")
            logger.debug(f"Renaming: {file.name} -> {new_name.name}")
            file.rename(new_name)


# -----------------------------
# SNAKEMAKE OPTIONS
# -----------------------------

_DRY_RUN_FLAGS = ("-n", "--dry-run", "--dryrun")
# Modes that run no jobs. A guard dry run carrying them would itself unlock or
# clean, so they go straight to Snakemake.
_MAINTENANCE_FLAGS = ("--unlock", "--cleanup-metadata", "--cm")


def snakemake_options(user_options: list[str], stop_on_error: bool) -> list[str]:
    """The options every Snakemake call gets, before the user's own.

    `--keep-going` is the default so one failed genome does not stop the other
    genomes' jobs overnight. The user's options come last and keep their order,
    because `--forcerun RULE...` swallows everything after it.
    """
    if stop_on_error:
        return list(user_options)
    return ["--keep-going", *user_options]


def is_dry_run(user_options: list[str]) -> bool:
    """Whether the user asked Snakemake for a dry run."""
    return any(flag in _DRY_RUN_FLAGS for flag in user_options)


def is_maintenance(user_options: list[str]) -> bool:
    """Whether the call is an --unlock or --cleanup-metadata, which runs no jobs."""
    return any(flag in _MAINTENANCE_FLAGS for flag in user_options)


def snakemake_command(targets: list[str], cores: int, options: list[str]) -> list[str]:
    """The Snakemake command line for `targets`: the same for the run and the guard.

    Snakemake always speaks in full; the launcher decides what reaches the screen
    and keeps everything in the run log.
    """
    return [
        "snakemake",
        *targets,
        "--cores",
        str(cores),
        "--rerun-incomplete",
        *options,
    ]


def capture_dry_run(targets: list[str], options: list[str]) -> tuple[int, str]:
    """A Snakemake dry run of `targets`, captured for the guard: (exit code, text)."""
    cmd = snakemake_command(targets, defaults.NUM_CORES, ["-n", *options])
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    return result.returncode, result.stdout + result.stderr


# -----------------------------
# WHAT THE RUN LOOKS LIKE
# -----------------------------


def _commit() -> str:
    """The checked-out commit, so a run log says which code produced it."""
    result = subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"],
        capture_output=True,
        text=True,
        check=False,
        cwd=Path(__file__).resolve().parent,
    )
    return result.stdout.strip() or "unknown commit"


def banner(chosen: list[stages.Stage], verbosity: str, run_log: Path) -> list[str]:
    """The lines that open a run: what runs, on what, and where its record goes."""
    names = " > ".join(s.flag.lstrip("-").replace("-", " ") for s in chosen)
    return [
        f"RetroSeek {_commit()}",
        f"  config     {defaults.CONFIG_FILE}",
        f"  genomes    {len(defaults.SPECIES)}    cores  {defaults.NUM_CORES}"
        f"    verbosity  {verbosity}",
        f"  stages     {names}",
        f"  run log    {run_log}",
    ]


def _append(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as out:
        out.write(text if text.endswith("\n") else text + "\n")


def write_warnings(tally: console.Tally, path: Path) -> None:
    """Every warning of the run, one line each, next to the run log."""
    _append(
        path,
        "".join(
            f"{e.time} WARN {e.step} {e.genome} | {e.message}\n" for e in tally.warnings
        ),
    )


# -----------------------------
# CLI ENTRYPOINT
# -----------------------------


def run_checks(args: argparse.Namespace, chosen: list[stages.Stage]) -> None:
    """Preflight (always) and the slow validation (unless -skp); exit 1 on failure."""
    standardize_fasta_extensions(defaults.PATH_DICT["SPECIES_DB"])
    if not preflight(chosen):
        logger.error("Preflight checks failed; nothing was run.")
        sys.exit(1)
    species_paths = [
        str(Path(defaults.PATH_DICT["SPECIES_DB"]) / f"{species}.fa")
        for species in defaults.SPECIES
    ]
    if not args.skip_validation and not validation_run(chosen, species_paths):
        sys.exit(1)
    logger.log(
        OK, "checks passed: config, tools%s", ", Pfam" if uses_pfam(chosen) else ""
    )


def run_guard(
    args: argparse.Namespace,
    chosen: list[stages.Stage],
    targets: list[str],
    options: list[str],
    screen: console.Screen,
    run_log: Path,
) -> None:
    """The heavy-rule guard. Exits when it stops the run, or after a dry run."""
    code, output = capture_dry_run(targets, options)
    _append(run_log, output)
    dry_run = is_dry_run(options)
    if code != 0 or dry_run:
        screen.lines(output.splitlines())
    if code != 0:
        logger.error(f"The dry run failed (exit {code}); nothing was run.")
        sys.exit(code)
    counts = guard.job_counts(output)
    heavy = guard.blocked(counts, stages.allowed_heavy(chosen))
    if heavy and not args.allow_heavy:
        logger.error(guard.report(heavy, counts, guard.first_reasons(output)))
        sys.exit(1)
    logger.log(
        OK, "guard: %s jobs to run, none of them heavy", f"{sum(counts.values()):,}"
    )
    if dry_run:
        logger.warning("Dry run only: nothing was run.")
        sys.exit(0)


def run_workflow(
    targets: list[str],
    options: list[str],
    screen: console.Screen,
    run_log: Path,
    started: float,
) -> int:
    """The one Snakemake call, streamed; then the summary. Returns its exit code."""
    screen.heading("Run")
    tally = console.Tally()
    cmd = snakemake_command(targets, defaults.NUM_CORES, options)
    code = console.stream(cmd, screen, tally, run_log)

    warnings_file = run_log.with_suffix(".warnings.txt")
    if tally.warnings:
        write_warnings(tally, warnings_file)
    status = {0: "done", 130: "interrupted"}.get(code, f"failed (exit {code})")
    screen.heading("Summary")
    summary = console.summary_lines(
        tally,
        status,
        time.monotonic() - started,
        str(run_log),
        str(warnings_file) if tally.warnings else None,
    )
    screen.lines(summary)
    _append(run_log, "\n".join(summary))
    return code


def cli_entry() -> None:
    """Main entrypoint for RetroSeek CLI.

    Order: banner, checks, the heavy-rule guard, then one Snakemake call for
    every requested stage and the summary. Exits non-zero whenever something
    stopped the run.
    """
    parser = stages.build_parser()
    args, unknown = parser.parse_known_args()
    chosen = stages.selected(args)
    if not chosen:
        parser.print_help()
        sys.exit(0)

    verbosity = args.verbosity or defaults.VERBOSITY
    # Every job inherits the environment, so log.py and log.R filter the same way.
    os.environ["RETROSEEK_VERBOSITY"] = verbosity
    started = time.monotonic()
    stamp = time.strftime("%Y-%m-%d_%H%M%S")
    run_log = Path(defaults.PATH_DICT["LOG_DIR"]) / "runs" / f"{stamp}.log"
    screen = console.Screen(verbosity)
    root = logging.getLogger()
    root.handlers = [console.ScreenHandler(screen, run_log)]
    root.setLevel(logging.DEBUG)

    screen.lines(banner(chosen, verbosity, run_log))
    screen.heading("Checks")
    run_checks(args, chosen)

    targets = stages.targets(chosen)
    options = snakemake_options(unknown, stop_on_error=args.stop_on_error)
    if not is_maintenance(unknown) and (is_dry_run(unknown) or not args.allow_heavy):
        run_guard(args, chosen, targets, options, screen, run_log)

    sys.exit(run_workflow(targets, options, screen, run_log, started))


# -----------------------------
# EXECUTION TRIGGER
# -----------------------------

if __name__ == "__main__":
    cli_entry()
