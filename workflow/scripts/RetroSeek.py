# -------------------
# DEPENDENCIES
# -------------------

import logging
import re
import subprocess
import sys
from pathlib import Path

import colored_logging
import defaults
import guard
import stages
from validator import preflight, validation_run

logger = logging.getLogger(__name__)

# -----------------------------
# FASTA EXTENSION STANDARDIZATION
# -----------------------------


def standardize_fasta_extensions(fasta_dir_path: str | Path) -> None:
    """
    Standardize extensions of all FASTA files in the provided directory to .fa.

    Parameters
    ----------
    fasta_dir_path : str
        Path to the directory containing FASTA files with various extensions (.fasta, .fna, .fas).
    """
    pattern = re.compile(r"\.(fasta|fna|fas)$", re.IGNORECASE)

    for file in Path(fasta_dir_path).iterdir():
        if file.is_file() and pattern.search(file.name):
            new_name: Path = file.with_name(f"{file.stem}.fa")
            logger.debug(f"Renaming: {file.name} -> {new_name.name}")
            file.rename(new_name)


# -----------------------------
# SNAKEMAKE RULE EXECUTION
# -----------------------------


def run_snakemake_rule(
    rule: str | list[str],
    num_cores: int,
    display_info: bool,
    snakemake_flags: list[str] | None = None,
) -> None:
    """
    Execute one or more Snakemake rules with specified options.

    Parameters
    ----------
        :param rule : str | list[str]
        Name(s) of the Snakemake rule(s) to execute. A list is passed to
        snakemake as multiple targets in a single invocation (one DAG).
        :param num_cores : int
        Number of cores to allocate for the rule.
        :param display_info : bool
        Whether to display detailed Snakemake command output.
        :param snakemake_flags:
    """
    if snakemake_flags is None:
        snakemake_flags = []
    rules = [rule] if isinstance(rule, str) else list(rule)
    rule_label = " ".join(rules)
    shell_cmd: list[str] = [
        "snakemake",
        *rules,
        "--cores",
        str(num_cores),
        "--rerun-incomplete",
        *snakemake_flags,
    ]

    if not display_info:
        shell_cmd.append("-q")

    try:
        result = subprocess.run(shell_cmd, check=False)
    except (FileNotFoundError, OSError) as exc:
        logger.error(f"Failed to invoke snakemake for rule(s) '{rule_label}': {exc}")
        sys.exit(1)

    if result.returncode != 0:
        logger.error(
            f"Snakemake rule(s) '{rule_label}' failed with exit code "
            f"{result.returncode}. See snakemake output above for details."
        )
        sys.exit(result.returncode)


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


def capture_dry_run(targets: list[str], options: list[str]) -> str:
    """A Snakemake dry run of `targets`, captured for the guard to read.

    Exits with Snakemake's own code and output if the dry run itself fails
    (a config or DAG error), since nothing sensible can be checked then.
    """
    cmd = [
        "snakemake",
        *targets,
        "--cores",
        str(defaults.NUM_CORES),
        "--rerun-incomplete",
        "-n",
        *options,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    output = result.stdout + result.stderr
    if result.returncode != 0:
        print(output)
        logger.error(f"The dry run failed with exit code {result.returncode}.")
        sys.exit(result.returncode)
    return output


# -----------------------------
# CLI ENTRYPOINT
# -----------------------------


def cli_entry() -> None:
    """
    Main entrypoint for RetroSeek CLI.

    Order: preflight (always), slow validation (unless -skp), the heavy-rule
    guard (a captured dry run), then one Snakemake call for every requested
    stage. Exits non-zero whenever something stopped the run.
    """
    colored_logging.colored_logging(log_file_name="RetroSeek_main.log")

    parser = stages.build_parser()
    args, unknown = parser.parse_known_args()
    chosen = stages.selected(args)
    if not chosen:
        parser.print_help()
        sys.exit(0)

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

    targets = stages.targets(chosen)
    options = snakemake_options(unknown, stop_on_error=args.stop_on_error)
    dry_run = is_dry_run(unknown)

    if not is_maintenance(unknown) and (dry_run or not args.allow_heavy):
        output = capture_dry_run(targets, options)
        if dry_run:
            print(output)
        counts = guard.job_counts(output)
        heavy = guard.blocked(counts, stages.allowed_heavy(chosen))
        if heavy and not args.allow_heavy:
            logger.error(guard.report(heavy, counts, guard.first_reasons(output)))
            sys.exit(1)
        if dry_run:
            return

    run_snakemake_rule(
        targets,
        num_cores=defaults.NUM_CORES,
        display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
        snakemake_flags=options,
    )
    logger.info("Finished execution.")


# -----------------------------
# EXECUTION TRIGGER
# -----------------------------

if __name__ == "__main__":
    cli_entry()
