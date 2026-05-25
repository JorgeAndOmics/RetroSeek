# -------------------
# DEPENDENCIES
# -------------------

import argparse
import logging
import re
import subprocess
import sys
from pathlib import Path

import colored_logging
import defaults
from validator import validation_run

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
# CLI ENTRYPOINT
# -----------------------------


def cli_entry() -> None:  # noqa: PLR0912, PLR0915
    """
    Main entrypoint for RetroSeek CLI.

    Handles argument parsing, input validation, and Snakemake rule dispatch.
    """
    colored_logging.colored_logging(log_file_name="RetroSeek_main.log")

    # -----------------------------
    # PARSE COMMAND-LINE ARGUMENTS
    # -----------------------------
    parser = argparse.ArgumentParser(
        description="RetroSeek: A tool for directed ERV detection and analysis."
    )

    parser.add_argument(
        "--download-genomes",
        action="store_true",
        help="Downloads genomes specified in .yaml file. "
        "Will skip download if the same files exist currently in the target directory",
    )

    parser.add_argument(
        "--download-hmm", action="store_true", help="Setup domain database"
    )

    parser.add_argument(
        "--suffix-arrays",
        action="store_true",
        help="Generate GenomeTools suffix arrays for genomes.",
    )

    parser.add_argument(
        "--ltr-candidates",
        action="store_true",
        help="Generates LTR candidate files for genomes via LTRHarvest.",
    )

    parser.add_argument(
        "--ltr-domains",
        action="store_true",
        help="Generates domain-enriched LTR candidate files for genomes via LTRDigest.",
    )

    parser.add_argument(
        "--blast-dbs",
        action="store_true",
        help="Generates BLAST databases from genomes.",
    )

    parser.add_argument(
        "--probe-extractor",
        action="store_true",
        help="Parses input CSV file and retrieves proviral sequences from provided NCBI accesion IDs.",
    )

    parser.add_argument(
        "--blast",
        action="store_true",
        help="Runs tBLASTn of provided sequences against genome databases.",
    )

    parser.add_argument(
        "--ranges-analysis",
        action="store_true",
        help="Performs genomic range integration of BLAST and LTR candidate results.",
    )

    parser.add_argument(
        "--generate-global-plots",
        action="store_true",
        help="Generate plots from the analysis.",
    )

    parser.add_argument(
        "--generate-circle-plots",
        action="store_true",
        help="Generate circle plots from the species in analysis.",
    )

    parser.add_argument(
        "--hotspot-detection",
        action="store_true",
        help="Perform ERV hotspot detection.",
    )

    parser.add_argument(
        "--pair-detection", action="store_true", help="Perform probe-pair detection."
    )

    parser.add_argument(
        "--solo-ltr-detection",
        action="store_true",
        help="Solo-LTR detection via LTR_retriever; retroviral-only pre-filter "
        "via valid_ranges; multi-label probe-label propagation; emits "
        "solo_ltr/{genome}.gff3 + solo_intact_ratio/{genome}.csv + "
        "solo_intact_ratio/all_species.csv.",
    )

    parser.add_argument(
        "--skip-validation", "-skp", action="store_true", help="Skip input validation."
    )

    args, unknown = parser.parse_known_args()

    # -----------------------------
    # STANDARDIZE FASTA EXTENSIONS
    # -----------------------------
    standardize_fasta_extensions(defaults.PATH_DICT["SPECIES_DB"])

    # -----------------------------
    # VALIDATE INPUT & EXECUTE RULES
    # -----------------------------
    species_paths: list[str] = [
        str(Path(defaults.PATH_DICT["SPECIES_DB"]) / f"{species}.fa")
        for species in defaults.SPECIES
    ]

    # Print help message if no arguments or wrong arguments are provided
    if not any(vars(args).values()):
        parser.print_help()
        sys.exit(0)

    validation = True if args.skip_validation else validation_run(species_paths)
    if validation:
        if args.download_genomes:
            run_snakemake_rule(
                "genome_downloader",
                num_cores=defaults.NUM_CORES,
                display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
                snakemake_flags=unknown,
            )

        if args.suffix_arrays:
            run_snakemake_rule(
                "ltr_index_generator",
                num_cores=defaults.NUM_CORES,
                display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
                snakemake_flags=unknown,
            )

        if args.download_hmm:
            run_snakemake_rule(
                "pfam_hmm_downloader",
                num_cores=defaults.NUM_CORES,
                display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
                snakemake_flags=unknown,
            )

        if args.ltr_candidates:
            run_snakemake_rule(
                "ltr_harvester",
                num_cores=defaults.NUM_CORES,
                display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
                snakemake_flags=unknown,
            )

        if args.ltr_domains:
            run_snakemake_rule(
                "ltr_digester",
                num_cores=defaults.NUM_CORES,
                display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
                snakemake_flags=unknown,
            )

        if args.probe_extractor:
            run_snakemake_rule(
                "probe_extractor",
                num_cores=defaults.NUM_CORES,
                display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
                snakemake_flags=unknown,
            )

        if args.blast_dbs:
            run_snakemake_rule(
                "blast_db_generator",
                num_cores=defaults.NUM_CORES,
                display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
                snakemake_flags=unknown,
            )

        if args.blast:
            run_snakemake_rule(
                "blast_pkl2parquet",
                num_cores=defaults.NUM_CORES,
                display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
                snakemake_flags=unknown,
            )

        if args.ranges_analysis:
            run_snakemake_rule(
                "ranges_analysis",
                num_cores=defaults.NUM_CORES,
                display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
                snakemake_flags=unknown,
            )

        if args.generate_global_plots:
            # --generate-global-plots produces the full panel: the provirus
            # panel (plot_generator final-tier + stage_plot_generator middle-
            # stage) plus the erv-like candidate panel (erv_like_plot_generator)
            # — one DAG, shared ranges_analysis upstream.
            run_snakemake_rule(
                ["plot_generator", "stage_plot_generator", "erv_like_plot_generator"],
                num_cores=defaults.NUM_CORES,
                display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
                snakemake_flags=unknown,
            )

        if args.generate_circle_plots:
            run_snakemake_rule(
                "circle_plot_generator",
                num_cores=defaults.NUM_CORES,
                display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
                snakemake_flags=unknown,
            )

        if args.hotspot_detection:
            run_snakemake_rule(
                "hotspot_detector",
                num_cores=defaults.NUM_CORES,
                display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
                snakemake_flags=unknown,
            )

        if args.pair_detection:
            run_snakemake_rule(
                "pair_detector",
                num_cores=defaults.NUM_CORES,
                display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
                snakemake_flags=unknown,
            )

        if args.solo_ltr_detection:
            run_snakemake_rule(
                "solo_ltr_detector",
                num_cores=defaults.NUM_CORES,
                display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
                snakemake_flags=unknown,
            )

    else:
        sys.exit(1)

    logger.info("Finished execution.")


# -----------------------------
# EXECUTION TRIGGER
# -----------------------------

if __name__ == "__main__":
    cli_entry()
