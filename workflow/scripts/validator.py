"""Input checks the launcher runs before a pipeline run.

The preflight (config schema, tools, Pfam library) always runs; the slower
validation (NCBI probe lookups, genome FASTAs) and the confirmation prompt are
skipped with `-skp`.
"""

# -------------------
# DEPENDENCIES
# -------------------

import logging
import os
import shutil
import sys
import time
from pathlib import Path

import pandas as pd
import pandera as pa
import yamale
from Bio import Entrez, SeqIO

import defaults
import stages

# The Pfam check reuses the subset builder's own accession logic, so the preflight
# and the rule cannot disagree about what "missing" means.
sys.path.insert(0, str(Path(__file__).resolve().parent / "domains"))
from subset_pfam import missing_accessions, missing_message, wanted_accessions

logger = logging.getLogger(__name__)

# Stages whose inputs come from NCBI: the slow probe check and the API-key prompt
# only serve these.
NCBI_STAGES = ("--download-genomes", "--probe-extractor", "--build-reference")
# Stages that read the curated Pfam subset: the preflight checks the library first.
PFAM_STAGES = ("--domain-scan", "--classify")

# -----------------------------
# YAML VALIDATION
# -----------------------------

# Keys that older configs carry, and what replaced them. The strict schema rejects
# them anyway; this turns "unexpected key" into something a user can act on.
RETIRED_KEYS: dict[tuple[str, ...], str] = {
    (
        "display",
        "display_snakemake_info",
    ): "display.verbosity (quiet | normal | verbose)",
    (
        "display",
        "display_operation_info",
    ): "display.verbosity: verbose shows the detail",
    ("display", "display_requests_warning"): (
        "display.verbosity: retries show at verbose, a failed fetch is always an error"
    ),
    ("logging",): "nothing: the terminal colours are fixed now (ADR-021)",
    (
        "plots",
        "circle_plot_bitscore_threshold",
    ): "nothing: the circle-plot stage was removed",
}


def _value_at(config: dict[str, object], path: tuple[str, ...]) -> object:
    """The value at a key path, or None when any key on the way is missing."""
    node: object = config
    for key in path:
        node = node.get(key) if isinstance(node, dict) else None
        if node is None:
            return None
    return node


def retired_key_messages(config: dict[str, object]) -> list[str]:
    """One message per retired key present in `config`, naming its replacement.

    A retired key set to null (a bare `key:` in YAML) counts as absent.
    """
    return [
        f"Config key `{'.'.join(path)}` was retired; delete it. Replaced by: "
        f"{replacement}. See docs/configuration.md."
        for path, replacement in RETIRED_KEYS.items()
        if _value_at(config, path) is not None
    ]


def yaml_validator(yaml_file: str | Path, yaml_schema: str) -> bool:
    """Validates the YAML configuration file against a Yamale schema.

    Parameters
    ----------
    yaml_file : str
        Path to the YAML configuration file.
    yaml_schema : str
        Path to the YAML schema file.

    Returns:
    -------
    bool
        True if the YAML is valid, False otherwise.
    """
    try:
        schema = yamale.make_schema(yaml_schema)
        data = yamale.make_data(yaml_file)
        retired = retired_key_messages(data[0][0])
        for message in retired:
            logger.error(message)
        if retired:
            return False
        yamale.validate(schema, data)
        logger.info("YAML configuration file is valid.")
        return True

    except yamale.YamaleError as e:
        for results in e.results:
            for error in results.errors:
                logger.warning(f"Configuration error: {error}")
        return False


# -----------------------------
# CSV VALIDATION
# -----------------------------


def csv_validator(csv_file: str) -> bool:
    """Validates the CSV file for required fields and checks NCBI accession validity.

    Parameters
    ----------
    csv_file : str
        Path to the CSV file.

    Returns:
    -------
    bool
        True if CSV passes validation, False otherwise.
    """

    def check_ncbi(series: pd.Series) -> pd.Series:
        """Check if accession IDs exist in the NCBI protein database."""
        logger.debug("Checking NCBI entries...")

        def ncbi_exists(acc: str) -> bool:
            time.sleep(0.3)
            Entrez.email = defaults.ENTREZ_EMAIL
            try:
                with Entrez.esearch(db="protein", term=acc) as handle:  # type: ignore[no-untyped-call]
                    record = Entrez.read(handle)  # type: ignore[no-untyped-call]
                    return int(record["Count"]) > 0
            except Exception:
                return False

        return series.apply(ncbi_exists)

    try:
        df = pd.read_csv(csv_file)

        schema = pa.DataFrameSchema(  # type: ignore[no-untyped-call]
            {
                "Label": pa.Column(str, nullable=False),
                "Name": pa.Column(str, nullable=False),
                "Abbreviation": pa.Column(str, nullable=False),
                "Probe": pa.Column(str, nullable=False),
                "Accession": pa.Column(
                    str, nullable=False, checks=pa.Check(check_ncbi, element_wise=False)
                ),
            }
        )

        schema.validate(df)
        logger.info("CSV input file is valid.")
        return True

    except (pa.errors.SchemaError, FileNotFoundError, KeyError) as e:
        logger.warning(f"CSV input error: {e}")
        return False


# -----------------------------
# FASTA VALIDATION
# -----------------------------


def fasta_validator(fasta_file: str) -> bool:
    """Validates a FASTA file for content and headers.

    Parameters
    ----------
    fasta_file : str
        Path to the FASTA file.

    Returns:
    -------
    bool
        True if valid, False otherwise.
    """
    if not Path(fasta_file).exists():
        logger.warning(f"FASTA file does not exist: {fasta_file}")
        return False

    try:
        records = list(SeqIO.parse(fasta_file, "fasta"))  # type: ignore[no-untyped-call]
        if not records:
            logger.warning("FASTA file is empty or has no valid records.")
            return False

        for i, record in enumerate(records):
            if not record.id:
                logger.warning(f"Record {i + 1} is missing a header.")
                return False

        logger.info(f"FASTA file {fasta_file} is valid.")
        return True

    except Exception as e:
        logger.warning(f"Error parsing FASTA: {e}")
        return False


# -----------------------------
# PREFLIGHT: FAST CHECKS THAT ALWAYS RUN
# -----------------------------


def missing_tools(tools: list[str]) -> list[str]:
    """The executables in `tools` that are not on PATH."""
    return [tool for tool in tools if shutil.which(tool) is None]


def pfam_problem(hmm_path: Path, classes_tsv: Path) -> str | None:
    """Why the Pfam library cannot serve the curated table, or None if it can.

    A library that is not there yet is not a problem: the pinned release will be
    downloaded (and the heavy-rule guard asks for --download-hmm first).
    """
    if not hmm_path.exists():
        return None
    missing = missing_accessions(hmm_path, wanted_accessions(classes_tsv))
    return missing_message(missing, hmm_path) if missing else None


def uses_pfam(chosen: list[stages.Stage]) -> bool:
    """Whether any chosen stage reads the curated Pfam subset."""
    return any(stage.flag in PFAM_STAGES for stage in chosen)


def preflight(chosen: list[stages.Stage]) -> bool:
    """Checks that take seconds and save hours; `-skp` does not skip them.

    The config against its schema, the tools the chosen stages call, and, for
    stages that read the curated Pfam subset, the Pfam library itself.
    """
    ok = yaml_validator(
        yaml_schema=str(Path(defaults.PATH_DICT["CONFIG_DIR"]) / "schema.yaml"),
        yaml_file=defaults.CONFIG_FILE,
    )

    absent = missing_tools(["snakemake", *stages.tools(chosen)])
    if absent:
        logger.error(
            f"Not installed: {', '.join(absent)}. Activate the RetroSeek conda "
            "environment, or add them with `make env-update`."
        )
        ok = False

    if uses_pfam(chosen):
        logger.info("Checking the Pfam library against the curated table...")
        problem = pfam_problem(
            Path(defaults.PATH_DICT["HMM_PROFILE_DIR"]) / "Pfam-A.hmm",
            Path(defaults.PFAM_DOMAIN_CLASSES),
        )
        if problem:
            logger.error(problem)
            ok = False

    return ok


# -----------------------------
# INTERACTIVE PROMPTS
# -----------------------------


def ask(question: str, default: str = "") -> str:
    """Asks the user a question, falling back to `default` when nobody answers.

    Unattended runs (CI, an agent, `nohup`) have no terminal attached, so a
    bare `input()` raises EOFError and takes the whole pipeline down before
    Snakemake is ever reached. An empty answer means the same thing as no
    answer at all: use the default.

    Parameters
    ----------
    question : str
        Prompt shown to the user.
    default : str
        Answer to assume when the user just hits enter, or when there is no
        terminal to ask.

    Returns:
    -------
    str
        The user's answer, or `default`.
    """
    try:
        return input(question) or default
    except EOFError:
        logger.info("No terminal attached. Continuing with the default answer.")
        return default


# -----------------------------
# NCBI API KEY VALIDATION
# -----------------------------


def validate_ncbi_key() -> None:
    """Make sure an NCBI API key is set in the environment, prompting if it is not.

    An empty answer skips it with a warning: NCBI retrievals then run slower.
    """
    if "NCBI_API_KEY" not in os.environ:
        logger.warning(
            'No NCBI API key found. You can set it with: export NCBI_API_KEY="your_api_key".'
        )
        if api_key := ask("Enter your NCBI API key [Leave empty to skip]: "):
            os.environ["NCBI_API_KEY"] = api_key
            logger.info("NCBI API key set in environment variables.")
        else:
            logger.warning("No NCBI API key provided. Expect slower NCBI retrievals.")
    else:
        logger.info("NCBI API key is set in the environment variables.")


# -----------------------------
# MASTER VALIDATION ENTRYPOINT
# -----------------------------


def main_validator(fasta_files: list[str] | None, chosen: list[stages.Stage]) -> bool:
    """Run the slow checks that `-skp` skips.

    These are NCBI lookups of every probe accession and the API-key prompt (only
    for stages that talk to NCBI), and the genome FASTAs.

    Parameters
    ----------
    fasta_files : list of str or None
        List of FASTA file paths to validate.
    chosen : list of Stage
        The stages about to run.

    Returns:
    -------
    bool
        True if all checks pass, False otherwise.
    """
    logger.debug("Starting input validation process...")
    uses_ncbi = any(stage.flag in NCBI_STAGES for stage in chosen)

    csv_ok = (
        csv_validator(csv_file=defaults.config["input"]["probe_csv"])
        if uses_ncbi
        else True
    )

    if not defaults.USE_SPECIES_DICT and fasta_files:
        fasta_results = [fasta_validator(f) for f in fasta_files]
        fasta_ok = all(fasta_results)
    else:
        fasta_ok = True

    if uses_ncbi:
        validate_ncbi_key()

    return all([csv_ok, fasta_ok])


# -----------------------------
# INTERACTIVE CONFIRMATION
# -----------------------------


def green_light(all_valid: bool) -> bool:
    """Asks user to confirm whether to proceed if all validations passed.

    Parameters
    ----------
    all_valid : bool
        Whether all validations were successful.

    Returns:
    -------
    bool
        True if user wants to proceed, False otherwise.
    """
    if not all_valid:
        logger.error(
            "Validation failed; the reasons are listed above. Nothing was run."
        )
        return False

    logger.info("Validation passed.")
    time.sleep(0.1)

    proceed = ask("Proceed [Y/n]: ", default="Y")
    if proceed.upper() == "Y":
        logger.info(
            "RetroSeek started. Depending on your system, this may take some time."
        )
        return True
    if proceed.upper() == "N":
        logger.warning("Workflow aborted by user.")
        return False
    logger.warning("Invalid input. Please try again.")
    return green_light(all_valid)


# -----------------------------
# ENTRYPOINT
# -----------------------------


def validation_run(
    chosen: list[stages.Stage], fasta_files: list[str] | None = None
) -> bool:
    """The slow, skippable validation, then the confirmation prompt.

    Parameters
    ----------
    chosen : list of Stage
        The stages about to run.
    fasta_files : list of str, optional
        List of FASTA file paths to validate.

    Returns:
    -------
    bool
        True if user confirms execution after passing validation, False otherwise.
    """
    all_valid = main_validator(fasta_files=fasta_files or [], chosen=chosen)

    return green_light(all_valid=all_valid)
