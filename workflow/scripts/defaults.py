"""
Defaults Configuration Script
=============================

This script loads configuration settings from a YAML file and sets up various
constants and directory paths used throughout the project.

Modules:
    - `yaml`: For parsing the YAML configuration file.
    - `pathlib.Path`: For handling file and directory paths.

Configuration:
    - The configuration file is expected to be located at `../data/config/config.yaml`.
    - Various constants and directory paths are initialized based on the configuration file.

Usage:
    This script is intended to be imported as a module and not run directly.
"""

import os
from pathlib import Path

import yaml

# Config-file resolution:
#   1. ``RETROSEEK_CONFIG`` env var, if set (absolute path or relative to repo root).
#   2. ``data/config/config.yaml`` (the committed default; portable repo-relative paths).
#
# Why an env var rather than a CLI flag: ``defaults.py`` is imported at
# Snakemake parse time, before any rule body or shell expansion runs.
# Snakemake's ``--configfile`` only affects ``config[...]`` lookups inside
# the Snakefile; it cannot override which file ``defaults.py`` reads.
# The env var lets ``./RetroSeek --configfile`` (or a manual snakemake
# invocation) declare the config before the import side effects fire.
_REPO_ROOT_FOR_CONFIG = Path(__file__).resolve().parents[2]
_env_config = os.environ.get("RETROSEEK_CONFIG")
if _env_config:
    _candidate = Path(_env_config)
    CONFIG_FILE = (
        _candidate if _candidate.is_absolute() else _REPO_ROOT_FOR_CONFIG / _candidate
    )
else:
    CONFIG_FILE = _REPO_ROOT_FOR_CONFIG / "data" / "config" / "config.yaml"

with CONFIG_FILE.open() as f:
    config = yaml.safe_load(f)

# BLAST
E_VALUE = config["blast"]["e_value"]
ACCESSION_ID_REGEX = r"[A-Z]{2,}_?[0-9]+\.[0-9]{1,2}"
PROBE_MIN_LENGTH = config["parameters"]["probe_min_length"]

# Logging
LEVEL_STYLES = config["logging"]["level_styles"]
FIELD_STYLES = config["logging"]["field_styles"]

# Anchor relative paths against the repo root so a fresh-clone run
# from any working directory still resolves the same way.
# `parents[2]` walks: defaults.py → scripts → workflow → repo root.
_REPO_ROOT = Path(__file__).resolve().parents[2]


def _anchor(value: str | Path, fallback: str | Path) -> Path:
    """Resolve a config path: absolute paths win, relative anchor at repo root."""
    raw = Path(value if value is not None else fallback)
    if raw.is_absolute():
        return raw.resolve()
    return (_REPO_ROOT / raw).resolve()


# Directories
PATH_DICT = {"ROOT": _anchor(config["root"].get("db_root_folder"), "data/species")}

# === Root Directories ===
PATH_DICT["DATA_DIR"] = _anchor(config["root"].get("data_root_folder"), "data")
PATH_DICT["RESULTS_DIR"] = _anchor(config["root"].get("results_root_folder"), "results")
PATH_DICT["LOG_DIR"] = _anchor(config["root"].get("logs_root_folder"), "logs")
PATH_DICT["WORKFLOW_DIR"] = Path(__file__).parent

# === Workflow Directories ===
PATH_DICT["SCRIPTS_DIR"] = (PATH_DICT["WORKFLOW_DIR"] / "scripts").resolve()

# === Database Directories ===
PATH_DICT["ROOT_DB"] = (PATH_DICT["ROOT"]).resolve()
PATH_DICT["SPECIES_DB"] = (PATH_DICT["ROOT_DB"]).resolve()
PATH_DICT["ACCESSORY_DB"] = (PATH_DICT["ROOT"] / "accessory").resolve()

# === Data Subdirectories ===
PATH_DICT["CONFIG_DIR"] = (PATH_DICT["DATA_DIR"] / "config").resolve()
PATH_DICT["SPECIES_DIR"] = (PATH_DICT["DATA_DIR"] / "species").resolve()
# User-provided input tables (e.g. the probe CSV). Kept under a dedicated
# `_input/` subdir so it doesn't sit loose alongside the pipeline's own
# data/tables/<name>/ output subdirectories.
PATH_DICT["TABLE_INPUT_DIR"] = (PATH_DICT["DATA_DIR"] / "tables" / "_input").resolve()
PATH_DICT["PICKLE_DIR"] = (PATH_DICT["DATA_DIR"] / "pickles").resolve()
PATH_DICT["TMP_DIR"] = (PATH_DICT["DATA_DIR"] / "tmp").resolve()
PATH_DICT["TBLASTN_PICKLE_DIR"] = (PATH_DICT["PICKLE_DIR"] / "tblastn").resolve()

# === Results - Tables ===
PATH_DICT["TABLE_OUTPUT_DIR"] = (PATH_DICT["RESULTS_DIR"] / "tables").resolve()


def table_dirs(name: str) -> tuple[Path, Path]:
    """Resolve the (parquet, csv) directory pair for a named table group.

    Pipeline-internal parquet copies live under ``data/tables/<name>/``;
    user-facing CSV copies under ``results/tables/<name>/``. The caller is
    expected to store both in PATH_DICT so they are auto-created below.
    """
    return (
        (PATH_DICT["DATA_DIR"] / "tables" / name).resolve(),
        (PATH_DICT["TABLE_OUTPUT_DIR"] / name).resolve(),
    )


# Every per-table directory is a (parquet, csv) pair via table_dirs(): the
# pipeline-internal Parquet copies under data/tables/<name>/, the user-facing
# CSV copies under results/tables/<name>/. No table file sits loose in tables/.
(
    PATH_DICT["HOTSPOT_PARQUET_DIR"],
    PATH_DICT["HOTSPOT_CSV_DIR"],
) = table_dirs("hotspots")
(
    PATH_DICT["OVERLAP_MATRIX_PARQUET_DIR"],
    PATH_DICT["OVERLAP_MATRIX_CSV_DIR"],
) = table_dirs("overlap_matrix")
(
    PATH_DICT["SEGMENTED_SPECIES_PARQUET_DIR"],
    PATH_DICT["SEGMENTED_SPECIES_CSV_DIR"],
) = table_dirs("segmented_species")
# Per-genome ranges-analysis tables (final_loci / homology_loci / ltr_structure
# / reduction_multiplicity / counts).
(
    PATH_DICT["RANGES_ANALYSIS_TABLES_PARQUET_DIR"],
    PATH_DICT["RANGES_ANALYSIS_TABLES_CSV_DIR"],
) = table_dirs("ranges_analysis")
(
    PATH_DICT["PROBE_PAIRS_PARQUET_DIR"],
    PATH_DICT["PROBE_PAIRS_CSV_DIR"],
) = table_dirs("probe_pairs")
(
    PATH_DICT["SOLO_INTACT_PARQUET_DIR"],
    PATH_DICT["SOLO_INTACT_CSV_DIR"],
) = table_dirs("solo_intact_ratio")
(
    PATH_DICT["PROBE_DICT_PARQUET_DIR"],
    PATH_DICT["PROBE_DICT_CSV_DIR"],
) = table_dirs("probe_dict")
(
    PATH_DICT["FULL_GENOME_BLAST_PARQUET_DIR"],
    PATH_DICT["FULL_GENOME_BLAST_CSV_DIR"],
) = table_dirs("full_genome_blast")
# Run manifest — provenance metadata (generator, timestamp, input md5s,
# resolved parameters, seed), not a table; lives directly under results/.
PATH_DICT["MANIFEST_DIR"] = (PATH_DICT["RESULTS_DIR"] / "manifest").resolve()
# LTRharvest screen-format (.scn) intermediate — consumed by LTR_retriever.
# Lives under /data (not /results) because it's a working format, not an output.
PATH_DICT["LTR_SCN_DIR"] = (PATH_DICT["DATA_DIR"] / "ltr_scn").resolve()
# LTR_RETRIEVER_DIR is defined below, after TRACK_DIR is set up.

# === Results - Plots ===
PATH_DICT["PLOT_DIR"] = (PATH_DICT["RESULTS_DIR"] / "plots").resolve()
# Provirus panel: the per-probe / stage plots about individual proviral hits and
# their LTR integration (plot2sort + stage_plot_generator) live under provirus/.
PATH_DICT["PROVIRUS_PLOT_DIR"] = (PATH_DICT["PLOT_DIR"] / "provirus").resolve()
# ERV-like panel: plots about assembled erv_like candidates (composition,
# completeness, gene order). Hyphenated per the output contract.
PATH_DICT["ERV_LIKE_PLOT_DIR"] = (PATH_DICT["PLOT_DIR"] / "erv-like").resolve()
PATH_DICT["CIRCLE_PLOT_DIR"] = (PATH_DICT["PLOT_DIR"] / "circle_plots").resolve()
PATH_DICT["HOTSPOT_PDF_DIR"] = (PATH_DICT["PLOT_DIR"] / "hotspot_pdfs").resolve()

# === Results - Tracks ===
PATH_DICT["TRACK_DIR"] = (PATH_DICT["RESULTS_DIR"] / "tracks").resolve()
PATH_DICT["TRACK_ORIGINAL_DIR"] = (PATH_DICT["TRACK_DIR"] / "original").resolve()
PATH_DICT["TRACK_CANDIDATES_DIR"] = (PATH_DICT["TRACK_DIR"] / "candidates").resolve()
PATH_DICT["TRACK_VALID_DIR"] = (PATH_DICT["TRACK_DIR"] / "valid").resolve()
# ERV-like assembly tier — composite candidates chained from valid main-probe
# loci (see assemble_erv_like in erv_assembly.R). Additive to the valid tier.
PATH_DICT["TRACK_ERV_LIKE_DIR"] = (PATH_DICT["TRACK_DIR"] / "erv_like").resolve()
PATH_DICT["TRACK_HOTSPOTS_DIR"] = (PATH_DICT["TRACK_DIR"] / "hotspots").resolve()

# === Results - LTR ===
PATH_DICT["LTRHARVEST_DIR"] = (PATH_DICT["TRACK_DIR"] / "ltrharvest").resolve()
PATH_DICT["LTRDIGEST_DIR"] = (PATH_DICT["TRACK_DIR"] / "ltrdigest").resolve()
# LTR_retriever output directory (intact-ERV filtered list, solo-LTR list,
# consensus library — all the files LTR_retriever emits per genome).
PATH_DICT["LTR_RETRIEVER_DIR"] = (PATH_DICT["TRACK_DIR"] / "ltr_retriever").resolve()
PATH_DICT["SOLO_LTR_DIR"] = (PATH_DICT["TRACK_DIR"] / "solo_ltr").resolve()
PATH_DICT["FLANKING_LTR_DIR"] = (PATH_DICT["TRACK_DIR"] / "flanking_ltr").resolve()

# === Logs & Workflow ===
PATH_DICT["DOWNLOAD_LOG"] = (PATH_DICT["LOG_DIR"] / "download_log.log").resolve()

# === Accessory Tools ===
PATH_DICT["HMM_PROFILE_DIR"] = (PATH_DICT["ACCESSORY_DB"] / "hmm_profiles").resolve()


# Directory generation
for value in PATH_DICT.values():
    value.mkdir(parents=True, exist_ok=True)

# Execution and requests
NUM_CORES = config["execution"].get("num_cores", 1)
USE_SPECIES_DICT = config.get("execution", False).get("use_species_dict", False)
RETRIVAL_TIME_LAG = config["execution"].get("retrieval_time_lag", 0.3)
MAX_RETRIEVAL_ATTEMPTS = config["execution"].get("max_retrieval_attempts", 3)
MAX_THREADPOOL_WORKERS = config["execution"].get("max_threadpool_workers", 1)
ENTREZ_EMAIL = config["execution"].get("entrez_email", "")

# Display
DISPLAY_SNAKEMAKE_INFO: bool = config["display"].get("display_snakemake_info", False)
DISPLAY_REQUESTS_WARNING: bool = config["display"].get(
    "display_requests_warning", False
)
DISPLAY_OPERATION_INFO: bool = config["display"].get("display_operation_info", False)

# INPUT
PROBE_CSV = _anchor(config["input"].get("probe_csv"), "data/tables/_input/probes.csv")

# Genomes
SPECIES_DICT: dict[str, str] = config.get("species", {})

SPECIES: list[str]
if not USE_SPECIES_DICT:
    # Discover genomes by scanning SPECIES_DB. Accept any of the FASTA
    # extension variants the genome_fasta_normalizer rule canonicalises
    # to .fa — otherwise a fresh machine with only .fna files would see
    # SPECIES = [] before the normalizer ever runs.
    _FASTA_EXTS = {".fa", ".fna", ".fasta", ".ffn"}
    SPECIES = sorted(
        {f.stem for f in PATH_DICT["SPECIES_DB"].iterdir() if f.suffix in _FASTA_EXTS}
    )
else:
    SPECIES = list(SPECIES_DICT.keys())
