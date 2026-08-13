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
# `parents[2]` walks: defaults.py -> scripts -> workflow -> repo root.
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
# CONFIG_DIR is REPO-relative, not DATA_DIR-relative: the only two files read
# from it - `schema.yaml` (validator.py) and `erv_class.tsv` (the
# `taxonomy_reference` rule) - are source artifacts versioned with the code,
# not pipeline outputs. Deriving it from DATA_DIR worked by coincidence under
# the committed config (`data_root_folder: 'data'` is already inside the repo)
# but broke whenever the data root pointed elsewhere: the mkdir loop below
# would create an empty `<data_root>/config/` and validation died on a missing
# schema.yaml. Everything else here stays DATA_DIR-relative on purpose.
PATH_DICT["CONFIG_DIR"] = (_REPO_ROOT_FOR_CONFIG / "data" / "config").resolve()
PATH_DICT["SPECIES_DIR"] = (PATH_DICT["DATA_DIR"] / "species").resolve()
# User-provided input tables (e.g. the probe CSV). Kept under a dedicated
# `_input/` subdir so it doesn't sit loose alongside the pipeline's own
# data/tables/<name>/ output subdirectories.
PATH_DICT["TABLE_INPUT_DIR"] = (PATH_DICT["DATA_DIR"] / "tables" / "_input").resolve()
PATH_DICT["PICKLE_DIR"] = (PATH_DICT["DATA_DIR"] / "pickles").resolve()
PATH_DICT["TMP_DIR"] = (PATH_DICT["DATA_DIR"] / "tmp").resolve()
PATH_DICT["TBLASTN_PICKLE_DIR"] = (PATH_DICT["PICKLE_DIR"] / "tblastn").resolve()

# === Taxonomic-classification reference ===
# Pinned, build-once reference for per-locus genus calls: the genus-comprehensive
# protein set (.faa) + accession->genus/gene table (.csv), data-derived taxonomy.tsv,
# curated erv_class.tsv, the per-gene placement tree packages (trees/), and a
# provenance manifest. Lives under /data because it is a reusable input, not a
# per-run output. Rebuilt only when missing (or via `make reference`). The small
# blastx DB is built on the fly in each run's workdir (383 proteins - instant).
PATH_DICT["TAXONOMY_REFERENCE_DIR"] = (
    PATH_DICT["DATA_DIR"] / "taxonomy_reference"
).resolve()
PATH_DICT["TAXONOMY_TREES_DIR"] = (
    PATH_DICT["TAXONOMY_REFERENCE_DIR"] / "trees"
).resolve()

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
# Taxonomic-classification output table - per genome <g>.loci (per-locus genus
# calls; the genus-founded ERV assembly, with per-gene evidence + mosaic packed
# in-row). Plot summaries are derived from this table, so they stay concordant.
(
    PATH_DICT["TAXONOMY_TABLES_PARQUET_DIR"],
    PATH_DICT["TAXONOMY_TABLES_CSV_DIR"],
) = table_dirs("taxonomy_classification")
# Tree coordinates (ADR-011) - flat x/y segment + tip tables written by
# tree_layout.py so the R plot generators can draw the taxon and host-species
# trees with geom_segment, without an R tree library.
PATH_DICT["TAXONOMY_TREE_COORDS_DIR"] = PATH_DICT["TAXONOMY_TABLES_CSV_DIR"] / "trees"
# Per-segment deliverables (ADR-011) - the catalog split by the taxon each locus
# rolls up to at classification.segment_rank.
PATH_DICT["TAXONOMY_SEGMENTS_DIR"] = PATH_DICT["TAXONOMY_TABLES_CSV_DIR"] / "segments"
# Co-phylogeny tables (ADR-014): the ERV-composition tree, the KRD distance
# matrix behind it, and the congruence verdict against the host phylogeny.
PATH_DICT["COPHYLOGENY_DIR"] = PATH_DICT["TAXONOMY_TABLES_CSV_DIR"] / "cophylogeny"
# Loss-analysis tables - the unified per-stage loss funnel + per-genome novel
# candidates (valid loci with zero blastx homology). Built by loss_analysis.R
# from the ranges-analysis counts and the blastx-stage classification counts.
(
    PATH_DICT["LOSS_ANALYSIS_PARQUET_DIR"],
    PATH_DICT["LOSS_ANALYSIS_CSV_DIR"],
) = table_dirs("loss_analysis")
# Run manifest - provenance metadata (generator, timestamp, input md5s,
# resolved parameters, seed), not a table; lives directly under results/.
PATH_DICT["MANIFEST_DIR"] = (PATH_DICT["RESULTS_DIR"] / "manifest").resolve()
# LTRharvest screen-format (.scn) intermediate - consumed by LTR_retriever.
# Lives under /data (not /results) because it's a working format, not an output.
PATH_DICT["LTR_SCN_DIR"] = (PATH_DICT["DATA_DIR"] / "ltr_scn").resolve()
# LTR_RETRIEVER_DIR is defined below, after TRACK_DIR is set up.

# === Results - Plots ===
# Layout mirrors the pipeline stages so the filesystem is self-documenting:
#
#   plots/
#     ranges/                 the ranges_analysis stage (pre-classification)
#       homology/             plot2sort: per-probe/virus/species hit distributions
#       integration/          stage_plot_generator: hit<->LTR-element integration + reduction
#     classification/         the taxonomy stage (per-locus assembly + calls)
#       taxonomy/             taxonomy_plot_generator: calls, confidence, mosaic, tiers
#       structure/            structural views of the assembly (completeness, structure_class)
#       loss/                 loss_analysis: per-stage attrition funnel
#     circle/                 per-genome Circos overviews
#     hotspot/                integration-hotspot enrichment (Manhattan / QQ / karyotype)
PATH_DICT["PLOT_DIR"] = (PATH_DICT["RESULTS_DIR"] / "plots").resolve()

# --- Ranges stage ---
PATH_DICT["RANGES_PLOT_DIR"] = (PATH_DICT["PLOT_DIR"] / "ranges").resolve()
PATH_DICT["HOMOLOGY_PLOT_DIR"] = (PATH_DICT["RANGES_PLOT_DIR"] / "homology").resolve()
PATH_DICT["INTEGRATION_PLOT_DIR"] = (
    PATH_DICT["RANGES_PLOT_DIR"] / "integration"
).resolve()

# --- Classification stage ---
PATH_DICT["CLASSIFICATION_PLOT_DIR"] = (
    PATH_DICT["PLOT_DIR"] / "classification"
).resolve()
PATH_DICT["TAXONOMY_PLOT_DIR"] = (
    PATH_DICT["CLASSIFICATION_PLOT_DIR"] / "taxonomy"
).resolve()
# Phylogenetic-placement figures (ADR-014): per-genome heat-trees plus the
# cross-genome co-phylogeny comparison. Kept beside the taxonomy panels because
# they answer the same question - which lineage is where - from the placement
# evidence rather than from the assembled catalog.
PATH_DICT["PLACEMENT_PLOT_DIR"] = (
    PATH_DICT["CLASSIFICATION_PLOT_DIR"] / "placement"
).resolve()
# gappa writes every artifact of a command into one --out-dir, mixing figures,
# tables and trees. They are split back out by type here so the directory
# contract holds: figures under plots/, CSVs under tables/, and the Newick and
# Nexus trees beside the .jplace they were derived from.
PATH_DICT["PLACEMENT_TABLE_DIR"] = (
    PATH_DICT["TABLE_OUTPUT_DIR"] / "placement"
).resolve()
# Structural views of the genus-founded ERV assembly (completeness, canonical
# order, structure_class). Built by erv_like_plot_generator.R from the loci table;
# the retired erv_like *tier* (ADR-007) is why the panel is now named 'structure'.
PATH_DICT["STRUCTURE_PLOT_DIR"] = (
    PATH_DICT["CLASSIFICATION_PLOT_DIR"] / "structure"
).resolve()
PATH_DICT["LOSS_PLOT_DIR"] = (PATH_DICT["CLASSIFICATION_PLOT_DIR"] / "loss").resolve()

# --- Standalone analyses ---
PATH_DICT["CIRCLE_PLOT_DIR"] = (PATH_DICT["PLOT_DIR"] / "circle").resolve()
PATH_DICT["HOTSPOT_PLOT_DIR"] = (PATH_DICT["PLOT_DIR"] / "hotspot").resolve()

# === Results - Tracks ===
PATH_DICT["TRACK_DIR"] = (PATH_DICT["RESULTS_DIR"] / "tracks").resolve()
PATH_DICT["TRACK_ORIGINAL_DIR"] = (PATH_DICT["TRACK_DIR"] / "original").resolve()
PATH_DICT["TRACK_CANDIDATES_DIR"] = (PATH_DICT["TRACK_DIR"] / "candidates").resolve()
PATH_DICT["TRACK_VALID_DIR"] = (PATH_DICT["TRACK_DIR"] / "valid").resolve()
PATH_DICT["TRACK_HOTSPOTS_DIR"] = (PATH_DICT["TRACK_DIR"] / "hotspots").resolve()
# Taxonomic-classification tier - per-locus genus calls projected to genome
# coordinates (GFF3 + BED for IGV, colour-by-genus). Additive to the valid tier.
PATH_DICT["TRACK_TAXONOMY_DIR"] = (PATH_DICT["TRACK_DIR"] / "taxonomy").resolve()
# Published phylogenetic-placement evidence: one .jplace + labelled .newick per
# (genome, tier, placement gene). EPA-ng writes these into TMP_DIR, which is
# cleared between runs; they are the evidence behind every taxon_call and the
# interchange format iTOL/gappa read, so they are promoted here (ADR-014).
PATH_DICT["PLACEMENT_DIR"] = (PATH_DICT["TRACK_TAXONOMY_DIR"] / "placements").resolve()
# Orphans tier - non-LTR-associated hits recovered + classified by their own
# sequence (parallel to taxonomy; only orphans that earn a taxonomic call).
PATH_DICT["TRACK_ORPHANS_DIR"] = (PATH_DICT["TRACK_DIR"] / "orphans").resolve()

# === Results - LTR ===
PATH_DICT["LTRHARVEST_DIR"] = (PATH_DICT["TRACK_DIR"] / "ltrharvest").resolve()
PATH_DICT["LTRDIGEST_DIR"] = (PATH_DICT["TRACK_DIR"] / "ltrdigest").resolve()
# LTR_retriever output directory (intact-ERV filtered list, solo-LTR list,
# consensus library - all the files LTR_retriever emits per genome).
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
    # to .fa - otherwise a fresh machine with only .fna files would see
    # SPECIES = [] before the normalizer ever runs.
    _FASTA_EXTS = {".fa", ".fna", ".fasta", ".ffn"}
    SPECIES = sorted(
        {f.stem for f in PATH_DICT["SPECIES_DB"].iterdir() if f.suffix in _FASTA_EXTS}
    )
else:
    SPECIES = list(SPECIES_DICT.keys())
