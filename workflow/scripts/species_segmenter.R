# =============================================================================
# 0. Suppress Warnings and Load Libraries
# =============================================================================
options(warn = -1)
suppressMessages({
  library(yaml)
  library(argparse)
  library(arrow)
  library(dplyr)
})

# =============================================================================
# 0b. Source pure helpers (partition + species-list logic; unit-tested)
# =============================================================================
.resolve_script_dir <- function() {
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile)) return(dirname(ofile))
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- cmd_args[grepl("^--file=", cmd_args)]
  if (length(file_arg) > 0L) return(dirname(sub("^--file=", "", file_arg[1])))
  "."
}
.script_dir <- .resolve_script_dir()
source(file.path(.script_dir, "species_segmenter", "segment.R"))

# =============================================================================
# 1. Parse Command-Line Arguments
# =============================================================================
parser <- ArgumentParser(description = "Split a Parquet dataset by species")
parser$add_argument("--input_file", required = TRUE, help = "Path to input Parquet file.")
parser$add_argument("--parquet_dir", required = TRUE, help = "Directory for pipeline-internal Parquet outputs.")
parser$add_argument("--csv_dir", required = TRUE, help = "Directory for user-facing CSV outputs.")
parser$add_argument("--config_file", required = TRUE, help = "Path to configuration YAML file.")
parser$add_argument(
  "--species",
  required = TRUE,
  help = paste(
    "Comma-separated list of ALL expected species IDs (matching Snakefile's",
    "SPECIES list). The script writes an empty parquet for any listed species",
    "not present in the input data, so Snakemake's declared outputs for the",
    "species_segmenter_setup rule always materialise. Needed for the B1",
    "checkpoint-based DAG to resolve downstream dependencies correctly."
  )
)

args <- parser$parse_args()

# Write a data frame as both Parquet (pipeline-internal, data/tables/) and CSV
# (user-facing, results/tables/), under the matching per-table subdirectory.
write_both <- function(df, name) {
  arrow::write_parquet(df, file.path(args$parquet_dir, paste0(name, ".parquet")))
  utils::write.csv(df, file.path(args$csv_dir, paste0(name, ".csv")), row.names = FALSE)
}

# =============================================================================
# 2. Load Input File and Config
# =============================================================================
data <- arrow::read_parquet(args$input_file)
config <- yaml::read_yaml(args$config_file)

# Store probe names from config
main_probe_names <- config$parameters$main_probes

# =============================================================================
# 3. Write Main & Accessory Parquet Outputs (partition once)
# =============================================================================
segments <- segment_by_probe(data, main_probe_names)

message("Writing full main and accessory tables.")
write_both(segments$main, "all_main")
write_both(segments$accessory, "all_accessory")

# =============================================================================
# 4. Split by Species and Write Per-Species Files
# =============================================================================
species_list <- data %>%
  group_by(species) %>%
  group_split(.keep = TRUE)

# Write species-specific files by extracting name directly from each group
written <- character(0)
for (i in seq_along(species_list)) {
  sp <- unique(species_list[[i]]$species)
  message(paste("Writing hits for species:", sp))
  write_both(species_list[[i]], sp)
  written <- c(written, sp)
}

# =============================================================================
# 5. Ensure every configured species has a parquet file (empty if no hits)
# =============================================================================
# The B1 checkpoint-based DAG expects species_segmenter_setup to produce a
# {genome}.parquet for every SPECIES passed from the Snakefile — not just
# the species that happened to produce BLAST hits. Downstream aggregates
# consume only the subset with hits (via species_with_hits(wildcards) at
# runtime) but the per-genome rule outputs must exist for Snakemake to
# resolve the dependency graph.
all_species <- parse_species_list(args$species)
missing_species <- species_to_backfill(all_species, written)
if (length(missing_species) > 0L) {
  # Use the same schema as the input so readers don't need to branch on
  # missing-column edge cases. Zero-row slice preserves all column types.
  empty_df <- data[0L, , drop = FALSE]
  for (sp in missing_species) {
    message(paste("Writing empty table for species (no hits):", sp))
    write_both(empty_df, sp)
  }
}