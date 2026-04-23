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
# 1. Parse Command-Line Arguments
# =============================================================================
parser <- ArgumentParser(description = "Split a Parquet dataset by species")
parser$add_argument("--input_file", required = TRUE, help = "Path to input Parquet file.")
parser$add_argument("--output_dir", required = TRUE, help = "Directory to write output Parquet files.")
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

# =============================================================================
# 2. Load Input File and Config
# =============================================================================
data <- arrow::read_parquet(args$input_file)
config <- yaml::read_yaml(args$config_file)

# Store probe names from config
main_probe_names <- config$parameters$main_probes

# =============================================================================
# 3. Write Main & Accessory Parquet Outputs (Filter Once)
# =============================================================================
main_probes_df <- data %>% filter(probe %in% main_probe_names)
accessory_probes_df <- data %>% filter(!probe %in% main_probe_names)

message("Writing full main and accessory parquet files.")
arrow::write_parquet(main_probes_df, file.path(args$output_dir, "all_main.parquet"))
arrow::write_parquet(accessory_probes_df, file.path(args$output_dir, "all_accessory.parquet"))

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
  arrow::write_parquet(
    species_list[[i]],
    file.path(args$output_dir, paste0(sp, ".parquet"))
  )
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
all_species <- strsplit(args$species, ",", fixed = TRUE)[[1]]
all_species <- trimws(all_species)
all_species <- all_species[nzchar(all_species)]
missing_species <- setdiff(all_species, written)
if (length(missing_species) > 0L) {
  # Use the same schema as the input so readers don't need to branch on
  # missing-column edge cases. Zero-row slice preserves all column types.
  empty_df <- data[0L, , drop = FALSE]
  for (sp in missing_species) {
    message(paste("Writing empty parquet for species (no hits):", sp))
    arrow::write_parquet(
      empty_df,
      file.path(args$output_dir, paste0(sp, ".parquet"))
    )
  }
}