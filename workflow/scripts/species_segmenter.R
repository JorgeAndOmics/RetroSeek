# =============================================================================
# 0. Suppress Warnings and Load Libraries
# =============================================================================
options(warn = -1)
suppressMessages({
  library(yaml)
  library(argparse)
  library(arrow)
  library(dplyr)  # More explicit than tidyverse
})

# =============================================================================
# 1. Parse Command-Line Arguments
# =============================================================================
parser <- ArgumentParser(description = "Split a Parquet dataset by species")
parser$add_argument("--input_file", required = TRUE, help = "Path to input Parquet file.")
parser$add_argument("--output_dir", required = TRUE, help = "Directory to write output Parquet files.")
parser$add_argument("--config_file", required = TRUE, help = "Path to configuration YAML file.")

args <- parser$parse_args()

# =============================================================================
# 2. Load Input File and Config
# =============================================================================
data <- arrow::read_parquet(args$input_file)
config <- yaml::read_yaml(args$config_file)

# Store probe names from config (don't overwrite later)
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
species_list <- split(data, data$species)

species_names <- data %>%
  distinct(species) %>%
  pull(species)

names(species_list) <- species_names

# Write species-specific files
for (sp in species_names) {
  message(paste("Writing hits for species:", sp))
  arrow::write_parquet(
    species_list[[sp]],
    file.path(args$output_dir, paste0(sp, ".parquet"))
  )
}