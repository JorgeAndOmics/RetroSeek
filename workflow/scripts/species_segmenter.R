# =============================================================================
# 0. Suppress Warnings and Load Libraries
# =============================================================================
options(warn = -1)
suppressMessages({
  library(yaml)
  library(argparse)
  library(arrow)
  library(tidyverse)
})

# =============================================================================
# 1. Parse Command-Line Arguments
# =============================================================================
parser <- ArgumentParser(description = "Split a Parquet dataset by species")

parser$add_argument("--input_file", required = TRUE, help = "Path to input Parquet file.")
parser$add_argument("--output_dir", required = TRUE, help = "Directory to write output Parquet files.")
parser$add_argument("--config_file", required = TRUE, help = "Path to configuration YAML file.")

args <- parser$parse_args()

args.input_file <- args$input_file
args.output_dir <- args$output_dir
args.config_file <- args$config_file

# =============================================================================
# 2. Load Input File
# =============================================================================
data <- arrow::read_parquet(args.input_file)

# Load configuration file
config <- yaml::read_yaml(args.config_file)

# Load main probes from config
main_probes <- config$parameters$main_probes

# =============================================================================
# 3. Helper Function: Named Group Split
# =============================================================================
group_split_named <- function(df, groupvar) {
  # Tidy evaluation to use function arguments as column names; 
  # if not, R would look for a column called "groupvar"
  grouped_df <- df %>% group_by({{ groupvar }})
  splitted <- grouped_df %>% group_split(.keep = TRUE)
  grp_vals <- grouped_df %>% group_keys() %>% pull({{ groupvar }})
  splitted <- set_names(splitted, grp_vals)
}

# =============================================================================
# 4. Perform Group Splitting
# =============================================================================
all_probes <- group_split_named(data, species)

main_probes <- data %>%
  filter(probe %in% main_probes)

accessory_probes <- data %>%
    filter(!probe %in% main_probes)

species_groups <- group_split_named(data, species)

species_names <- data %>%
  arrange(species) %>%
  distinct(species) %>%
  pull(species)

# =============================================================================
# 5. Write Main & Accessory Parquet Outputs
# =============================================================================
arrow::write_parquet(main_probes, file.path(args.output_dir, "all_main.parquet"))
arrow::write_parquet(accessory_probes, file.path(args.output_dir, "all_accessory.parquet"))

# =============================================================================
# 6. Write Species-Specific Outputs
# =============================================================================
for (sp in species_names) {
  message(paste("Splitting hits for species:", sp))
  
  # ALL
  if (sp %in% names(all_probes)) {
    arrow::write_parquet(all_probes[[sp]], file.path(args.output_dir, paste0(sp, ".parquet")))
  } else {
    arrow::write_parquet(data[0, ], file.path(args.output_dir, paste0(sp, ".parquet")))
  }
}