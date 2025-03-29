# =============================================================================
# 0. Suppress Warnings and Load Libraries
# =============================================================================
options(warn = -1)
suppressMessages({
  library(argparse)
  library(arrow)
  library(tidyverse)
})

# =============================================================================
# 1. Parse Command-Line Arguments
# =============================================================================
parser <- ArgumentParser(description = "Split a Parquet dataset by species and probe class.")

parser$add_argument("--input_file", required = TRUE, help = "Path to input Parquet file.")
parser$add_argument("--output_dir", required = TRUE, help = "Directory to write output Parquet files.")

args <- parser$parse_args()

input_file <- args$input_file
output_dir <- args$output_dir

# =============================================================================
# 2. Load Input File
# =============================================================================
data <- arrow::read_parquet(input_file)

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

main_probes <- data %>% filter(toupper(probe) %in% c("ENV", "POL", "GAG"))
accessory_probes <- data %>% filter(!toupper(probe) %in% c("ENV", "POL", "GAG"))

species_groups_main <- group_split_named(main_probes, species)
species_groups_accessory <- group_split_named(accessory_probes, species)

species_names <- data %>%
  arrange(species) %>%
  distinct(species) %>%
  pull(species)

# =============================================================================
# 5. Write Main & Accessory Parquet Outputs
# =============================================================================
message("Splitting hits for main and accessory probes")
arrow::write_parquet(main_probes, file.path(output_dir, "all_main.parquet"))
arrow::write_parquet(accessory_probes, file.path(output_dir, "all_accessory.parquet"))

# =============================================================================
# 6. Write Species-Specific Outputs
# =============================================================================
for (sp in species_names) {
  message(paste("Splitting hits for species:", sp))
  
  # ALL
  if (sp %in% names(all_probes)) {
    arrow::write_parquet(all_probes[[sp]], file.path(output_dir, paste0(sp, "_full.parquet")))
  } else {
    arrow::write_parquet(data[0, ], file.path(output_dir, paste0(sp, "_full.parquet")))
  }
  
  # MAIN
  if (sp %in% names(species_groups_main)) {
    arrow::write_parquet(species_groups_main[[sp]], file.path(output_dir, paste0(sp, "_main.parquet")))
  } else {
    arrow::write_parquet(main_probes[0, ], file.path(output_dir, paste0(sp, "_main.parquet")))
  }
  
  # ACCESSORY
  if (sp %in% names(species_groups_accessory)) {
    arrow::write_parquet(species_groups_accessory[[sp]], file.path(output_dir, paste0(sp, "_accessory.parquet")))
  } else {
    arrow::write_parquet(accessory_probes[0, ], file.path(output_dir, paste0(sp, "_accessory.parquet")))
  }
}
