# Installing libraries
library(arrow)
library(tidyverse)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Please provide exactly two arguments: input file and output directory.")
}

input_file <- args[1]
output_dir <- args[2]

# Loading data
data <- arrow::read_parquet(input_file)

# Ensure output directory exists
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Segment by species
all_probes <- data %>%
  group_by(species) %>%
  group_split()

# Segment main and accessory probes
main_probes <- data %>%
  filter(probe %in% c("ENV", "POL", "GAG"))

accessory_probes <- data %>%
  filter(!probe %in% c("ENV", "POL", "GAG"))

# Segmenting species in main and accessory probes
species_groups_main <- main_probes %>%
  group_by(species) %>%
  group_split()

species_groups_accessory <- accessory_probes %>%
  group_by(species) %>%
  group_split()

# Assign names to each split dataframe based on the species
species_names <- data %>%
  arrange(species) %>%
  distinct(species) %>%
  pull(species)

names(all_probes) <- species_names
names(species_groups_main) <- species_names
names(species_groups_accessory) <- species_names

# Create Parquet files for main and accessory probes
print("Processing hits for main and accessory probes")
arrow::write_parquet(main_probes,  file.path(output_dir, "all_main.parquet"))
arrow::write_parquet(accessory_probes, file.path(output_dir, "all_accessory.parquet"))

# Create Parquet files for all categories
for (species_name in species_names) {
  print(paste("Processing hits for ", species_name))
  # Save all probes
  if (!is.null(all_probes[[species_name]])) {
    output_file_all <- file.path(output_dir, paste0(species_name, "_full", ".parquet"))
    arrow::write_parquet(all_probes[[species_name]], output_file_all)
  }
  
  # Save main probes
  if (!is.null(species_groups_main[[species_name]])) {
    output_file_main <- file.path(output_dir, paste0(species_name, "_main", ".parquet"))
    arrow::write_parquet(species_groups_main[[species_name]], output_file_main)
  }
  
  # Save accessory probes
  if (!is.null(species_groups_accessory[[species_name]])) {
    output_file_accessory <- file.path(output_dir, paste0(species_name, "_accessory", ".parquet"))
    arrow::write_parquet(species_groups_accessory[[species_name]], output_file_accessory)
  }
}
