# Turn off warnings for a cleaner console (optional)
options(warn = -1)

# Quietly load the necessary packages
suppressMessages({
  library(arrow)
  library(tidyverse)
})

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Please provide exactly two arguments: input file and output directory.")
}

# Input/output paths
input_file <- args[1]
output_dir <- args[2]

# Read the Parquet file
data <- arrow::read_parquet(input_file)

# Define a helper function to split a data frame by a given grouping variable
# while automatically setting the resulting list's names to match each group.
group_split_named <- function(df, groupvar) {
  # 1) Group by the specified column
  grouped_df <- df %>% group_by({{ groupvar }})
  
  # 2) Split into a list of data frames, retaining the group var column
  splitted <- grouped_df %>% group_split(.keep = TRUE)
  
  # 3) Extract the distinct group values
  grp_vals <- grouped_df %>% group_keys() %>% pull({{ groupvar }})
  
  # 4) Assign these group values as names of the list elements
  splitted <- set_names(splitted, grp_vals)
}

# Split the full dataset by species (named properly)
all_probes <- group_split_named(data, species)

# Create main and accessory subsets
main_probes <- data %>% filter(probe %in% c("ENV", "POL", "GAG"))
accessory_probes <- data %>% filter(!probe %in% c("ENV", "POL", "GAG"))

# Split main and accessory subsets by species, also with proper names
species_groups_main <- group_split_named(main_probes, species)
species_groups_accessory <- group_split_named(accessory_probes, species)

# Collect a sorted vector of all species seen in the data
species_names <- data %>%
  arrange(species) %>%
  distinct(species) %>%
  pull(species)

# Write single Parquet files for all main and all accessory
message("Splitting hits for main and accessory probes")
arrow::write_parquet(main_probes, file.path(output_dir, "all_main.parquet"))
arrow::write_parquet(accessory_probes, file.path(output_dir, "all_accessory.parquet"))

# For each species, write out full, main, accessory
for (sp in species_names) {
  message(paste("Splitting hits for species:", sp))
  
  # ALL probes for sp
  if (sp %in% names(all_probes)) {
    arrow::write_parquet(
      all_probes[[sp]],
      file.path(output_dir, paste0(sp, "_full.parquet"))
    )
  } else {
    # Write an empty file if the species doesn't exist in all_probes
    arrow::write_parquet(data[0, ], file.path(output_dir, paste0(sp, "_full.parquet")))
  }
  
  # MAIN probes for sp
  if (sp %in% names(species_groups_main)) {
    arrow::write_parquet(
      species_groups_main[[sp]],
      file.path(output_dir, paste0(sp, "_main.parquet"))
    )
  } else {
    # Write an empty file if there are no main probes for sp
    arrow::write_parquet(main_probes[0, ], file.path(output_dir, paste0(sp, "_main.parquet")))
  }
  
  # ACCESSORY probes for sp
  if (sp %in% names(species_groups_accessory)) {
    arrow::write_parquet(
      species_groups_accessory[[sp]],
      file.path(output_dir, paste0(sp, "_accessory.parquet"))
    )
  } else {
    # Write an empty file if there are no accessory probes for sp
    arrow::write_parquet(accessory_probes[0, ], file.path(output_dir, paste0(sp, "_accessory.parquet")))
  }
}
