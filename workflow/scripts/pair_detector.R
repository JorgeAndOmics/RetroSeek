# -------------------
# DEPENDENCIES
# -------------------

# Suppress warnings so that the output is cleaner; warnings won't appear on-screen.
options(warn = -1)

# Suppress messages from loading libraries, so the console remains uncluttered.
suppressMessages({
  library(yaml)           # For reading YAML configuration file
  library(arrow)          # Provides tools for reading and writing Parquet files
  library(tidyverse)      # Collection of R packages for data manipulation and visualization
  library(waffle)         # For creating waffle charts
  library(argparse)       # Command-line argument parsing
  library(GenomicRanges)  # Genomic interval operations
  library(plyranges)      # "Tidyverse"-style GRanges operations
  library(rtracklayer)    # Reading/writing genome annotation files (e.g., GFF, BED)
})

# -------------------------------
# 1b. SOURCE PURE PAIRING HELPERS (get_5prime + find_pairs; unit-tested)
# -------------------------------
.resolve_script_dir <- function() {
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile)) return(dirname(ofile))
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- cmd_args[grepl("^--file=", cmd_args)]
  if (length(file_arg) > 0L) return(dirname(sub("^--file=", "", file_arg[1])))
  "."
}
.script_dir <- .resolve_script_dir()
source(file.path(.script_dir, "pair_detector", "pairing.R"))

# -------------------------------
# 2. PARSE COMMAND-LINE ARGUMENTS
# -------------------------------
parser <- ArgumentParser(description = 'Process tBLASTn and LTRdigest integration overlaps')

# Define expected command-line arguments
parser$add_argument("--config", required=TRUE, help="Configuration YAML file")
parser$add_argument("--gff", required=TRUE, help="GFF3 input: Ranges file to analyse")
parser$add_argument("--parquet_dir", required=TRUE, help="Directory for pipeline-internal Parquet output")
parser$add_argument("--csv_dir", required=TRUE, help="Directory for user-facing CSV output")

args <- parser$parse_args()

# Print which tBLASTn file is being processed
file.basename <- tools::file_path_sans_ext(basename(args$gff))
print(paste0("Processing pairs for ", file.basename, "..."))

# -----------------------------
# 3. CONFIGURATION PARAMETERS
# -----------------------------
# Read YAML configuration file
config <- yaml::read_yaml(args$config)

# Import thresholds and merging behavior from config
probe_to_pair <- config$parameters$probe_to_pair

# -----------------------------
# 3. GFF3 FILE IMPORT
# -----------------------------
# Import the GFF3 file containing validated hits
ranges <- rtracklayer::import(args$gff, format="gff3")

# Assert the existence of the required fields
required_fields <- c("probe", "label", "virus")
missing_fields <- setdiff(required_fields, names(mcols(ranges)))

if (length(missing_fields) > 0) {
  stop(paste("The following required fields are missing from the GFF3 file:", paste(missing_fields, collapse=", ")))
}

# Assert the existence of the required probe from config within the `probe` field
if (!any(grepl(probe_to_pair, ranges$probe))) {
  stop(paste("The probe", probe_to_pair, "is not present in the GFF3 file."))
}

# -----------------------------
# 4. DATA PREPARATION
# -----------------------------
# Remove duplicates from ranges
ranges <- ranges[!duplicated(ranges)]

# Maximum gap (bp) allowed between paired ranges, from config
max_gap <- as.integer(config$parameters$pair_max_gap)

# Pairing logic (get_5prime + find_pairs) lives in pair_detector/pairing.R
pair_df <- find_pairs(ranges, probe_to_pair, max_gap)


# ------------------------------
# 5. PLOT GENERATION
# ------------------------------
# Create a waffle plot to visualize the distribution of pairs
pair_counts <- pair_df %>%
  group_by(other.probe) %>%
  summarise(count = dplyr::n(), .groups = 'drop') %>%
  mutate(type = ifelse(toupper(other.probe) %in% c("ENV", "GAG", "POL"), "main", "accessory")) %>%
  arrange(desc(count))

# Create a waffle plot
waffle_plot <- pair_counts %>%
  ggplot(aes(fill = other.probe, values = count)) +
  geom_waffle(size = 0.1, color = "white") +
  labs(title = paste("Distribution of pairs for", config$parameters$probe_to_pair),
       subtitle = paste("Probe:", probe_to_pair)) +
  theme_minimal() +
  facet_grid(~type) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),  # removes x-axis numbers
        axis.ticks.x = element_blank(),  # removes x-axis ticks
        axis.ticks.y = element_blank(),  # removes y-axis ticks
        axis.text.y = element_blank(),  # removes y-axis numbers
  )


# ------------------------------
# EXPORT RESULTS
# ------------------------------
write.csv(pair_df, file.path(args$csv_dir, paste0(file.basename, ".csv")), row.names = FALSE)
arrow::write_parquet(pair_df, file.path(args$parquet_dir, paste0(file.basename, ".parquet")))
