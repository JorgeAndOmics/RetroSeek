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
# 2. PARSE COMMAND-LINE ARGUMENTS
# -------------------------------
parser <- ArgumentParser(description = 'Process tBLASTn and LTRdigest integration overlaps')

# Define expected command-line arguments
parser$add_argument("--config", required=TRUE, help="Configuration YAML file")
parser$add_argument("--gff", required=TRUE, help="GFF3 input: Ranges file to analyse")
parser$add_argument("--table_output_dir", required=TRUE, help="Output directory for table results")

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
# function to find all "self vs. other" pairs within max_gap base pairs
max_gap <- as.integer(config$parameters$pair_max_gap)

find_pairs <- function(ranges, max_gap) {
  
  # Find overlaps between the "probe" and other labels
  self_ranges <- ranges[ranges$probe == probe_to_pair]
  other_ranges <- ranges[ranges$probe != probe_to_pair]
  
  # 2. Find all overlaps (or up-to max_gap apart)
  ov <- GenomicRanges::findOverlaps(self_ranges, other_ranges, maxgap = max_gap)
  
  # 3. Materialize a data.frame/tibble of exact coordinates
  #    Note: seqnames() and strand() return Rle objects, so coerce to character.
  pairs_df <- tibble(
    self.seqname   = as.character(seqnames(self_ranges)[queryHits(ov)]),
    other.seqname  = as.character(seqnames(other_ranges)[subjectHits(ov)]),
    
    self.start     = start(self_ranges)[queryHits(ov)],
    other.start    = start(other_ranges)[subjectHits(ov)],
    
    self.end       = end(self_ranges)[queryHits(ov)],
    other.end      = end(other_ranges)[subjectHits(ov)],
    
    self.strand    = as.character(strand(self_ranges)[queryHits(ov)]),
    other.strand   = as.character(strand(other_ranges)[subjectHits(ov)]),
    
    self.label     = mcols(self_ranges)$label[queryHits(ov)],
    other.label    = mcols(other_ranges)$label[subjectHits(ov)],
    
    self.virus     = mcols(self_ranges)$virus[queryHits(ov)],
    other.virus    = mcols(other_ranges)$virus[subjectHits(ov)],
    
    self.probe     = mcols(self_ranges)$probe[queryHits(ov)],
    other.probe    = mcols(other_ranges)$probe[subjectHits(ov)],
  )
  
  # 4. Conditionally add `name` and `ID` columns if they exist in the GFF3 file
  if ("name" %in% names(mcols(self_ranges))) {
    pairs_df$self.domain <- mcols(self_ranges)$name[queryHits(ov)]
    pairs_df$other.domain <- mcols(other_ranges)$name[subjectHits(ov)]
  }
  
  if ("ID" %in% names(mcols(self_ranges))) {
    pairs_df$self.ID <- mcols(self_ranges)$ID[queryHits(ov)]
    pairs_df$other.ID <- mcols(other_ranges)$ID[subjectHits(ov)]
  }
  
  # 5. Add a column to indicate the distance between the two ranges
  pairs_df$distance <- abs(pairs_df$self.start - pairs_df$other.start)
  
  return(pairs_df)
}

pair_df <- find_pairs(ranges, max_gap)


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
write.csv(pair_df, file.path(args$table_output_dir, paste0(file.basename, ".csv")), row.names = FALSE)
