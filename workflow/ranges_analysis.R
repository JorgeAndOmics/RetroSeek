
# Dependencies


options(warn = -1)  # Suppress warnings
suppressMessages({  # Suppress messages
  library(arrow)
  library(tidyverse)
  library(docopt)
  library(GenomicRanges)
  library(plyranges)
  library(rtracklayer)
})


# Define the command-line interface (CLI) description
doc <- "
Usage:
  script.R --enERVate=<file> --LTRharvest=<file> --candidate_ranges=<file> --validated_ranges=<file> --overlap_matrix=<file> --overlap_percent_matrix=<file>

Options:
  --enERVate=<file>              Path to the enERVate Parquet input file
  --LTRharvest=<file>            Path to the LTRharvest GFF input file
  --candidate_ranges=<file>      Path to the Candidate Ranges output file
  --validated_ranges=<file>      Path to the Validated Ranges output file
  --overlap_matrix=<file>        Path to the Overlap Matrix output file
  --overlap_percent_matrix=<file> Path to the Overlap Percent Matrix output file
"

# Parse command-line arguments
args <- docopt(doc)


# Message
print(paste0('Processing ', basename(args$enERVate), '...'))


# enERVate input file
data <- arrow::read_parquet(args$enERVate)


# LTRharvest input file
ltr_data <- rtracklayer::import(args$LTRharvest, format = "gff3")


## ENERVATE DATA

# Generate GRanges Object
gr <- GRanges(seqnames = data$accession, 
              ranges = IRanges(
                start = data$hsp_sbjct_start,
                end = data$hsp_sbjct_end),
              strand = data$strand)

mcols(gr)$family <- data$family
mcols(gr)$virus <- data$virus
mcols(gr)$bitscore <- data$hsp_bits
mcols(gr)$identity <- (data$hsp_identity / data$hsp_align_length) * 100
mcols(gr)$species <- data$species
mcols(gr)$probe <- data$probe


# Filter ranges by bitscore
bitscore_threshold <- 200
identity_threshold <- 40
gr <- gr %>%
         filter(bitscore > bitscore_threshold,
                identity > identity_threshold)


# Reduce overlapping ranges for the single species. Collapsing #1
reducing.gr <- function (gr) {
  gr %>%
    group_by(probe) %>%
    reduce_ranges_directed(
      species = species[1],
      virus = unique(virus),
      family = unique(family),
      mean_bitscore = mean(bitscore),
      mean_identity = mean(identity),
      type = "proviral_sequence") %>%
    arrange(.by_group = start)
}

reduced_gr <- reducing.gr(gr)


# Improving readibility of the reduced ranges
reduced_gr$virus <- as.list(reduced_gr$virus)
reduced_gr$family <- as.list(reduced_gr$family)

# Naming the reduced ranges by probe findings
named_reduced_gr <- reduced_gr %>% 
  group_by(probe) %>%
  mutate(ID = paste0(probe, "_", seq_along(probe))) %>%
  ungroup()


## LTRHARVEST DATA

# Filter only LTR_retrotransposon full regions
ltr_data <- ltr_data %>%
  filter(type == "LTR_retrotransposon")


# Find overlaps between enERVate and LTRharvest
ov.E2L <- findOverlaps(named_reduced_gr, ltr_data)
ov.L2E <- findOverlaps(ltr_data, named_reduced_gr)

EonL <- named_reduced_gr[unique(queryHits(ov.E2L))]
LonE <- ltr_data[unique(queryHits(ov.L2E))]
EoutsideL <- named_reduced_gr[-unique(queryHits(ov.E2L))]
LoutsideE <- ltr_data[-unique(queryHits(ov.L2E))]

percent_EonL <- length(EonL) / length(named_reduced_gr) * 100
percent_LonE <- length(LonE) / length(ltr_data) * 100
percent_EoutsideL <- length(EoutsideL) / length(named_reduced_gr) * 100
percent_LoutsideE <- length(LoutsideE) / length(ltr_data) * 100


# Calculate numbers and percentages of overlap
overlap_df <- data.frame(matrix(NA, nrow = 2, ncol = 2))
rownames(overlap_df) <- c('enERVate.query', 'LTRharvest.query')
colnames(overlap_df) <- c('In', 'Out')
overlap_df[1, 1] <- length(EonL)
overlap_df[1, 2] <- length(EoutsideL)
overlap_df[2, 1] <- length(LonE)
overlap_df[2, 2] <- length(LoutsideE)

percent_overlap_df <- data.frame(matrix(NA, nrow = 2, ncol = 2))
rownames(percent_overlap_df) <- c('enERVate.query', 'LTRharvest.query')
colnames(percent_overlap_df) <- c('In', 'Out')
percent_overlap_df[1, 1] <- paste0(round(percent_EonL, 2), '%')
percent_overlap_df[1, 2] <- paste0(round(percent_EoutsideL, 2), '%')
percent_overlap_df[2, 1] <- paste0(round(percent_LonE, 2), '%')
percent_overlap_df[2, 2] <- paste0(round(percent_LoutsideE, 2), '%')


# Validate hits via LTRharvest
valid_hits <- named_reduced_gr[queryHits(ov.E2L)]


# Export hits to GFF3
rtracklayer::export(named_reduced_gr, args$candidate_ranges, format = "gff3")
rtracklayer::export(valid_hits, args$validated_ranges, format = "gff3")


# Export matrices to CSV
write.csv(overlap_df, args$overlap_matrix, row.names = TRUE)
write.csv(percent_overlap_df, args$overlap_percent_matrix, row.names = TRUE)