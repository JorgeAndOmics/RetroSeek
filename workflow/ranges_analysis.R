
# Dependencies
library(arrow)
library(tidyverse)
library(esquisse)
library(ggsci)
library(ggalluvial)
library(ggdist)

library(plyranges)
library(GenomicRanges)
library(rtracklayer)
library(AnnotationHub)
library(Biostrings)
library(igvR)


# enERVate input file
data <- arrow::read_parquet(
  file.path(
    'C:',
    'Users',
    'Lympha',
    'Documents',
    'Repositories',
    'enERVate',
    'results',
    'tables',
    'segmented_species',
    'Hipposideros_larvatus_main.parquet'
  )
)


# LTRharvest input file
ltr_data <- rtracklayer::import(file.path(
  'C:',
  'Users',
  'Lympha',
  'Documents',
  'Repositories',
  'enERVate',
  'results',
  'ltrharvest',
  'Hipposideros_larvatus_sorted.gff'), format = "gff3")


## ENERVATE DATA

# Columns for Bioproject and Biosample
data <- data %>%
  mutate(bioproject = map_chr(genbank_dbxrefs, ~str_extract(.[1], "(?<=BioProject:)\\w+")),
         biosample = map_chr(genbank_dbxrefs, ~str_extract(.[2], "(?<=BioSample:)\\w+")))

head(data$biosample)


# Calculating distributions
mean_bit <- mean(data$hsp_bits)
q1_bit <- quantile(data$hsp_bits, 0.25)
median_bit <- quantile(data$hsp_bits, 0.5)
q3_bit <- quantile(data$hsp_bits, 0.75)

q1_bit
median_bit
q3_bit


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

gr


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


rtracklayer::export(named_reduced_gr, "Hipposideros_larvatus_main.gff3", format = "gff3")

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

rtracklayer::export(valid_hits, "Hipposideros_larvatus_main_valid_hits.gff3", format = "gff3")
