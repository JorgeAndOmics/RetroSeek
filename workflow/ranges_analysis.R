
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

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10) {
  stop("Please provide exactly ten arguments as follows:
  Options:
  fasta=<file>               Path to the genome FASTA input file
  enervate=<file>            Path to the enERVate Parquet input file
  ltrdigest=<file>           Path to the LTRdigest GFF input file
  bitscore_threshold=<val>   Bitscore threshold for filtering
  identity_threshold=<val>   Identity threshold for filtering
  ltr-resize=<val>           Amount to resize LTRdigest ranges
  original_ranges=<file>     Path to the Original Ranges output file
  candidate_ranges=<file>    Path to the Candidate Ranges output file
  valid_ranges=<file>        Path to the Validated Ranges output file
  overlap_matrix=<file>      Path to the Overlap Matrix output file")
}

args.fasta <- args[1]
args.enervate <- args[2]
args.ltrdigest <- args[3]
args.bitscore_threshold <- as.numeric(args[4])
args.identity_threshold <- as.numeric(args[5])
args.ltr_resize <- as.numeric(args[6])
args.original_ranges <- args[7]
args.candidate_ranges <- args[8]
args.valid_ranges <- args[9]
args.overlap_matrix <- args[10]



# Message
print(paste0('Processing ranges for ', basename(args.enervate), '...'))


# Define if it's a main file
is_main <- grepl("_main", args.enervate) & !grepl("_accessory", args.enervate)

# Original FASTA file
fa_file <- Biostrings::readDNAStringSet(args.fasta)
chrom_lengths <- setNames(width(fa_file), names(fa_file))

# Clean up chrom_lengths
names(chrom_lengths) <- stringr::str_extract(names(chrom_lengths), "^[A-Za-z]+_?[0-9]+\\.[0-9]{1,2}")

# enERVate input file
data <- arrow::read_parquet(args.enervate)


# LTRdigest input file
ltr_data <- rtracklayer::import(args.ltrdigest, format = "gff3")


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
bitscore_threshold <- as.numeric(args.bitscore_threshold)
identity_threshold <- as.numeric(args.identity_threshold)

gr <- gr %>%
         filter(bitscore > bitscore_threshold,
                identity > identity_threshold)

seqinfo(gr) <- Seqinfo(seqnames = seqlevels(gr), genome = gr$species[1])

seqlengths(gr) <- chrom_lengths[names(seqlengths(gr))]


# Reduce overlapping ranges for the single species. Collapsing #1
reducing.gr <- function (gr) {
  gr %>%
    group_by(probe) %>%
    reduce_ranges_directed(
      species = paste(unique(species), collapse = "; "),
      virus = paste(sort(unique(virus)), collapse = "; "),
      family = paste(sort(unique(family)), collapse = "; "),
      mean_bitscore = mean(bitscore),
      mean_identity = mean(identity),
      type = "proviral_sequence") %>%
    arrange(.by_group = start)
}

reduced_gr <- reducing.gr(gr)

# Improving readibility of the reduced ranges
reduced_gr$virus <- map(reduced_gr$virus, ~as.list(.x))
reduced_gr$family <- map(reduced_gr$family, ~as.list(.x))

# Naming the reduced ranges by probe findings
named_reduced_gr <- reduced_gr %>% 
  group_by(probe) %>%
  mutate(ID = paste0(probe, "_", seq_along(probe))) %>%
  ungroup()


## LTRDIGEST DATA

# Filter only LTR_retrotransposon full regions
ltr_retro <- ltr_data %>%
  filter(type == "LTR_retrotransposon")

# Expand LTR Retrotrasposon ranges by an amount to rescue potential hits
flank_resize <- as.numeric(args.ltr_resize)

ltr_retro <- ltr_retro %>%
  resize(width = width(ltr_retro) + flank_resize, fix = "center")


# Filter protein domains
ltr_domain <- ltr_data %>%
  filter(!is.na(name)) %>%
  select(-source, -phase, -ID, -ltr_similarity, -seq_number, -reading_frame)

ltr_domain <- ltr_domain %>%
  mutate(probe = case_when(
    grepl("ase|RVT_1|RVT_2|RVT_thumb|rve|IN_DBD_C", name) ~ "POL",
    grepl("Gag|gag|GAG|zf|PTAP|YPXL", name) ~ "GAG",
    grepl("coat|FP|HR1|HR2", name) ~ "ENV",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(probe))  # Keep if `name` has NAs that need filtering


# Find overlaps between enERVate and LTRdigest
ov.E2L <- findOverlaps(named_reduced_gr, ltr_retro)
ov.L2E <- findOverlaps(ltr_retro, named_reduced_gr)

EonL <- named_reduced_gr[unique(queryHits(ov.E2L))]
LonE <- ltr_retro[unique(queryHits(ov.L2E))]
EoutsideL <- named_reduced_gr[-unique(queryHits(ov.E2L))]
LoutsideE <- ltr_retro[-unique(queryHits(ov.L2E))]

percent_EonL <- round(length(EonL) / length(named_reduced_gr) * 100, 2)
percent_LonE <- round(length(LonE) / length(ltr_retro) * 100, 2)
percent_EoutsideL <- round(length(EoutsideL) / length(named_reduced_gr) * 100, 2)
percent_LoutsideE <- round(length(LoutsideE) / length(ltr_retro) * 100, 2)


# Calculate numbers and percentages of overlap for full enERVate
enervate_overlap_df <- data.frame(matrix(NA, nrow = 1, ncol = 4))
rownames(enervate_overlap_df) <- c('enERVate')
colnames(enervate_overlap_df) <- c('In', 'Out', 'In%', 'Out%')
enervate_overlap_df[1, 1] <- length(EonL)
enervate_overlap_df[1, 2] <- length(EoutsideL)
enervate_overlap_df[1, 3] <- percent_EonL
enervate_overlap_df[1, 4] <- percent_EoutsideL


# Calculate numbers and percentages of overlap for full LTRdigest
ltr.full_overlap_df <- data.frame(matrix(NA, nrow = 1, ncol = 4))
rownames(ltr.full_overlap_df) <- c('LTRdigest')
colnames(ltr.full_overlap_df) <- c('In', 'Out', 'In%', 'Out%')
ltr.full_overlap_df[1, 1] <- length(LonE)
ltr.full_overlap_df[1, 2] <- length(LoutsideE)
ltr.full_overlap_df[1, 3] <- percent_LonE
ltr.full_overlap_df[1, 4] <- percent_LoutsideE


# Calculate numbers and percentages of overlap for probes over LTRdigest
probe_overlap_calculator <- function(gr, ltr_retro){
  
  overlap_df_probes <- data.frame(matrix(NA, nrow = length(unique(gr$probe)), ncol = 4))
  rownames(overlap_df_probes) <- c(unique(gr$probe))
  colnames(overlap_df_probes) <- c('In', 'Out', 'In%', 'Out%')
  
  for(pr in unique(gr$probe)){
    ov.E2L <- findOverlaps(gr[gr$probe == pr], ltr_retro)
    ov.L2E <- findOverlaps(ltr_retro, gr[gr$probe == pr])
    
    EonL <- gr[gr$probe == pr][unique(queryHits(ov.E2L))]
    LonE <- ltr_retro[unique(queryHits(ov.L2E))]
    EoutsideL <- gr[gr$probe == pr][-unique(queryHits(ov.E2L))]
    LoutsideE <- ltr_retro[-unique(queryHits(ov.L2E))]
    
    percent_EonL <- round(length(EonL) / length(gr[gr$probe == pr]) * 100, 2)
    percent_LonE <- round(length(LonE) / length(ltr_retro) * 100, 2)
    percent_EoutsideL <- round(length(EoutsideL) / length(gr[gr$probe == pr]) * 100, 2)
    percent_LoutsideE <- round(length(LoutsideE) / length(ltr_retro) * 100, 2)
    
    overlap_df_probes[pr, 1] <- length(EonL)
    overlap_df_probes[pr, 2] <- length(EoutsideL)
    overlap_df_probes[pr, 3] <- percent_EonL
    overlap_df_probes[pr, 4] <- percent_EoutsideL
  }
  
  return(overlap_df_probes)
  
}
probe_overlap_df <- probe_overlap_calculator(named_reduced_gr, ltr_retro)


# Calculate numbers and percentages of overlap for LTRdigest over probes
ltrdigest_overlap_calculator <- function(gr, ltr_retro){
  
  overlap_df_ltr <- data.frame(matrix(NA, nrow = length(unique(gr$probe)), ncol = 4))
  rownames(overlap_df_ltr) <- c(paste0('LTR_', unique(gr$probe)))
  colnames(overlap_df_ltr) <- c('In', 'Out', 'In%', 'Out%')
  
  for(pr in unique(gr$probe)){
    overlap <- paste0('LTR_', pr)
    
    ov.E2L <- findOverlaps(gr[gr$probe == pr], ltr_retro)
    ov.L2E <- findOverlaps(ltr_retro, gr[gr$probe == pr])
    
    EonL <- gr[gr$probe == pr][unique(queryHits(ov.E2L))]
    LonE <- ltr_retro[unique(queryHits(ov.L2E))]
    EoutsideL <- gr[gr$probe == pr][-unique(queryHits(ov.E2L))]
    LoutsideE <- ltr_retro[-unique(queryHits(ov.L2E))]
    
    percent_EonL <- round(length(EonL) / length(gr[gr$probe == pr]) * 100, 2)
    percent_LonE <- round(length(LonE) / length(ltr_retro) * 100, 2)
    percent_EoutsideL <- round(length(EoutsideL) / length(gr[gr$probe == pr]) * 100, 2)
    percent_LoutsideE <- round(length(LoutsideE) / length(ltr_retro) * 100, 2)
    
    overlap_df_ltr[overlap, 1] <- length(LonE)
    overlap_df_ltr[overlap, 2] <- length(LoutsideE)
    overlap_df_ltr[overlap, 3] <- percent_LonE
    overlap_df_ltr[overlap, 4] <- percent_LoutsideE
  }
  
  return(overlap_df_ltr)
  
}

ltr_overlap_df <- ltrdigest_overlap_calculator(named_reduced_gr, ltr_retro)

# Combine all overlap dataframes
overlap_df <- rbind(enervate_overlap_df, ltr.full_overlap_df, probe_overlap_df, ltr_overlap_df)


## VALID HITS

# Validate hits via LTR Retrotrasposons
ltr_valid_hits <- named_reduced_gr[queryHits(ov.E2L)] %>% arrange(start)

# Find overlaps between enERVate and LTR domains in main files
if (is_main) {
  domain_valid_hits <- gr %>%
    
    # 1) Overlap join with ltr_valid_hits on same strand
    join_overlap_inner_directed(
      ltr_valid_hits,
      suffix = c(".gr", ".ltr_valid")) %>%
    
    # 2) Only keep rows whose probe matches after the merge
    filter(probe.gr == probe.ltr_valid) %>%
    mutate(probe = probe.gr) %>%
    
    # 3) Overlap join the result with ltr_domain on same strand
    join_overlap_inner_directed(
      ltr_domain,
      suffix = c(".merged", ".ltr_domain")) %>%
    
    # 4) Again, keep only rows with matching probe
    filter(probe.merged == probe.ltr_domain) %>%
    
    group_by(probe.ltr_valid) %>%
    
    # 5) Finally, reduce
    reduce_ranges_directed(
      probe         = paste(unique(probe.ltr_valid), collapse = "; "),
      species       = paste(unique(species.gr), collapse = "; "),
      virus         = paste(sort(unique(virus.gr)), collapse = "; "),
      family        = paste(sort(unique(family.gr)), collapse = "; "),
      type.x        = "proviral_sequence",
      type.y        = "protein_match",
      name          = paste(sort(unique(name)), collapse = "; "),
      ID            = paste(sort(unique(ID)), collapse = "; "),
      mean_bitscore = paste(unique(mean_bitscore), collapse = "; "),
      mean_identity = paste(unique(mean_identity), collapse = "; ")
    ) %>%
    select(-probe.ltr_valid)
}

# Export hits to GFF3
rtracklayer::export(gr, args.original_ranges, format = "gff3")
rtracklayer::export(ltr_valid_hits, args.candidate_ranges, format = "gff3")
if (is_main) {
rtracklayer::export(domain_valid_hits, args.valid_ranges, format = "gff3")
}

# Export matrices to CSV
write.csv(overlap_df, args.overlap_matrix, row.names = TRUE)