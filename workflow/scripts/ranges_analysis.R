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
parser$add_argument("--fasta", required=TRUE, help="Genome FASTA input file")
parser$add_argument("--blast", required=TRUE, help="tBLASTn Parquet input file")
parser$add_argument("--ltrdigest", required=TRUE, help="LTRdigest GFF3 input file")
parser$add_argument("--probes", required=TRUE, help="Probes metadata file (CSV)")
parser$add_argument("--config", required=TRUE, help="Configuration YAML file")
parser$add_argument("--original_ranges", required=TRUE, help="GFF3 output: original reduced and merged hits")
parser$add_argument("--candidate_ranges", required=TRUE, help="GFF3 output: hits overlapping LTRs")
parser$add_argument("--valid_ranges", required=TRUE, help="GFF3 output: domain-validated hits")
parser$add_argument("--solo_ltr_ranges", required=FALSE, help="GFF3 output: solo LTRs")
parser$add_argument("--flanking_ltr_ranges", required=FALSE, help="GFF3 output: flanking LTRs")
parser$add_argument("--overlap_matrix", required=TRUE, help="CSV summary of overlaps")
parser$add_argument("--plot_dataframe", required=TRUE, help="Parquet output for plotting")

args <- parser$parse_args()

# Print which tBLASTn file is being processed
print(paste0("Processing ranges for ", tools::file_path_sans_ext(basename(args$blast)), "..."))

# -----------------------------
# 3. CONFIGURATION PARAMETERS
# -----------------------------
# Read YAML configuration file
config <- yaml::read_yaml(args$config)

# Import thresholds and merging behavior from config
probe_min_length   <- unlist(config$parameters$probe_min_length)
bitscore_threshold <- as.numeric(config$parameters$bitscore_threshold) %||% 0
identity_threshold <- as.numeric(config$parameters$identity_threshold) %||% 0
ltr_resize         <- as.numeric(config$parameters$ltr_resize) %||% 0
merge_options      <- config$parameters$merge_options %||% "virus"

# -----------------------------
# 4. DETERMINE FILE TYPE
# -----------------------------
# Whether this is a _main file (used for optional domain validation and ltr parsing)
is_main <- grepl("_main", args$blast) & !grepl("_accessory", args$blast)

# -----------------------------
# 5. LOAD FASTA AND CHR LENGTHS
# -----------------------------
fa_file <- Biostrings::readDNAStringSet(args$fasta)
chrom_lengths <- setNames(width(fa_file), names(fa_file))
names(chrom_lengths) <- stringr::str_extract(names(chrom_lengths), "^[A-Za-z]+_?[0-9]+\\.[0-9]{1,2}")

# -----------------------------
# 6. LOAD INPUT FILES
# -----------------------------
data <- arrow::read_parquet(args$blast)
ltr_data <- rtracklayer::import(args$ltrdigest, format = "gff3")

# -----------------------------
# 7. BUILD GRanges OBJECT
# -----------------------------
gr <- GRanges(
  seqnames = data$accession,
  ranges   = IRanges(start = data$hsp_sbjct_start, end = data$hsp_sbjct_end),
  strand   = data$strand
)

# Attach metadata to GRanges object
mcols(gr)$label    <- data$label
mcols(gr)$virus    <- data$virus
mcols(gr)$bitscore <- data$hsp_bits
mcols(gr)$identity <- (data$hsp_identity / data$hsp_align_length) * 100
mcols(gr)$species  <- data$species
mcols(gr)$probe    <- data$probe

# -----------------------------
# 8. FILTER LOW-QUALITY RANGES
# -----------------------------
gr <- gr %>% filter(
  width(.) > ifelse(!is.na(probe_min_length[as.character(probe)]), probe_min_length[as.character(probe)], 0),
  bitscore > bitscore_threshold,
  identity > identity_threshold
)

# -----------------------------
# 9. ADD SEQUENCE METADATA
# -----------------------------
seqinfo(gr) <- Seqinfo(seqnames = seqlevels(gr), genome = gr$species[1])
seqlengths(gr) <- chrom_lengths[names(seqlengths(gr))]

# -----------------------------
# 10. REDUCE BY PROBE + GAP WIDTH
# -----------------------------
gap_vals <- ifelse(!is.na(probe_min_length[as.character(gr$probe)]), probe_min_length[as.character(gr$probe)], 0)
mcols(gr)$min_gapwidth <- gap_vals
gr_list <- split(gr, ~ probe)

# Group-specific reduction using reduce_ranges_directed
if (merge_options == "virus") {
  gr_list_reduced <- sapply(gr_list, function(sub_gr) {
    gap_val <- unique(sub_gr$min_gapwidth)
    sub_gr %>% group_by(probe, virus) %>% reduce_ranges_directed(
      min.gapwidth = gap_val,
      virus    = as.character(unique(virus)),
      species  = as.character(unique(species)),
      label   = as.character(unique(label)),
      bitscore = max(bitscore),
      identity = max(identity)
    )
  })
} else if (merge_options == "label") {
  gr_list_reduced <- sapply(gr_list, function(sub_gr) {
    gap_val <- unique(sub_gr$min_gapwidth)
    sub_gr %>% group_by(probe, label) %>% reduce_ranges_directed(
      min.gapwidth = gap_val,
      virus    = paste(unique(virus), collapse = "; "),
      species  = as.character(unique(species)),
      bitscore = max(bitscore),
      identity = max(identity)
    )
  })
}

# Combine all reduced probe groups
gr <- bind_ranges(gr_list_reduced)

# -----------------------------
# 11. GLOBAL REDUCTION ACROSS PROBES
# -----------------------------
reducing.gr <- function(gr) {
  gr %>%
    group_by(probe) %>%
    reduce_ranges_directed(
      species        = paste(unique(species), collapse = "; "),
      virus          = paste(sort(unique(virus)), collapse = "; "),
      label          = paste(sort(unique(label)), collapse = "; "),
      mean_bitscore  = if (length(bitscore) > 1) mean(bitscore) else bitscore,
      mean_identity  = if (length(identity) > 1) mean(identity) else identity,
      type           = "proviral_sequence"
    ) %>%
    arrange(.by_group = start)
}

reduced_gr <- reducing.gr(gr)
named_reduced_gr <- reduced_gr %>% group_by(probe) %>% mutate(ID = paste0(probe, "_", seq_along(probe))) %>% ungroup()

# -----------------------------
# 12. PROCESS LTR RETRO ELEMENTS
# -----------------------------
ltr_retro <- ltr_data %>% filter(type == "repeat_region")
flank_resize <- as.numeric(ltr_resize)
ltr_retro <- ltr_retro %>% resize(width = width(ltr_retro) + flank_resize, fix = "center")


# ------------------------------------------
# 13. PROCESS LTR PROTEIN DOMAINS FROM LTRDIGEST
# ------------------------------------------
# Step 1: Filter LTRdigest GFF3 entries to keep only those with a 'name' attribute.
# These correspond to protein domains (e.g., RVT, Integrase, etc).
# Remove metadata columns that are not required for downstream analysis.
ltr_domain <- ltr_data %>%
  filter(!is.na(name)) %>%  # Keep only rows that have a named domain
  select(
    -source, -phase,        # Remove unused GFF fields
    -ID,                    # We will aggregate and reassign IDs later
    -ltr_similarity,        # Remove LTR similarity metadata
    -seq_number,            # Redundant field
    -reading_frame          # Not needed for domain-level comparisons
  )

# Step 2: Create a domain-to-probe mapping dictionary.
# Each probe in config$domains is a vector of regular expressions to look for; we collapse each into a single pattern string (x|y|z).
# This creates a named character vector where each name is a probe and value is a regex pattern (POL: "rvt|integrase|aspartic").
domain_map <- purrr::map_chr(config$domains, ~ paste(.x, collapse = "|"))

# Step 3: Function to assign a 'probe' label to a given domain name.
# It searches the domain name against all domain_map regex patterns.
# If multiple probes match, it chooses the longest matching probe name (favoring specificity).
assign_probe <- function(domain_name) {
  matched_probe <- purrr::keep(names(domain_map), function(probe) {
    grepl(domain_map[[probe]], domain_name, ignore.case = TRUE)
  })
  
  if (length(matched_probe) > 0) {
    # Prefer the most specific probe (longest name)
    matched_probe[order(nchar(matched_probe), decreasing = TRUE)[1]]
  } else {
    # Return NA if no match is found
    NA_character_
  }
}

# Step 4: Apply the probe assignment function to each domain name.
# This creates a new 'probe' column identifying which probe the domain matches.
# Then we remove rows that couldn't be assigned to any probe (i.e., unrecognized domains).
ltr_domain <- ltr_domain %>%
  mutate(probe = purrr::map_chr(name, assign_probe)) %>%
  filter(!is.na(probe))  # Retain only domains successfully mapped to a probe

# -----------------------------
# 13B. EXTRACT SOLO LTRs
# -----------------------------
# Main and Accessory will provide the same tracks, as they don't depend on BLAST
# results, so we only perform it on main tracks
if (is_main) {
  # Step 1: Extract all LTR features and their Parent IDs
  ltr_seqs <- ltr_data[ltr_data$type == "long_terminal_repeat"]
  ltr_seqs$ParentID <- sapply(mcols(ltr_seqs)$Parent, function(x) sub(".*Parent=([^;]+).*", "\\1", x))
  
  # Step 2: Extract all LTR_retrotransposon parent IDs (composite ERVs)
  ervs <- ltr_data[ltr_data$type == "LTR_retrotransposon"]
  erv_ids <- unique(mcols(ervs)$ID)

  # Step 3: Filter out LTRs that are part of full ERVs (i.e., retain only solo LTRs)
  solo_ltr <- ltr_seqs[!ltr_seqs$ParentID %in% erv_ids]
  
  # Step 4: Assign unique ID to each solo LTR
  if (length(solo_ltr) > 0) {
  solo_ltr$ID <- paste0("soloLTR_", seq_along(solo_ltr))
  } 
} else {
  solo_ltr <- GRanges()  # Empty GRanges if no solo LTRs
}

# -----------------------------
# 13C. EXTRACT FLANKING LTRs FROM ERVs
# -----------------------------
if (is_main) {
  # Step 1: Keep only LTRs that belong to full ERVs
  ltr_flanking <- ltr_seqs[ltr_seqs$ParentID %in% erv_ids]
  
  # Step 2: Join strand/start from parent ERVs for orientation
  erv_meta <- data.frame(
    ID = mcols(ervs)$ID,
    strand = as.character(strand(ervs)),
    start = start(ervs)
  )
  
  # Step 3: Add orientation and relative position info
  mcols(ltr_flanking)$strand_erv <- erv_meta$strand[match(ltr_flanking$ParentID, erv_meta$ID)]
  mcols(ltr_flanking)$start_erv <- erv_meta$start[match(ltr_flanking$ParentID, erv_meta$ID)]
  
  # Step 4: Label LTR side (left/right) heuristically
  ltr_flanking$LTR_side <- ifelse(
    (ltr_flanking$strand_erv == "+" & start(ltr_flanking) <= ltr_flanking$start_erv + config$parameters$ltr_flank_margin),
    "left",
    ifelse(
      (ltr_flanking$strand_erv == "-" & start(ltr_flanking) >= ltr_flanking$start_erv - config$parameters$ltr_flank_margin),
      "left",
      "right"
    )
  )
  
  # Step 5: Assign ID based on ERV and LTR side
  if (length(ltr_flanking) > 0) {
    ltr_flanking$ID <- paste0("flankLTR_", ltr_flanking$ParentID, "_", ltr_flanking$LTR_side)
  } 
} else {
  ltr_flanking <- GRanges()   # Empty GRanges if no flanking LTRs
}

# -----------------------------
# 14. OVERLAP DETECTION: tBLASTn vs. LTRs
# -----------------------------

# Compute overlaps between tBLASTn hits and LTR retrotransposons
ov.B2L <- findOverlaps(named_reduced_gr, ltr_retro)
ov.L2B <- findOverlaps(ltr_retro, named_reduced_gr)

# Extract overlapping/non-overlapping hits
BonL       <- named_reduced_gr[unique(queryHits(ov.B2L))]    # tBLASTn hits overlapping LTRs
BoutsideL  <- named_reduced_gr[-unique(queryHits(ov.B2L))]   # tBLASTn hits outside LTRs
LonB       <- ltr_retro[unique(queryHits(ov.L2B))]           # LTRs overlapping tBLASTn hits
LoutsideB  <- ltr_retro[-unique(queryHits(ov.L2B))]          # LTRs with no tBLASTn overlap

# Compute percentages of overlap
percent_BonL       <- round(length(BonL) / length(named_reduced_gr) * 100, 2)
percent_BoutsideL  <- round(length(BoutsideL) / length(named_reduced_gr) * 100, 2)
percent_LonB       <- round(length(LonB) / length(ltr_retro) * 100, 2)
percent_LoutsideB  <- round(length(LoutsideB) / length(ltr_retro) * 100, 2)


# -----------------------------
# 15. OVERLAP SUMMARY TABLES
# -----------------------------

# Summary: tBLASTn intervals
blast_overlap_df <- data.frame(
  In   = length(BonL),
  Out  = length(BoutsideL),
  InP  = percent_BonL,
  OutP = percent_BoutsideL,
  row.names = "tBLASTn"
)

# Summary: LTR intervals
ltr.full_overlap_df <- data.frame(
  In   = length(LonB),
  Out  = length(LoutsideB),
  InP  = percent_LonB,
  OutP = percent_LoutsideB,
  row.names = "LTRdigest"
)


# -----------------------------
# 16. PER-PROBE OVERLAP ANALYSIS
# -----------------------------

probe_overlap_calculator <- function(gr, ltr_retro) {
  probes <- unique(gr$probe)
  overlap_df <- data.frame(matrix(NA, nrow = length(probes), ncol = 4))
  rownames(overlap_df) <- probes
  colnames(overlap_df) <- c("In", "Out", "InP", "OutP")
  
  for (pr in probes) {
    gr_pr <- gr[gr$probe == pr]
    ov <- findOverlaps(gr_pr, ltr_retro)
    
    in_hits <- gr_pr[unique(queryHits(ov))]
    out_hits <- gr_pr[-unique(queryHits(ov))]
    
    overlap_df[pr, ] <- c(
      length(in_hits),
      length(out_hits),
      round(length(in_hits) / length(gr_pr) * 100, 2),
      round(length(out_hits) / length(gr_pr) * 100, 2)
    )
  }
  
  return(overlap_df)
}

probe_overlap_df <- probe_overlap_calculator(named_reduced_gr, ltr_retro)


# -----------------------------
# 17. PER-LTR OVERLAP ANALYSIS PER PROBE
# -----------------------------

ltrdigest_overlap_calculator <- function(gr, ltr_retro) {
  probes <- unique(gr$probe)
  overlap_df <- data.frame(matrix(NA, nrow = length(probes), ncol = 4))
  rownames(overlap_df) <- paste0("LTR_", probes)
  colnames(overlap_df) <- c("In", "Out", "InP", "OutP")
  
  for (pr in probes) {
    gr_pr <- gr[gr$probe == pr]
    ov <- findOverlaps(ltr_retro, gr_pr)
    
    in_hits <- ltr_retro[unique(queryHits(ov))]
    out_hits <- ltr_retro[-unique(queryHits(ov))]
    
    overlap_df[paste0("LTR_", pr), ] <- c(
      length(in_hits),
      length(out_hits),
      round(length(in_hits) / length(ltr_retro) * 100, 2),
      round(length(out_hits) / length(ltr_retro) * 100, 2)
    )
  }
  
  return(overlap_df)
}

ltr_overlap_df <- ltrdigest_overlap_calculator(named_reduced_gr, ltr_retro)


# -----------------------------
# 18. COMBINE OVERLAP RESULTS
# -----------------------------

# Merge global and per-probe overlap dataframes
overlap_df <- rbind(blast_overlap_df, ltr.full_overlap_df, probe_overlap_df, ltr_overlap_df)


# -----------------------------
# 19. CONSTRUCT PLOTTING DATAFRAME
# -----------------------------

# Read virus metadata
probe_df <- read.csv2(args$probes, header = TRUE, sep = ",")
probe_df_sum <- probe_df %>% distinct(Name, Label, Abbreviation)

# Prepare plot-friendly dataframe from reduced GRanges
plot_df <- as.data.frame(reduced_gr) %>%
  tidyr::separate_rows(virus, sep = "; ") %>%
  mutate(
    label = probe_df_sum$Label[match(virus, probe_df_sum$Name)],
    abbreviation = probe_df_sum$Abbreviation[match(virus, probe_df_sum$Name)]
  )


# -----------------------------
# 20. IDENTIFY VALID tBLASTn HITS
# -----------------------------

# Candidate hits are those overlapping LTR regions
ltr_valid_hits <- named_reduced_gr[queryHits(ov.B2L)] %>% arrange(start)


# -----------------------------
# 21. OPTIONAL DOMAIN VALIDATION
# -----------------------------

domain_valid_hits <- gr %>%
  
  # Step 1: Overlap tBLASTn hits with valid LTR hits (strand-aware)
  join_overlap_inner_directed(
    ltr_valid_hits,
    suffix = c(".gr", ".ltr_valid")
  ) %>%
  
  # Step 2: Keep only pairs where probes match between tBLASTn and LTR
  filter(probe.gr == probe.ltr_valid) %>%
  mutate(probe = probe.gr) %>%
  
  # Step 3: Overlap the result with LTR domains (strand-aware again)
  join_overlap_inner_directed(
    ltr_domain,
    suffix = c(".merged", ".ltr_domain")
  ) %>%
  
  # Step 4: Filter only cases where probe labels match again
  filter(probe.merged == probe.ltr_domain) %>%
  mutate(Parent = as.character(Parent)) %>%
  
  # Step 5: Group by probe for reduction
  group_by(probe.ltr_valid) %>%
  
  # Step 6: Reduce overlapping domains, aggregating relevant fields
  reduce_ranges_directed(
    probe         = paste(unique(probe.ltr_valid), collapse = "; "),
    species       = paste(unique(species.gr), collapse = "; "),
    virus         = paste(sort(unique(virus.gr)), collapse = "; "),
    label        = paste(sort(unique(label.gr)), collapse = "; "),
    name          = paste(sort(unique(name)), collapse = "; "),
    ID            = paste(sort(unique(ID)), collapse = "; "),
    Parent        = paste(unique(Parent), collapse = "; "),
    mean_bitscore = paste(unique(mean_bitscore), collapse = "; "),
    mean_identity = paste(unique(mean_identity), collapse = "; "),
  ) %>%
  
  # Step 7: Clean up and order result
  select(-probe.ltr_valid) %>%
  arrange(.by_group = start)


# -----------------------------
# 22. EXPORT OUTPUT FILES
# -----------------------------
# Export GRanges in GFF3 format. If empty GRanges, it will create a GFF3 file with only the header.

track_exporter <- function(track, path) {
  if (length(track) > 0) {
    rtracklayer::export(track, path, format = "gff3")
  } else {
    writeLines("##gff-version 3", con = path)
  }
}

track_exporter(named_reduced_gr, args$original_ranges)
track_exporter(ltr_valid_hits, args$candidate_ranges)
track_exporter(domain_valid_hits, args$valid_ranges)
if (is_main) {    # Export solo LTRs and flanking LTRs only if this is a main file
  track_exporter(solo_ltr, args$solo_ltr_ranges)
  track_exporter(ltr_flanking, args$flanking_ltr_ranges)
}

# Export overlap summary and plotting data
arrow::write_parquet(plot_df, args$plot_dataframe)
write.csv(overlap_df, args$overlap_matrix, row.names = TRUE)