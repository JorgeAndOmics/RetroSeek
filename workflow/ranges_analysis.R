################################################################################
# ============================
# EXTENSIVELY COMMENTED SCRIPT
# ============================
# Below is the original script with additional comprehensive comments explaining
# the purpose and functionality of each section and line. No changes were made
# to the original code itself (no modifications to logic, no renaming of
# variables, etc.). All that follows are clarifying remarks and commentary.
################################################################################

# -------------------
# 1. DEPENDENCIES
# -------------------

# Suppress warnings so that the output is cleaner; warnings won't appear on-screen.
options(warn = -1)

# Suppress messages from loading libraries, so the console remains uncluttered.
suppressMessages({
  library(arrow)         # Provides tools for reading and writing Parquet files
  library(tidyverse)     # Collection of R packages for data manipulation and visualization
  library(docopt)        # Library for command-line interface argument parsing
  library(GenomicRanges) # For representing and manipulating genomic intervals
  library(plyranges)     # "Tidy"-style genomic data manipulation using GRanges
  library(rtracklayer)   # For reading and writing genome annotation files (e.g., GFF, BED)
})


# -------------------------------
# 2. PARSE COMMAND-LINE ARGUMENTS
# -------------------------------
# The script expects exactly 13 arguments. If fewer or more are found, it stops
# and prints an error message showing the user the correct usage.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 14) {
  stop("Please provide exactly fourteen arguments as follows:
  Options:
  fasta=<file>               Path to the genome FASTA input file
  enervate=<file>            Path to the enERVate Parquet input file
  ltrdigest=<file>           Path to the LTRdigest GFF input file
  probes=<file>              Path to the probes input file
  defaults=<path>            Path to the defaults.py directory
  bitscore_threshold=<val>   Bitscore threshold for filtering
  identity_threshold=<val>   Identity threshold for filtering
  ltr_resize=<val>           Amount to resize LTRdigest ranges
  merge_options=<val>        Options for merging ranges (species or virus)
  original_ranges=<file>     Path to the Original Ranges output file
  candidate_ranges=<file>    Path to the Candidate Ranges output file
  valid_ranges=<file>        Path to the Validated Ranges output file
  overlap_matrix=<file>      Path to the Overlap Matrix output file
  plot_dataframe=<file>      Path to the Plot DataFrame output file")
}

# Each command-line argument is assigned to a corresponding variable for clarity.
args.fasta <- args[1]
args.enervate <- args[2]
args.ltrdigest <- args[3]
args.probes <- args[4]
args.defaults <- args[5]
args.bitscore_threshold <- as.numeric(args[6])
args.identity_threshold <- as.numeric(args[7])
args.ltr_resize <- as.numeric(args[8])
args.merge_options <- args[9]
args.original_ranges <- args[10]
args.candidate_ranges <- args[11]
args.valid_ranges <- args[12]
args.overlap_matrix <- args[13]
args.plot_dataframe <- args[14]


# This line prints a progress message showing which enERVate file is being processed.
print(paste0('Processing ranges for ', basename(args.enervate), '...'))


# -------------------------------------------------
# 3. DETERMINE IF WE ARE DEALING WITH A MAIN FILE
# -------------------------------------------------
# Check if the filename for the enERVate file includes "_main" but does not
# include "_accessory".

is_main <- grepl("_main", args.enervate) & !grepl("_accessory", args.enervate)


# ----------------------------------------------
# 4. PREPARE FASTA SEQUENCES AND CHROM LENGTHS
# ----------------------------------------------
# Read the genome FASTA input file into a DNAStringSet.

fa_file <- Biostrings::readDNAStringSet(args.fasta)

# Gather the lengths of each sequence in the FASTA. The names of the 'chrom_lengths'
# vector will be used to match chrom names in the GRanges objects.
chrom_lengths <- setNames(width(fa_file), names(fa_file))

# Clean up the chromosome names by applying a regex that captures typical
# accession naming schemes like "ABC_123.45".
names(chrom_lengths) <- stringr::str_extract(names(chrom_lengths), "^[A-Za-z]+_?[0-9]+\\.[0-9]{1,2}")


# -------------------------------------
# 5. ENERVATE INPUT PARQUET PROCESSING
# -------------------------------------
# Read the enERVate data from a Parquet file using the 'arrow' library.

data <- arrow::read_parquet(args.enervate)


# ----------------------------------------------
# 6. LTRDIGEST INPUT FILE AND LTR RETRO ELEMENTS
# ----------------------------------------------
# Import the LTRdigest GFF3 file.

ltr_data <- rtracklayer::import(args.ltrdigest, format = "gff3")


# ------------------------------------------
# 7. CREATE GRANGES OBJECT FROM ENERVATE DATA
# ------------------------------------------
# Convert the enERVate table into a GRanges object where each row becomes
# a genomic interval with associated metadata.

gr <- GRanges(
  seqnames = data$accession,                # Chromosome/accession name
  ranges   = IRanges(
    start = data$hsp_sbjct_start,
    end   = data$hsp_sbjct_end
  ),
  strand   = data$strand                    # Strand information
)

# Assign relevant metadata (columns) to each range in the GRanges object:
mcols(gr)$family <- data$family
mcols(gr)$virus <- data$virus
mcols(gr)$bitscore <- data$hsp_bits
mcols(gr)$identity <- (data$hsp_identity / data$hsp_align_length) * 100
mcols(gr)$species <- data$species
mcols(gr)$probe <- data$probe


# ------------------------------------------------------
# 8. FILTER RANGES BY BITSCORE AND IDENTITY THRESHOLDS
# ------------------------------------------------------
# The user-specified bitscore_threshold and identity_threshold define how stringent
# we are in filtering out low-quality or partial hits in the genomic data.

bitscore_threshold <- as.numeric(args.bitscore_threshold)
identity_threshold <- as.numeric(args.identity_threshold)

# Use plyranges 'filter' to keep only those ranges with bitscore and identity
# exceeding the specified thresholds.
gr <- gr %>%
  filter(bitscore > bitscore_threshold,
         identity > identity_threshold)


# --------------------------------------------------------
# 9. ATTACH SEQUENCE INFO TO THE GENOMIC RANGES OBJECT
# --------------------------------------------------------
# The 'Seqinfo' in GRanges provides metadata about each chromosome: its length,
# genome version, etc. We add the species from the first row for the 'genome' field.

seqinfo(gr) <- Seqinfo(seqnames = seqlevels(gr), genome = gr$species[1])

# Match up the known chromosome lengths from the FASTA to the ranges in 'gr'.
seqlengths(gr) <- chrom_lengths[names(seqlengths(gr))]


# -----------------------------------------------------------------------
# 10. MERGE RANGES BASED ON PROBE-SPECIFIC DEFAULTS (GAP WIDTH)
# -----------------------------------------------------------------------
# We import 'defaults.py' using 'reticulate'.

defaults <- reticulate::import_from_path("defaults", args.defaults)
probe_min_length <- unlist(defaults$PROBE_MIN_LENGTH)

# For each row in the GRanges object, we lookup the appropriate gap width by
# matching the 'probe' field to 'probe_min_length' in the defaults.
gap_vals <- probe_min_length[as.character(gr$probe)]

# If any probe is missing in the dictionary, replace NA with 0 to avoid errors.
gap_vals[is.na(gap_vals)] <- 0

# Store the gap width in the metadata so we can use it later in reduce_ranges_directed.
mcols(gr)$min_gapwidth <- gap_vals

# -----------------------
# 10a. SPLIT BY PROBE
# -----------------------
# Split the GRanges object into a GRangesList, grouped by 'probe'. This lets us
# apply group-specific merging rules. A known Bioconductor limitation is that it
# let perform dynamically assigned gap widths in 'reduce_ranges_directed'.

gr_list <- split(gr, ~ probe)

# -----------------------
# 10b. REDUCE BY PROBE
# -----------------------
# For each subgroup (i.e., each probe), apply a function that merges overlapping
# or adjacent ranges based on the 'min_gapwidth'. We keep selected metadata fields,
# summarizing bitscore and identity by taking the max for each merging block. Two 
# ways to merge: By species or by family

if (args.merge_options == "species") {
  gr_list_reduced <- sapply(gr_list, function(sub_gr) {
    gap_val <- unique(sub_gr$min_gapwidth)  # Each subgroup has a single gap width
    sub_gr <- sub_gr %>% 
      group_by(probe, virus) %>%           # Group further by probe and virus
      reduce_ranges_directed(
        min.gapwidth = gap_val,
        virus    = as.character(unique(virus)),
        species  = as.character(unique(species)),
        family   = as.character(unique(family)),
        bitscore = max(bitscore),
        identity = max(identity)
      )
  })
}

if (args.merge_options == "virus") {
  gr_list_reduced <- sapply(gr_list, function(sub_gr) {
    gap_val <- unique(sub_gr$min_gapwidth)  # Each subgroup has a single gap width
    sub_gr <- sub_gr %>% 
      group_by(probe, family) %>%           # Group further by probe and virus
      reduce_ranges_directed(
        min.gapwidth = gap_val,
        virus    = paste(unique(virus), collapse = "; "),
        species  = as.character(unique(species)),
        bitscore = max(bitscore),
        identity = max(identity)
      )
  })
}

# -----------------------
# 10c. RECOMBINE RESULTS
# -----------------------
# Combine the reduced GRanges objects from each probe back into a single GRanges.

gr <- bind_ranges(gr_list_reduced)


# -------------------------------------------------
# 11. REDUCE OVERLAPPING RANGES ACROSS ALL PROBES
# -------------------------------------------------
# Use an internal function to reduce (merge) ranges across the entire dataset,
# grouping by 'probe'. Again, we combine overlapping intervals, aggregating some
# metadata by taking the mean, etc.

reducing.gr <- function (gr) {
  gr %>%
    group_by(probe) %>%
    reduce_ranges_directed(
      species       = paste(unique(species), collapse = "; "),
      virus         = as.character(paste(sort(unique(virus)), collapse = "; ")),
      family        = as.character(paste(sort(unique(family)), collapse = "; ")),
      mean_bitscore = mean(bitscore),
      mean_identity = mean(identity),
      type          = "proviral_sequence"
    ) %>%
    arrange(.by_group = start)
}

# We apply this function to our GRanges object to get a set of minimal,
# non-overlapping intervals that represent consolidated hits.
reduced_gr <- reducing.gr(gr)

# We name the reduced ranges by creating a unique ID in each group. For example,
# "probeX_1", "probeX_2", etc.
named_reduced_gr <- reduced_gr %>% 
  group_by(probe) %>%
  mutate(ID = paste0(probe, "_", seq_along(probe))) %>%
  ungroup()


# -----------------------------------------
# 12. LTRDIGEST DATA AND LTR RETROELEMENTS
# -----------------------------------------
# We'll focus on only the main 'LTR_retrotransposon' annotations, which represent
# the predicted retrotransposon intervals. Then, we may want to virtually expand
# (resize) those intervals to capture additional sequences.

ltr_retro <- ltr_data %>%
  filter(type == "LTR_retrotransposon")

# Expand LTR retrotransposons by 'flank_resize' on each side (center-based).
flank_resize <- as.numeric(args.ltr_resize)
ltr_retro <- ltr_retro %>%
  resize(width = width(ltr_retro) + flank_resize, fix = "center")


# -------------------------------------------------
# 13. PROCESS LTR PROTEIN DOMAINS FROM LTRDIGEST
# -------------------------------------------------
# Filter out only lines that have a 'name' (protein domain). We exclude certain
# columns we do not need. Then, assign a 'probe' label based on domain keywords.

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
  filter(!is.na(probe))  # Keep only recognized probe domains


# ---------------------------------------------------------------------
# 14. OVERLAP DETECTION BETWEEN enERVate AND LTR Retrotransposons
# ---------------------------------------------------------------------
# We find intervals from enERVate that overlap with intervals from LTRdigest,
# specifically the expanded 'ltr_retro'.

ov.E2L <- findOverlaps(named_reduced_gr, ltr_retro)
ov.L2E <- findOverlaps(ltr_retro, named_reduced_gr)

# EonL: enERVate intervals that do overlap with LTR intervals.
EonL <- named_reduced_gr[unique(queryHits(ov.E2L))]

# LonE: LTR intervals that do overlap with enERVate intervals.
LonE <- ltr_retro[unique(queryHits(ov.L2E))]

# EoutsideL: enERVate intervals that do NOT overlap with any LTR intervals.
EoutsideL <- named_reduced_gr[-unique(queryHits(ov.E2L))]

# LoutsideE: LTR intervals that do NOT overlap with any enERVate intervals.
LoutsideE <- ltr_retro[-unique(queryHits(ov.L2E))]

# Calculate simple percentages to show the proportion of intervals overlapping
# or not in each dataset.
percent_EonL <- round(length(EonL) / length(named_reduced_gr) * 100, 2)
percent_LonE <- round(length(LonE) / length(ltr_retro) * 100, 2)
percent_EoutsideL <- round(length(EoutsideL) / length(named_reduced_gr) * 100, 2)
percent_LoutsideE <- round(length(LoutsideE) / length(ltr_retro) * 100, 2)


# ----------------------------------------------------------
# 15. SUMMARIZE OVERLAP RESULTS FOR enERVate AND LTRdigest
# ----------------------------------------------------------
# We build small data frames to store absolute counts and percentages of the
# "In" vs "Out" (overlapping vs not-overlapping) intervals for both enERVate
# and LTRdigest.

# enERVate overlap summary
enervate_overlap_df <- data.frame(matrix(NA, nrow = 1, ncol = 4))
rownames(enervate_overlap_df) <- c('enERVate')
colnames(enervate_overlap_df) <- c('In', 'Out', 'In%', 'Out%')
enervate_overlap_df[1, 1] <- length(EonL)
enervate_overlap_df[1, 2] <- length(EoutsideL)
enervate_overlap_df[1, 3] <- percent_EonL
enervate_overlap_df[1, 4] <- percent_EoutsideL

# LTRdigest overlap summary
ltr.full_overlap_df <- data.frame(matrix(NA, nrow = 1, ncol = 4))
rownames(ltr.full_overlap_df) <- c('LTRdigest')
colnames(ltr.full_overlap_df) <- c('In', 'Out', 'In%', 'Out%')
ltr.full_overlap_df[1, 1] <- length(LonE)
ltr.full_overlap_df[1, 2] <- length(LoutsideE)
ltr.full_overlap_df[1, 3] <- percent_LonE
ltr.full_overlap_df[1, 4] <- percent_LoutsideE


# ------------------------------------------------------------
# 16. OVERLAP CALCULATOR FOR EACH PROBE VS. LTR RETROELEMENTS
# ------------------------------------------------------------
# This function calculates, for each unique probe in enERVate, how many intervals
# overlap with LTR retrotransposons, and how many do not, both in counts and as percentages.

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

# Apply it to 'named_reduced_gr' (the consolidated enERVate ranges) vs. 'ltr_retro'.
probe_overlap_df <- probe_overlap_calculator(named_reduced_gr, ltr_retro)


# ------------------------------------------------------------
# 17. OVERLAP CALCULATOR FOR LTRdigest VS. EACH PROBE
# ------------------------------------------------------------
# Similar to the previous function, but from the perspective of LTRdigest intervals
# instead of enERVate intervals. Summarizes how many LTR intervals overlap or not
# for each probe type in 'gr'.

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

# Apply it, then store the result in 'ltr_overlap_df'.
ltr_overlap_df <- ltrdigest_overlap_calculator(named_reduced_gr, ltr_retro)


# ---------------------------------------------------------------
# 18. COMBINE ALL OVERLAP DATAFRAMES INTO A SINGLE MATRIX
# ---------------------------------------------------------------
# We row-bind the summary for enERVate, LTRdigest, and the per-probe overlap results
# into one consolidated table for easy export or inspection.

overlap_df <- rbind(enervate_overlap_df, ltr.full_overlap_df, probe_overlap_df, ltr_overlap_df)


# ---------------------------------------------
# 19. CONSTRUCT A PLOTTING DATAFRAME
# ---------------------------------------------
# We import our probe info CSV file, then we combine that info with the reduced GRanges.

probe_df <- read.csv2(args.probes, header = TRUE, sep = ",")

# Keep only unique (Name, Family, Abbreviation) rows to avoid duplicates.
probe_df_sum <- probe_df %>%
  distinct(Name, Family, Abbreviation)

# Convert the reduced GRanges object to a data frame and separate rows that
# have multiple viruses listed, splitting them by "; ".
df_plot <- as.data.frame(reduced_gr) %>%
  tidyr::separate_rows(virus, sep = "; ")

# Then, map each virus to its family and abbreviation from 'probe_df_sum'.
df_plot <- df_plot %>%
  mutate(
    family       = probe_df_sum$Family[match(virus, probe_df_sum$Name)],
    abbreviation = probe_df_sum$Abbreviation[match(virus, probe_df_sum$Name)]
  )


# -----------------------------
# 20. IDENTIFY VALIDATED HITS
# -----------------------------
# "Validated hits" are those enERVate intervals that overlap with expanded LTR
# retrotransposons. We label these as the safer candidate integrations.

ltr_valid_hits <- named_reduced_gr[queryHits(ov.E2L)] %>% arrange(start)


# ------------------------------------------------------------------
# 21. OPTIONAL DOMAIN VALIDATION (IF MAIN FILE)
# ------------------------------------------------------------------
# We additionally check that overlaps exist between the enERVate intervals, 
# the validated LTR intervals, and the recognized LTR domains, all on the same strand.

if (is_main) {
  domain_valid_hits <- gr %>%
    
    # Overlap join enERVate data with 'ltr_valid_hits' on the same strand.
    join_overlap_inner_directed(
      ltr_valid_hits,
      suffix = c(".gr", ".ltr_valid")
    ) %>%
    
    # Keep only rows that match the same 'probe' after the merge.
    filter(probe.gr == probe.ltr_valid) %>%
    mutate(probe = probe.gr) %>%
    
    # Overlap join the result with the 'ltr_domain' data, again on the same strand.
    join_overlap_inner_directed(
      ltr_domain,
      suffix = c(".merged", ".ltr_domain")
    ) %>%
    
    # Again, filter by matching 'probe' fields after the second merge.
    filter(probe.merged == probe.ltr_domain) %>%
    
    mutate(Parent = as.character(Parent)) %>%
    
    group_by(probe.ltr_valid) %>%
    
    # Consolidate (reduce) these intervals, aggregating relevant fields such as
    # bitscore, identity, and domain names.
    reduce_ranges_directed(
      probe         = paste(unique(probe.ltr_valid), collapse = "; "),
      species       = paste(unique(species.gr), collapse = "; "),
      virus         = paste(sort(unique(virus.gr)), collapse = "; "),
      family        = paste(sort(unique(family.gr)), collapse = "; "),
      origin        = "proviral_sequence",
      type          = "protein_match",
      name          = paste(sort(unique(name)), collapse = "; "),
      ID            = paste(sort(unique(ID)), collapse = "; "),
      Parent        = paste(unique(Parent), collapse = "; "),
      mean_bitscore = paste(unique(mean_bitscore), collapse = "; "),
      mean_identity = paste(unique(mean_identity), collapse = "; ")
    ) %>%
    select(-probe.ltr_valid) %>%
    
    arrange(.by_group = start)
}


# -----------------------------------------------------
# 22. EXPORT RANGES TO GFF3 AND PARQUET FILES
# -----------------------------------------------------
# Write out the original ranges, the candidate ranges that overlap with LTR
# retrotransposons, and the domain-validated hits.

rtracklayer::export(gr, args.original_ranges, format = "gff3")
rtracklayer::export(ltr_valid_hits, args.candidate_ranges, format = "gff3")
if (is_main) {
  rtracklayer::export(domain_valid_hits, args.valid_ranges, format = "gff3")
}

# Convert 'df_plot' into a Parquet file.
arrow::write_parquet(df_plot, args.plot_dataframe)

# Export the overlap matrix as CSV.
write.csv(overlap_df, args.overlap_matrix, row.names = TRUE)
