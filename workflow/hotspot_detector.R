# =============================================================================
# 0. Suppress Warnings and Load Libraries
# =============================================================================
options(warn = -1)
suppressMessages({
  library(Biostrings)
  library(rtracklayer)
  library(GenomicRanges)
  library(regioneR)
  library(tidyverse)
  library(purrr)
  library(ggsci)
})

# TODO: DEDUPLICATOR FAMILY FUNCTION

# =============================================================================
# 1. Define File Paths and Species Name
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 10) {
  stop("Usage: Rscript hotspot_detector.R 
       <genome.fa> <ervs.gff3> <window_size> <mask_size> <mask_mismatch> 
       <perms> <pvalue_threshold> <csv_output_dir> <pdf_output_dir> 
       <hotspot_output_dir>")
}

args.genome <- file.path(args[1])
args.gff <- file.path(args[2])
args.window_size <- as.numeric(args[3])
args.mask_size <- as.numeric(args[4])
args.mask.mismatch <- as.numeric(args[5])
args.perms <- as.numeric(args[6])
args.pvalue_threshold <- as.numeric(args[7])
args.csv_output_dir <- file.path(args[8])
args.pdf_output_dir <- file.path(args[9])
args.hotspot_output_dir <- file.path(args[10])

species_name <- tools::file_path_sans_ext(basename(args.genome))

# =============================================================================
# 2. Load Genome and Create Genome Boundaries
# =============================================================================
# Initialization message
print(paste0("Initializing permutation analysis for ", species_name, " ..."))

# Load the genome as a DNAStringSet
subject_genome <- readDNAStringSet(args.genome)

# Create a GRanges object representing each chromosome's full range
genome_granges <- GRanges(
  seqnames = names(subject_genome),
  ranges   = IRanges(start = 1, end = width(subject_genome))
)

# Create a named numeric vector (chromosome lengths) for tiling
seqinfo.vec <- setNames(width(subject_genome), names(subject_genome))

# =============================================================================
# 3. Create Genomic Windows (Tiles) and Mask Some Unmappable Regions
# =============================================================================
print(paste0("Segmenting ", species_name, " genome in tiles and masking ..."))

window_size <- args.window_size  
genomic_windows <- tileGenome(
  seqlengths = seqinfo.vec,
  tilewidth = window_size,
  cut.last.tile.in.chrom = TRUE
)

mask_size <- args.mask_size
mask_mismatch <- args.mask.mismatch
mask_pattern <- DNAString(paste(rep("N", mask_size), collapse = ""))
mask_match <- vmatchPattern(mask_pattern, subject_genome, max.mismatch = mask_mismatch)  # Avoid some unmappable regions
genome_mask <- as(mask_match, "GRanges")

# =============================================================================
# 4. Define a Randomization Function for ERV Regions
# =============================================================================
# This function randomizes regions on a per-chromosome basis.
# The ellipsis (...) ensures extra arguments passed by permTest() are accepted.
randomize_ervs <- function(A, genome, ...) {
  randomizeRegions(A, 
                   genome = genome,
                   mask = genome_mask,
                   per.chromosome = TRUE)
}

# =============================================================================
# 5. Import and Split GFF3 Data by "family"
# =============================================================================
# Import the GFF3 file
suppressMessages({subject_genome_gff <- rtracklayer::import(args.gff, format = "gff3")})

# Check if the "family" metadata field exists; if not, stop with an error.
if (!"family" %in% colnames(mcols(subject_genome_gff))) {
  stop("The GFF3 file does not contain a 'family' metadata field.")
}

# Convert the imported GFF3 data to a GRanges object (keeping metadata)
subject_genome_gff_clean <- as(subject_genome_gff, "GRanges")

# Split the GRanges object by the "family" field into a GRangesList
gff_by_family <- split(subject_genome_gff_clean, mcols(subject_genome_gff_clean)$family)

# =============================================================================
# 6. Perform Permutation Tests for Each Family (for Overall Overlap)
# =============================================================================
# For each family, run a permutation test comparing overlaps between genomic windows and the family’s integration events.
n_perms = args.perms

perm_results <- map2(
  names(gff_by_family), gff_by_family,
  ~{
    message("Running permutation test for family: ", .x, "\n")
    result <- permTest(
      A = genomic_windows,
      B = .y,
      genome = genome_granges,
      randomize.function = randomize_ervs,
      evaluate.function   = numOverlaps,
      ntimes = n_perms
    )
    # Return a list containing the family name and the permTest result.
    list(family = .x, perm_result = result)
  }
)

# =============================================================================
# 7. Extract and Export Permutation Test Summary Results
# =============================================================================
# Extract summary statistics (observed overlaps, p-value, z-score, number of permutations, alternative) 
# and include the species name.
summary_df <- map_dfr(perm_results, function(res) {
  data.frame(
    Species      = species_name,
    Family       = res$family,
    Observed     = res$perm_result$numOverlaps$observed, 
    P_value      = res$perm_result$numOverlaps$pval, 
    Z_score      = res$perm_result$numOverlaps$zscore,
    Permutations = res$perm_result$numOverlaps$ntimes,
    Alternative  = res$perm_result$numOverlaps$alternative
  )
})

# Export the summary CSV for post-analysis.
write_csv(summary_df, file.path(args.csv_output_dir, paste0(species_name, ".csv")))

# =============================================================================
# 8. Export Permutation Test Plots (ggplot2) for Each Family
# =============================================================================
# Save all permutation test plots in a multipage PDF.
pdf(file.path(args.pdf_output_dir, (paste0(species_name, ".pdf"))), width = 15, height = 12)
for (res in perm_results) {
  # Extract permutation values and observed value
  perm_values    <- res$perm_result$numOverlaps$permuted
  observed_value <- res$perm_result$numOverlaps$observed
  
  # Extract additional info
  p_value     <- res$perm_result$numOverlaps$pval
  ntimes      <- res$perm_result$numOverlaps$ntimes
  alternative <- res$perm_result$numOverlaps$alternative
  
  # Create a data frame for the null distribution
  df <- data.frame(Overlaps = perm_values)
  
  # Build the ggplot histogram with a vertical line for the observed value.
  p <- ggplot(df, aes(x = Overlaps)) +
    geom_histogram(binwidth = 10, fill = "lightblue", color = "black") +
    geom_vline(xintercept = observed_value, color = "red", linetype = "dashed", size = 1) +
    labs(
      title = paste("Permutation Test -", species_name, "\nFamily:", res$family),
      subtitle = paste("Observed =", observed_value,
                       "| Permuted =", perm_values,
                       "| p-value =", p_value,
                       "| Permutations =", ntimes,
                       "| Alternative:", alternative),
      x = "Number of Overlaps",
      y = "Frequency"
    ) +
    theme_minimal() +
    scale_fill_nejm()
  
  print(p)
}
suppressMessages({dev.off()})

# =============================================================================
# 9. Integration Hotspot Analysis per Family (Window-Based)
# =============================================================================
# Set the number of permutations for the hotspot analysis.
n_perm_hotspots <- n_perms

# Initialize a list to store hotspot results for each family.
hotspots_list <- list()

# Loop over each family in the GRangesList (gff_by_family).
for (fam in names(gff_by_family)) {
  cat("Performing hotspot permutation analysis for family:", fam, "\n")
  
  # Get integration events for this family.
  events <- gff_by_family[[fam]]
  
  # Count observed integration events per genomic window.
  obs_counts <- countOverlaps(genomic_windows, events)
  
  # Create a matrix to store null counts for each window across n_perm_hotspots iterations.
  null_matrix <- matrix(NA, nrow = length(genomic_windows), ncol = n_perm_hotspots)
  
  # Perform randomization n_perm_hotspots times.
  for (i in 1:n_perm_hotspots) {
    randomized_events <- randomize_ervs(events, genome_granges)
    null_matrix[, i] <- countOverlaps(genomic_windows, randomized_events)
  }
  
  # For each window, compute an empirical p-value: the fraction of permutations with a count >= the observed count.
  pvals <- apply(null_matrix, 1, function(x) mean(x >= obs_counts))
  
  # Flag windows as hotspots if p-value < threshold.
  pvalue_threshold <- args.pvalue_threshold
  hotspot_windows <- genomic_windows[pvals < pvalue_threshold]
  
  # Store the results for this family.
  hotspots_list[[fam]] <- list(
    observed_counts = obs_counts,
    pvalues         = pvals,
    hotspots        = hotspot_windows
  )
}

# -----------------------------------------------------------------------------
# 10.  Export All Hotspot Regions as a Single GFF3 File
# -----------------------------------------------------------------------------

# Initialize an empty GRanges object
all_hotspots <- GRanges()

# Loop through each family in hotspots_list
for (fam in names(hotspots_list)) {
  hs <- hotspots_list[[fam]]$hotspots
  if (length(hs) > 0) {
    # Tag the hotspots with the family name
    mcols(hs)$family <- fam
    # Concatenate with previous results
    all_hotspots <- c(all_hotspots, hs)
  }
}

# Define the output file path
outfile <- file.path(args.hotspot_output_dir, paste0(species_name, ".gff3"))

# Export the combined hotspots as a single GFF3 file
rtracklayer::export(all_hotspots, outfile, format = "gff3")
