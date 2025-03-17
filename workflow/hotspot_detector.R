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
  library(ggdist)
  library(viridis)
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
print(paste0("Initializing permutation analysis for ", species_name, " ..."))

subject_genome <- readDNAStringSet(args.genome)

genome_granges <- GRanges(
  seqnames = names(subject_genome),
  ranges   = IRanges(start = 1, end = width(subject_genome))
)

# Clean names in case of chromosomes with description
chrom.names.filtered <- stringr::str_extract(names(subject_genome), "^[A-Za-z]+_?[0-9]+\\.[0-9]{1,2}")

seqinfo.vec <- setNames(width(subject_genome), chrom.names.filtered)

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
mask_match <- vmatchPattern(mask_pattern, subject_genome, max.mismatch = mask_mismatch)
genome_mask <- as(mask_match, "GRanges")

# =============================================================================
# 4. Define a Randomization Function for ERV Regions
# =============================================================================
randomize_ervs <- function(A, genome, ...) {
  randomizeRegions(A, 
                   genome = genome,
                   mask = genome_mask,
                   per.chromosome = TRUE)
}

# =============================================================================
# 5. Import and Split GFF3 Data by "family"
# =============================================================================
suppressMessages({
  subject_genome_gff <- rtracklayer::import(args.gff, format = "gff3")
})

if (!"family" %in% colnames(mcols(subject_genome_gff))) {
  stop("The GFF3 file does not contain a 'family' metadata field.")
}

subject_genome_gff_clean <- as(subject_genome_gff, "GRanges")
gff_by_family <- split(subject_genome_gff_clean, mcols(subject_genome_gff_clean)$family)

# =============================================================================
# 6. Perform Permutation Tests for Each Family (for Overall Overlap)
# =============================================================================
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
    list(family = .x, perm_result = result)
  }
)

# =============================================================================
# 7. Extract and Export Permutation Test Summary Results
# =============================================================================
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
summary_df$Adjusted_pvalue <- p.adjust(summary_df$P_value, method = "BH")

write_csv(summary_df, file.path(args.csv_output_dir, paste0(species_name, ".csv")))

# =============================================================================
# 8. Export Permutation Test Plots (Histogram) for Each Family
# =============================================================================
pdf(file.path(args.pdf_output_dir, paste0(species_name, "_histogram.pdf")), width = 15, height = 12)
for (res in perm_results) {
  perm_values    <- res$perm_result$numOverlaps$permuted   # Recompute for easier plotting within a single pdf
  observed_value <- res$perm_result$numOverlaps$observed
  ntimes         <- res$perm_result$numOverlaps$ntimes
  alternative    <- res$perm_result$numOverlaps$alternative
  p_value        <- res$perm_result$numOverlaps$pval
  adjusted_p     <- summary_df %>%
    filter(Family == res$family) %>%
    pull(Adjusted_pvalue)
  
  df <- data.frame(Overlaps = perm_values)
  perm_mean <- mean(perm_values)
  perm_sd   <- sd(perm_values)
  
  p_hist <- ggplot(df, aes(x = Overlaps)) +
    geom_histogram(binwidth = 10, fill = "lightblue", color = "black") +
    geom_vline(xintercept = observed_value, color = "red", linetype = "dashed", size = 1) +
    labs(
      title = paste("Permutation Test -", species_name, "\nFamily:", res$family),
      subtitle = paste(
        "Observed overlaps =", observed_value,
        "| Permuted overlaps (mean ± SD) =",
        sprintf("%.2f ± %.2f", perm_mean, perm_sd),
        "| Adjusted p-value (BH) =", sprintf("%.3g", adjusted_p),
        "| Permutations =", ntimes,
        "| Alternative:", alternative
      ),
      x = "Number of Overlaps",
      y = "Frequency"
    ) +
    theme_minimal() +
    scale_fill_nejm()
  
  print(p_hist)
}
dev.off()

# =============================================================================
# 9. Integration Hotspot Analysis per Family (Window-Based)
# =============================================================================
n_perm_hotspots <- n_perms
hotspots_list <- list()

for (fam in names(gff_by_family)) {
  cat("Performing hotspot empirical permutation analysis for family:", fam, "\n")
  
  events <- gff_by_family[[fam]]   # Factual observations
  obs_counts <- countOverlaps(genomic_windows, events)
  
  null_matrix <- matrix(NA, nrow = length(genomic_windows), ncol = n_perm_hotspots)
  for (i in seq_len(n_perm_hotspots)) {
    randomized_events <- randomize_ervs(events, genome_granges)
    null_matrix[, i] <- countOverlaps(genomic_windows, randomized_events)
  }
  
  pvals <- apply(null_matrix, 1, function(x) mean(x >= obs_counts))
  pvals_adjusted <- p.adjust(pvals, method = "BH")
  
  pvalue_threshold <- args.pvalue_threshold
  hotspot_windows <- genomic_windows[pvals_adjusted < pvalue_threshold]
  
  hotspots_list[[fam]] <- list(
    observed_counts = obs_counts,
    pvalues         = pvals_adjusted,
    hotspots        = hotspot_windows
  )
}

all_hotspots <- GRanges()
for (fam in names(hotspots_list)) {
  hs <- hotspots_list[[fam]]$hotspots
  if (length(hs) > 0) {
    mcols(hs)$family <- fam
    all_hotspots <- c(all_hotspots, hs)
  }
}

outfile <- file.path(args.hotspot_output_dir, paste0(species_name, ".gff3"))
rtracklayer::export(all_hotspots, outfile, format = "gff3")

# =============================================================================
# A) EMPIRICAL DENSITY PLOTS (Permutation Distribution)
# =============================================================================
density_pdf_path <- file.path(args.pdf_output_dir, paste0(species_name, "_density.pdf"))
pdf(density_pdf_path, width = 15, height = 12)
for (res in perm_results) {
  family_name    <- res$family
  perm_values    <- res$perm_result$numOverlaps$permuted
  observed_value <- res$perm_result$numOverlaps$observed
  perm_mean      <- mean(perm_values)
  perm_sd        <- sd(perm_values)
  
  adjusted_pval <- summary_df %>%
    filter(Family == family_name) %>%
    pull(Adjusted_pvalue)
  
  df <- data.frame(Overlaps = perm_values)
  
  p_density <- ggplot(df, aes(x = Overlaps)) +
    geom_density(fill = "lightblue", alpha = 0.5) +
    geom_vline(xintercept = observed_value, color = "red", linetype = "dashed", size = 1) +
    labs(
      title = paste0("Empirical Density of Permuted Overlaps - ", species_name),
      subtitle = paste(
        "Family:", family_name,
        "| Observed =", observed_value,
        "| Perm. Mean ± SD =", sprintf("%.2f ± %.2f", perm_mean, perm_sd),
        "| Adj. p-value =", sprintf("%.3g", adjusted_pval)
      ),
      x = "Number of Overlaps",
      y = "Density"
    ) +
    theme_minimal()
  
  print(p_density)
}
dev.off()


# =============================================================================
# B) HEATMAP OF OBSERVED COUNTS
# =============================================================================
heatmap_pdf_path <- file.path(args.pdf_output_dir, paste0(species_name, "_heatmap.pdf"))
pdf(heatmap_pdf_path, width = 15, height = 12)

for (fam in names(hotspots_list)) {
  obs_counts <- hotspots_list[[fam]]$observed_counts
  fam_df <- data.frame(
    seqnames     = as.character(seqnames(genomic_windows)),
    start        = start(genomic_windows),
    OverlapCount = obs_counts
  )
  
  p_heatmap <- ggplot(fam_df, aes(x = seqnames, y = start, fill = OverlapCount)) +
    geom_tile() +
    scale_fill_bs5("pink") +
    labs(
      title = paste("Heatmap of Observed Overlap Counts -", species_name),
      subtitle = paste("Family:", fam),
      x = "Chromosome",
      y = "Window Start (bp)",
      fill = "Overlap\nCount"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=90, hjust=1))
  
  print(p_heatmap)
}
dev.off()
