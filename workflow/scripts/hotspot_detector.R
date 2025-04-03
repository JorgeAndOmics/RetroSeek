# =============================================================================
# Suppress Warnings and Load Libraries
# =============================================================================
options(warn = -1)
suppressMessages({
  library(argparse)       # Argument parsing
  library(Biostrings)     # DNA sequence utilities
  library(rtracklayer)    # GFF I/O
  library(GenomicRanges)  # Genomic interval operations
  library(plyranges)      # Genomic data manipulation
  library(regioneR)       # Permutation testing for genomic intervals
  library(tidyverse)      # Data manipulation and visualization
  library(purrr)          # Functional programming tools
  library(yaml)           # YAML config parsing
  library(ggsci)          # Scientific color palettes
})

# =============================================================================
# 1. Parse Command-Line Arguments
# =============================================================================
parser <- ArgumentParser(description = "ERV Hotspot Permutation Analysis")

parser$add_argument("--fasta", required=TRUE, help="Path to genome FASTA file")
parser$add_argument("--gff", required=TRUE, help="Path to GFF3 with ERV annotations")
parser$add_argument("--config", required=TRUE, help="YAML config file")
parser$add_argument("--csv_output_dir", required=TRUE, help="Directory to write CSV output")
parser$add_argument("--pdf_output_dir", required=TRUE, help="Directory to write PDF output")
parser$add_argument("--hotspot_output_dir", required=TRUE, help="Directory to write GFF3 hotspot output")

args <- parser$parse_args()

args.genome <- file.path(args$fasta)
args.gff <- file.path(args$gff)
args.yaml <- file.path(args$config)
args.csv_output_dir <- file.path(args$csv_output_dir)
args.pdf_output_dir <- file.path(args$pdf_output_dir)
args.hotspot_output_dir <- file.path(args$hotspot_output_dir)

config <- yaml::read_yaml(args.yaml)

species <- tools::file_path_sans_ext(basename(args.genome))
species_name <- ifelse(!is.na(config$species[species]), config$species[species], species)

# =============================================================================
# 2. Load Genome and Create Genome Boundaries
# =============================================================================
print(paste0("Initializing permutation analysis for ", species_name, "..."))

subject_genome <- readDNAStringSet(args.genome)

genome_granges <- GRanges(
  seqnames = names(subject_genome),
  ranges   = IRanges(start = 1, end = width(subject_genome))
)

chrom.names.filtered <- stringr::str_extract(names(subject_genome), "^[A-Za-z]+_?[0-9]+\\.[0-9]{1,2}")

seqlevels(genome_granges) <- chrom.names.filtered
seqnames(genome_granges) <- chrom.names.filtered

seqinfo.vec <- setNames(width(subject_genome), chrom.names.filtered)

# =============================================================================
# 3. Create Genomic Windows and Mask Regions
# =============================================================================
print(paste0("Segmenting ", species_name, " genome in tiles and masking ..."))

window_size <- as.integer(config$parameters$hotspot_window_size)
genomic_windows <- tileGenome(
  seqlengths = seqinfo.vec,
  tilewidth = window_size,
  cut.last.tile.in.chrom = TRUE
)

mask_size <- as.integer(config$parameters$hotspot_mask_size)
mask_mismatch <- as.integer(config$parameters$hotspot_mask_mismatch)
mask_pattern <- DNAString(paste(rep("N", mask_size), collapse = ""))
mask_match <- vmatchPattern(mask_pattern, subject_genome, max.mismatch = mask_mismatch)
genome_mask <- as(mask_match, "GRanges")

# =============================================================================
# 4. Define Randomization Function
# =============================================================================
randomize_ervs <- function(A, genome, ...) {
  randomizeRegions(A, genome = genome, mask = genome_mask, per.chromosome = TRUE)
}

# =============================================================================
# 5. Import and Split GFF3 by Label
# =============================================================================
subject_genome_gff <- rtracklayer::import(args.gff, format = "gff3")

if (!"label" %in% colnames(mcols(subject_genome_gff))) {
  stop("The GFF3 file does not contain a 'label' metadata field.")
}

subject_genome_gff_clean <- as(subject_genome_gff, "GRanges")

gff_list <- if (config$parameters$hotspot_group_split) {
  split(subject_genome_gff_clean, mcols(subject_genome_gff_clean)$label)
} else {
  list("Ungrouped" = subject_genome_gff_clean)
}

# =============================================================================
# 6. Permutation Tests
# =============================================================================
n_perms <- config$parameters$hotspot_permutations

perm_results <- map2(
  names(gff_list), gff_list,
  ~{
    message("Running permutation test: ", .x)
    result <- permTest(
      A = .y,
      B = genomic_windows,
      genome = genome_granges,
      randomize.function = randomize_ervs,
      evaluate.function = numOverlaps,
      ntimes = n_perms
    )
    list(label = .x, perm_result = result)
  }
)

# =============================================================================
# 7. Summarize Permutation Results
# =============================================================================
summary_df <- map_dfr(perm_results, function(res) {
  data.frame(
    Species      = species_name,
    Label       = res$label,
    Observed     = res$perm_result$numOverlaps$observed,
    P_value      = res$perm_result$numOverlaps$pval,
    Z_score      = res$perm_result$numOverlaps$zscore,
    Permutations = res$perm_result$numOverlaps$ntimes,
    Alternative  = res$perm_result$numOverlaps$alternative
  )
})
summary_df$Adjusted_pvalue <- p.adjust(summary_df$P_value, method = "BH")

write_csv(summary_df, file.path(args.csv_output_dir, paste0(species, ".csv")))

# =============================================================================
# 8. Export Permutation Test Plots (Histogram) for Each Label
# =============================================================================
pdf(file.path(args.pdf_output_dir, paste0(species, "_histogram.pdf")), width = 15, height = 12)
for (res in perm_results) {
  perm_values    <- res$perm_result$numOverlaps$permuted   # Recompute for easier plotting within a single pdf
  observed_value <- res$perm_result$numOverlaps$observed
  ntimes         <- res$perm_result$numOverlaps$ntimes
  alternative    <- res$perm_result$numOverlaps$alternative
  p_value        <- res$perm_result$numOverlaps$pval
  adjusted_p     <- summary_df %>%
    filter(Label == res$label) %>%
    pull(Adjusted_pvalue)
  
  df <- data.frame(Overlaps = perm_values)
  perm_mean <- mean(perm_values)
  perm_sd   <- sd(perm_values)
  
  p_hist <- ggplot(df, aes(x = Overlaps)) +
    geom_histogram(binwidth = 10, fill = "lightblue", color = "black") +
    geom_vline(xintercept = observed_value, color = "red", linetype = "dashed", size = 1) +
    labs(
      title = paste("Permutation Test -", species_name, "\nLabel:", res$label),
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
# 9. Integration Hotspot Analysis per Label (Window-Based)
# =============================================================================
n_perm_hotspots <- n_perms
hotspots_list <- list()

for (fam in names(gff_list)) {
  cat("Performing hotspot empirical permutation analysis for label:", fam, "\n")
  
  events <- gff_list[[fam]]   # Factual observations
  obs_counts <- countOverlaps(genomic_windows, events)
  
  null_matrix <- matrix(NA, nrow = length(genomic_windows), ncol = n_perm_hotspots)
  for (i in seq_len(n_perm_hotspots)) {
    randomized_events <- randomize_ervs(events, genome_granges)
    null_matrix[, i] <- countOverlaps(genomic_windows, randomized_events)
  }
  
  pvals <- apply(null_matrix, 1, function(x) mean(x >= obs_counts))
  pvals_adjusted <- p.adjust(pvals, method = "BH")
  
  pvalue_threshold <- config$parameters$hotspot_pvalue_threshold
  
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
    mcols(hs)$label <- fam
    all_hotspots <- c(all_hotspots, hs)
  }
}

outfile <- file.path(args.hotspot_output_dir, paste0(species, ".gff3"))
rtracklayer::export(all_hotspots, outfile, format = "gff3")

# =============================================================================
# A) EMPIRICAL DENSITY PLOTS (Permutation Distribution)
# =============================================================================
density_pdf_path <- file.path(args.pdf_output_dir, paste0(species, "_density.pdf"))
pdf(density_pdf_path, width = 15, height = 12)
for (res in perm_results) {
  label_name    <- res$label
  perm_values    <- res$perm_result$numOverlaps$permuted
  observed_value <- res$perm_result$numOverlaps$observed
  perm_mean      <- mean(perm_values)
  perm_sd        <- sd(perm_values)
  
  adjusted_pval <- summary_df %>%
    filter(Label == label_name) %>%
    pull(Adjusted_pvalue)
  
  df <- data.frame(Overlaps = perm_values)
  
  p_density <- ggplot(df, aes(x = Overlaps)) +
    geom_density(fill = "lightblue", alpha = 0.5) +
    geom_vline(xintercept = observed_value, color = "red", linetype = "dashed", size = 1) +
    labs(
      title = paste0("Empirical Density of Permuted Overlaps - ", species_name),
      subtitle = paste(
        "Label:", label_name,
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
heatmap_pdf_path <- file.path(args.pdf_output_dir, paste0(species, "_heatmap.pdf"))
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
      subtitle = paste("Label:", fam),
      x = "Chromosome",
      y = "Window Start (bp)",
      fill = "Overlap\nCount"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=90, hjust=1))
  
  print(p_heatmap)
}
dev.off()
