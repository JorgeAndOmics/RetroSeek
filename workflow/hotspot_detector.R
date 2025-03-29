# =============================================================================
# Suppress Warnings and Load Libraries
# =============================================================================
options(warn = -1)
suppressMessages({
  library(argparse)       # Argument parsing
  library(Biostrings)     # DNA sequence utilities
  library(rtracklayer)    # GFF I/O
  library(GenomicRanges)  # Genomic interval operations
  library(regioneR)       # Permutation testing for genomic intervals
  library(tidyverse)      # Data manipulation and visualization
  library(purrr)          # Functional programming tools
  library(yaml)           # YAML config parsing
})

# =============================================================================
# 1. Parse Command-Line Arguments
# =============================================================================
parser <- ArgumentParser(description = "ERV Hotspot Permutation Analysis")

parser$add_argument("--genome", required=TRUE, help="Path to genome FASTA file")
parser$add_argument("--gff", required=TRUE, help="Path to GFF3 with ERV annotations")
parser$add_argument("--config", required=TRUE, help="YAML config file")
parser$add_argument("--csv_output_dir", required=TRUE, help="Directory to write CSV output")
parser$add_argument("--pdf_output_dir", required=TRUE, help="Directory to write PDF output")
parser$add_argument("--hotspot_output_dir", required=TRUE, help="Directory to write GFF3 hotspot output")

args <- parser$parse_args()

args.genome <- file.path(args$genome)
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
# 5. Import and Split GFF3 by Family
# =============================================================================
subject_genome_gff <- rtracklayer::import(args.gff, format = "gff3")

if (!"family" %in% colnames(mcols(subject_genome_gff))) {
  stop("The GFF3 file does not contain a 'family' metadata field.")
}

subject_genome_gff_clean <- as(subject_genome_gff, "GRanges")

gff_list <- if (config$parameters$hotspot_group_split) {
  split(subject_genome_gff_clean, mcols(subject_genome_gff_clean)$family)
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
    list(family = .x, perm_result = result)
  }
)

# =============================================================================
# 7. Summarize Permutation Results
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
