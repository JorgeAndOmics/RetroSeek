# ==========================
#   CIRCULAR GENOME PLOT
# ==========================
# This script reads genomic data (FASTA and GFF files), processes it into 
# a GenomicRanges format, and visualizes it using ggbio in a circular plot.
# The final output is saved as a high-resolution image.

# ------------------------------
# 1. LOAD REQUIRED LIBRARIES
# ------------------------------

# Suppress warnings
options(warn = -1)

suppressMessages({
  library(ggbio)          # For genomic visualization (circular plots)
  library(rtracklayer)    # For importing GFF annotation files
  library(Biostrings)     # For handling DNA sequences
  library(GenomicRanges)  # For managing genomic intervals
  library(plyranges)      # For genomic data manipulation
  library(ggsci)          # For scientific color palettes
  library(ggnewscale)     # For adding new color scales
})

# ------------------------------
# 2. PARSE COMMAND-LINE ARGUMENTS
# ------------------------------
# The script expects 4 input arguments:
#   1. FASTA file (genome sequence)
#   2. GFF custom annotation file
#   3. GFF LTRDigest output file
#   4. Bitscore threshold for plotting
#   5. Output file path for saving the plot, without extension


args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Please provide exactly five arguments: fasta file, gff custom file, and gff LTRDigest file, bitscore threshold and output file")
}

# Assign arguments to variables
args.fasta_file <- args[1]
args.gff_custom_file <- args[2]
args.gff_ltrdigest_file <- args[3]
args.bitscore_threshold <- as.numeric(args[4])
args.output_file <- args[5]

# Print progress message
print(paste0("Processing circle plots for ", basename(args.gff_custom_file), "..."))

# ------------------------------
# 3. LOAD FASTA FILE (GENOME SEQUENCES)
# ------------------------------
# Reads the FASTA genome file and extracts chromosome names and lengths

genome <- readDNAStringSet(args.fasta_file)

# Extract chromosome names (Read accession via regex in case they contain extensive information instead of just the name)
chr_names <- stringr::str_extract(names(genome), "^[A-Za-z]+_?[0-9]+\\.[0-9]{1,2}")
chr_lengths <- width(genome)  # Get sequence lengths

# Create a GenomicRanges object for chromosomes
chromosomes <- GRanges(seqnames = chr_names,
                       ranges   = IRanges(start = 1, end = chr_lengths))

# Create a Seqinfo object and assign it to chromosomes
my_seqinfo <- Seqinfo(seqnames   = chr_names,
                      seqlengths = chr_lengths)
seqinfo(chromosomes) <- my_seqinfo  # Attach metadata

# ------------------------------
# 4. IMPORT GFF FILES (ANNOTATIONS)
# ------------------------------
# Reads two GFF3 annotation files:
# - gff_custom: Custom feature annotations
# - gff_ltrdigest: LTR retrotransposon predictions

suppressMessages({
  gff_custom <- rtracklayer::import(args.gff_custom_file)
  gff_ltrdigest <- rtracklayer::import(args.gff_ltrdigest_file)
}) 

# Subset LTRDigest results to retain only overlapping features from the custom GFF
gff_ltrdigest <- subsetByOverlaps(gff_ltrdigest, gff_custom)

# ------------------------------
# 5. FILTER CHROMOSOMES IN ALL DATASETS
# ------------------------------
# Ensure that all datasets only contain common chromosome names

common_levels <- intersect(seqlevels(chromosomes), seqlevels(gff_custom))
common_levels <- intersect(common_levels, seqlevels(gff_ltrdigest))

chromosomes   <- keepSeqlevels(chromosomes,  common_levels, pruning.mode = "coarse")
gff_custom    <- keepSeqlevels(gff_custom,   common_levels, pruning.mode = "coarse")
gff_ltrdigest <- keepSeqlevels(gff_ltrdigest, common_levels, pruning.mode = "coarse")

# ------------------------------
# 6. DEFINE FEATURE TYPE MAPPING
# ------------------------------
# This assigns numerical values to different genomic features for plotting

type_mapping <- c(
  "RR_tract"                = 1,
  "protein_match"           = 2,
  "long_terminal_repeat"    = 3,
  "target_site_duplication" = 4,
  "repeat_region"           = 5,
  "LTR_retrotransposon"     = 6  # Outer ring (most significant feature)
)

# ------------------------------
# 7. CLEAN AND TRANSFORM DATA
# ------------------------------
# Convert certain columns from factor to numeric where necessary

gff_custom <- gff_custom %>%
  mutate(mean_bitscore = as.numeric(mean_bitscore),
         mean_identity = as.numeric(mean_identity))

gff_ltrdigest <- gff_ltrdigest %>%
  mutate(
    ltr_similarity = as.numeric(ltr_similarity),
    type_num = type_mapping[as.character(type)]
  )

# Ensure that sequence information matches across datasets
seqinfo(gff_custom) <- seqinfo(chromosomes)
seqinfo(gff_ltrdigest) <- seqinfo(chromosomes)

# ------------------------------
# 8. PLOT CIRCULAR GENOME MAP
# ------------------------------
p <- ggplot() +
  
  # (A) Chromosome ring: Represents the genome
  layout_circle(chromosomes, geom = "ideo", fill = "white", color = "grey", alpha = 0.5, radius = 9) +
  
  # (B) Feature ring: Displays genomic features colored by probe type
  layout_circle(
    gff_custom,
    geom = "rect",
    aes(
      color  = factor(probe),
      fill   = factor(probe),              
      alpha  = as.numeric(mean_bitscore)  
    ),
    radius = 9
  ) +
  
  # (C) Outermost ring: Plots probe identities if mean_bitscore > threshold
  layout_circle(
    subset(gff_custom, as.numeric(mean_bitscore) > args.bitscore_threshold),
    aes(y = mean_bitscore, color = probe, alpha = mean_identity),
    geom = "point",
    radius = 14,
    size = 2
  ) + 
  
  # (D) Chromosome labels
  layout_circle(
    data = chromosomes,
    geom = "text",
    aes(label = seqnames, size = width),
    radius    = 8.5,
    trackWidth= 1,
    hjust     = 0.5,
    vjust     = 0.5,
    check_overlap = TRUE
  ) +
  
  # Color scale for probe features
  scale_colour_futurama(aesthetics = c("colour", "fill"),
                        guide = guide_legend(title = "Probe", order = 1)) + 
  
  # ------------------------------
# 8.1 RESET COLOR SCALES FOR STRAND INFO
# ------------------------------

new_scale_color() +
  new_scale_fill() +
  
  # (E) LTR Digest Ring: Displays LTR annotations (outermost track)
  layout_circle(gff_ltrdigest, geom = "point",
                aes(y = type_num, shape = type, colour= strand, fill = strand, size = width),
                radius = 3, trackWidth = 5, alpha = 0.4) +
  
  # Strand-specific color scheme
  scale_colour_manual(
    aesthetics = c("colour", "fill"),
    name = "Strand",
    values = c("+" = "hotpink", "-" = "deepskyblue"),
    guide = guide_legend(order = 2)
  ) +
  
  # Additional styling settings
  scale_alpha(range = c(0.3, 1), guide = "none") +
  scale_shape_manual(
    name   = "Feature",
    values = c(
      "LTR_retrotransposon"  = 22,
      "protein_match"        = 20,
      "target_site_duplication" = 23,
      "long_terminal_repeat" = 24,
      "repeat_region"        = 25,
      "RR_tract"             = 3
    ),
    guide = guide_legend(order = 3)
  ) +
  
  scale_size_continuous(range = c(1, 3), guide = "none") +
  
  # Remove background and set properties
  theme_void() +
  theme(aspect.ratio = 1,
        text = element_text(face = "bold"),
        legend.position = "right"
  )

# ------------------------------
# 9. SAVE PLOT AS IMAGE FILE
# ------------------------------

ggplot2::ggsave(filename=paste0(args.output_file, ".png"), plot = p, width = 20, height = 20, dpi = 500)
ggplot2::ggsave(filename=paste0(args.output_file, ".pdf"), plot = p, width = 20, height = 20, dpi = 500)