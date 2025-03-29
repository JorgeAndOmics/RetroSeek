# ==========================
#   CIRCULAR GENOME PLOT
# ==========================
# This script reads genomic data (FASTA and GFF files), processes it into 
# a GenomicRanges format, and visualizes it using ggbio in a circular plot.
# The final output is saved as a high-resolution image.

# ------------------------------
# 1. LOAD REQUIRED LIBRARIES
# ------------------------------

options(warn = -1)  # Suppress warnings

suppressMessages({
  library(argparse)        # For argument parsing
  library(ggbio)           # For genomic visualization (circular plots)
  library(rtracklayer)     # For importing GFF annotation files
  library(Biostrings)      # For handling DNA sequences
  library(GenomicRanges)   # For managing genomic intervals
  library(plyranges)       # For genomic data manipulation
  library(ggsci)           # For scientific color palettes
  library(ggnewscale)      # For adding new color scales
  library(stringr)         # For regex-based string extraction
})

# ------------------------------
# 2. PARSE COMMAND-LINE ARGUMENTS
# ------------------------------

parser <- ArgumentParser(description = "Create circular genome plot from FASTA and GFF annotations")

parser$add_argument("--fasta", required = TRUE, help = "FASTA file containing the genome")
parser$add_argument("--gff_custom", required = TRUE, help = "Custom GFF annotation file")
parser$add_argument("--gff_ltrdigest", required = TRUE, help = "LTRdigest GFF annotation file")
parser$add_argument("--threshold", required = TRUE, type = "double", help = "Bitscore threshold for outer ring")
parser$add_argument("--output", required = TRUE, help = "Output prefix path (no extension)")

args <- parser$parse_args()

args.fasta_file        <- args$fasta
args.gff_custom_file   <- args$gff_custom
args.gff_ltrdigest_file<- args$gff_ltrdigest
args.bitscore_threshold<- args$threshold
args.output_file       <- args$output

# Inform user
message("Processing circle plots for: ", basename(args.gff_custom_file))

# ------------------------------
# 3. LOAD FASTA FILE (GENOME SEQUENCES)
# ------------------------------
genome <- readDNAStringSet(args.fasta_file)

# Extract chromosome names and lengths
chr_names   <- str_extract(names(genome), "^[A-Za-z]+_?[0-9]+\\.[0-9]{1,2}")
chr_lengths <- width(genome)

chromosomes <- GRanges(seqnames = chr_names, ranges = IRanges(start = 1, end = chr_lengths))
seqinfo(chromosomes) <- Seqinfo(seqnames = chr_names, seqlengths = chr_lengths)

# ------------------------------
# 4. IMPORT GFF FILES (ANNOTATIONS)
# ------------------------------
gff_custom    <- rtracklayer::import(args.gff_custom_file)
gff_ltrdigest <- rtracklayer::import(args.gff_ltrdigest_file)

# Filter LTRdigest to only overlapping entries with custom annotations
gff_ltrdigest <- subsetByOverlaps(gff_ltrdigest, gff_custom)

# ------------------------------
# 5. FILTER CHROMOSOMES IN ALL DATASETS
# ------------------------------
common_levels <- Reduce(intersect, list(seqlevels(chromosomes), seqlevels(gff_custom), seqlevels(gff_ltrdigest)))

chromosomes   <- keepSeqlevels(chromosomes,   common_levels, pruning.mode = "coarse")
gff_custom    <- keepSeqlevels(gff_custom,    common_levels, pruning.mode = "coarse")
gff_ltrdigest <- keepSeqlevels(gff_ltrdigest, common_levels, pruning.mode = "coarse")

# ------------------------------
# 6. DEFINE FEATURE TYPE MAPPING
# ------------------------------
type_mapping <- c(
  "RR_tract"                = 1,
  "protein_match"           = 2,
  "long_terminal_repeat"    = 3,
  "target_site_duplication" = 4,
  "repeat_region"           = 5,
  "LTR_retrotransposon"     = 6
)

# ------------------------------
# 7. CLEAN AND TRANSFORM DATA
# ------------------------------
gff_custom <- gff_custom %>%
  mutate(
    mean_bitscore = as.numeric(mean_bitscore),
    mean_identity = as.numeric(mean_identity)
  )

gff_ltrdigest <- gff_ltrdigest %>%
  mutate(
    ltr_similarity = as.numeric(ltr_similarity),
    type_num = type_mapping[as.character(type)]
  )

seqinfo(gff_custom)    <- seqinfo(chromosomes)
seqinfo(gff_ltrdigest) <- seqinfo(chromosomes)

# ------------------------------
# 8. PLOT CIRCULAR GENOME MAP
# ------------------------------
p <- ggplot() +
  layout_circle(chromosomes, geom = "ideo", fill = "white", color = "grey", alpha = 0.5, radius = 9) +
  
  layout_circle(
    gff_custom,
    geom = "rect",
    aes(
      color = factor(probe),
      fill  = factor(probe),
      alpha = as.numeric(mean_bitscore)
    ),
    radius = 9
  ) +
  
  layout_circle(
    subset(gff_custom, as.numeric(mean_bitscore) > args.bitscore_threshold),
    aes(y = mean_bitscore, color = probe, alpha = mean_identity),
    geom = "point",
    radius = 14,
    size = 2
  ) +
  
  layout_circle(
    data = chromosomes,
    geom = "text",
    aes(label = seqnames, size = width),
    radius = 8.5,
    trackWidth = 1,
    hjust = 0.5,
    vjust = 0.5,
    check_overlap = TRUE
  ) +
  
  scale_colour_futurama(aesthetics = c("colour", "fill"), guide = guide_legend(title = "Probe", order = 1)) +
  
  new_scale_color() +
  new_scale_fill() +
  
  layout_circle(
    gff_ltrdigest,
    geom = "point",
    aes(y = type_num, shape = type, colour = strand, fill = strand, size = width),
    radius = 3,
    trackWidth = 5,
    alpha = 0.4
  ) +
  
  scale_colour_manual(
    aesthetics = c("colour", "fill"),
    name = "Strand",
    values = c("+" = "hotpink", "-" = "deepskyblue"),
    guide = guide_legend(order = 2)
  ) +
  
  scale_alpha(range = c(0.3, 1), guide = "none") +
  
  scale_shape_manual(
    name = "Feature",
    values = c(
      "LTR_retrotransposon"     = 22,
      "protein_match"           = 20,
      "target_site_duplication" = 23,
      "long_terminal_repeat"    = 24,
      "repeat_region"           = 25,
      "RR_tract"                = 3
    ),
    guide = guide_legend(order = 3)
  ) +
  
  scale_size_continuous(range = c(1, 3), guide = "none") +
  
  theme_void() +
  theme(
    aspect.ratio = 1,
    text = element_text(face = "bold"),
    legend.position = "right"
  )

# ------------------------------
# 9. SAVE PLOT AS IMAGE FILE
# ------------------------------
ggplot2::ggsave(filename = paste0(args.output, ".png"), plot = p, width = 20, height = 20, dpi = 500)
ggplot2::ggsave(filename = paste0(args.output, ".pdf"), plot = p, width = 20, height = 20, dpi = 500)
