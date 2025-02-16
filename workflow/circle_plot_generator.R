options(warn = -1)

suppressMessages({
library(ggbio)
library(rtracklayer)
library(Biostrings)
library(GenomicRanges)
library(plyranges)
library(ggsci)}
)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Please provide exactly three arguments: fasta file, gff custom file, and gff LTRDigest file and output file")
}

args.fasta_file = args[1]
args.gff_custom_file = args[2]
args.gff_ltrdigest_file = args[3]
args.output_file = args[4]

print(paste0("Processing circle plots for ", basename(args.gff_custom_file), "..."))

# Fasta file reading
genome     <- readDNAStringSet(args.fasta_file)

chr_names   <- names(genome)
chr_lengths <- width(genome)

# Create GRanges for chromosomes
chromosomes <- GRanges(seqnames = chr_names,
                       ranges   = IRanges(start = 1, end = chr_lengths))

# Create a Seqinfo object and assign it to chromosomes
my_seqinfo <- Seqinfo(seqnames   = chr_names,
                      seqlengths = chr_lengths)
seqinfo(chromosomes) <- my_seqinfo

# LTRDigest imports
suppressMessages({
  gff_custom      <- rtracklayer::import(args.gff_custom_file)
  gff_ltrdigest <- rtracklayer::import(args.gff_ltrdigest_file)
}) 

# Subset hits
gff_ltrdigest <- subsetByOverlaps(gff_ltrdigest, gff_custom)


# Keep only the common seqlevels in both objects
common_levels <- intersect(seqlevels(chromosomes), seqlevels(gff_custom))
common_levels <- intersect(common_levels, seqlevels(gff_ltrdigest))
chromosomes   <- keepSeqlevels(chromosomes,  common_levels, pruning.mode = "coarse")
gff_custom    <- keepSeqlevels(gff_custom,   common_levels, pruning.mode = "coarse")
gff_ltrdigest <- keepSeqlevels(gff_ltrdigest, common_levels, pruning.mode = "coarse")
type_mapping <- c(
  "RR_tract"              = 1,
  "protein_match"         = 2,
  "long_terminal_repeat"  = 3,
  "target_site_duplication" = 4,
  "repeat_region"         = 5,
  "LTR_retrotransposon"   = 6  # largest => outer ring
)


# Clean data formats
gff_custom <- gff_custom %>%
  mutate(mean_bitscore = as.numeric(mean_bitscore),
         mean_identity = as.numeric(mean_identity))

gff_ltrdigest <- gff_ltrdigest %>%
  mutate(
    ltr_similarity = as.numeric(ltr_similarity),
    type_num = type_mapping[as.character(type)]
  )


# Assign matching seqinfo
seqinfo(gff_custom) <- seqinfo(chromosomes)
seqinfo(gff_ltrdigest) <- seqinfo(chromosomes)



# Plot
p <- ggplot() +
  # (A) Chromosome ring
  layout_circle(chromosomes, geom = "ideo", fill = "white", color = "grey", alpha = 0.5, radius = 9) +
  
  
  # (B) Main feature ring with dynamic color by probe
  layout_circle(
    gff_custom,
    geom = "rect",
    aes(
      color  = factor(probe),
      fill = factor(probe),              
      alpha = as.numeric(mean_bitscore)  
    ),
    radius = 9
  ) +
  
  
  # (C) A separate ring for probe and identity
  layout_circle(
    subset(gff_custom, as.numeric(mean_bitscore) > 200),
    aes(y = mean_bitscore, color = probe, alpha = mean_identity),
    geom = "point",
    radius = 14,     # outward ring
    size = 2
  ) + 
  
 # LTRdigest  info
  layout_circle(gff_ltrdigest, geom = "point", aes(y = type_num, shape = type, color = strand, fill = strand),
         radius = 3, trackWidth = 5) +
  
  # --- Scales ---
  # Combine colour and fill into one discrete scale
  
  scale_colour_futurama(aesthetics = c("colour", "fill"),
                        guide = guide_legend(title = "Feature", order = 1)) + 
  
  
  # Remove legends for alpha and size
  scale_alpha(range = c(0.3, 1), guide = "none") +
  
  scale_size_continuous(range = c(1, 3), guide = "none") +
  
  scale_shape_manual(
    name   = "Sequence",
    values = c(
      "LTR_retrotransposon"  = 22,
      "protein_match"        = 20,
      "target_site_duplication" = 23,
      "long_terminal_repeat" = 24,
      "repeat_region"        = 25,
      "RR_tract"             = 3
    ),
    breaks = c(
      "repeat_region",
      "RR_tract",
      "target_site_duplication",
      "long_terminal_repeat",
      "protein_match",
      "LTR_retrotransposon"
    ),
    guide = guide_legend(order = 2)
  ) +
  
  theme_void() +
  theme(aspect.ratio = 1,
        text = element_text(face = "bold"),
        legend.position = "right")

p


# Save plot
ggplot2::ggsave(filename=args.output_file, plot = p, width = 20, height = 20, dpi = 300)
