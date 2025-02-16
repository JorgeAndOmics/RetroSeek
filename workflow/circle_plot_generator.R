# If needed, install packages:
# BiocManager::install(c("ggbio", "rtracklayer", "Biostrings"))

library(ggbio)
library(rtracklayer)
library(Biostrings)
library(GenomicRanges)
library(plyranges)
library(ggrepel)
library(ggsci)

# ------------------------------------------------
# 1) Read genome and create a GRanges "ideogram"
# ------------------------------------------------

futurama_unlimited_palette <- function(input_colour_number = 12, output_colour_number) {
  planet_express <- pal_futurama("planetexpress")(input_colour_number)
  output_colour <- colorRampPalette(planet_express)(output_colour_number)
  return (output_colour)
}

fasta_file <- "V:/databases/local/blast_dbs/species/Hipposideros_larvatus.fa"
genome     <- readDNAStringSet(fasta_file)

chr_names   <- names(genome)
chr_lengths <- width(genome)

# Create GRanges for chromosomes
chromosomes <- GRanges(seqnames = chr_names,
                       ranges   = IRanges(start = 1, end = chr_lengths))

# Create a Seqinfo object and assign it to chromosomes
my_seqinfo <- Seqinfo(seqnames   = chr_names,
                      seqlengths = chr_lengths)
seqinfo(chromosomes) <- my_seqinfo

# ------------------------------------------------
# 2) Import your GFF3 annotations
# ------------------------------------------------
gff_custom_file <- "C:/Users/Lympha/Documents/Repositories/enERVate/results/tracks/validated/Hipposideros_larvatus_main.gff3"
gff_custom      <- import(gff_custom_file)

gff_ltr_digest_file <- "C:/Users/Lympha/Documents/Repositories/enERVate/results/ltrdigest/Hipposideros_larvatus.gff3"
gff_ltrdigest <- rtracklayer::import(gff_ltr_digest_file)
gff_ltrdigest <- subsetByOverlaps(gff_ltrdigest, gff_custom) %>% filter(type == "LTR_retrotransposon")

# Keep only the common seqlevels in both objects
common_levels <- intersect(seqlevels(chromosomes), seqlevels(gff_custom))
common_levels <- intersect(common_levels, seqlevels(gff_ltrdigest))
chromosomes   <- keepSeqlevels(chromosomes,  common_levels, pruning.mode = "coarse")
gff_custom    <- keepSeqlevels(gff_custom,   common_levels, pruning.mode = "coarse")
gff_ltrdigest <- keepSeqlevels(gff_ltrdigest, common_levels, pruning.mode = "coarse")

# Assign matching seqinfo
seqinfo(gff_custom) <- seqinfo(chromosomes)
seqinfo(gff_ltrdigest) <- seqinfo(chromosomes)

# ------------------------------------------------
# 3) Plot everything in one concentric circle
# ------------------------------------------------
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
    aes(color = probe, alpha = as.numeric(mean_identity)),
    geom = "point",
    radius = 12,     # outward ring
    size = 2
  ) + 
  
  # inner ring with chromosome names
  layout_circle(chromosomes, geom = "text", aes(label = seqnames, size = width(chromosomes)), radius = 6.5) +
  theme_void() +
  theme(aspect.ratio = 1,
        text = element_text(face = "bold"),
        legend.position = "None") +
  
  
  # (B) LTRDigest feature ring with dynamic color by probe
  layout_circle(
    gff_ltrdigest,
    geom = "rect",
    aes(
      color  = factor(strand),
      fill = factor(strand),              
      alpha = as.numeric(ltr_similarity),
    ),
    radius = 3
  ) +
  # --- Scales ---
  # Combine colour and fill into one discrete scale
  scale_colour_futurama(aesthetics = c("colour", "fill"),
                        guide = guide_legend(title = "Probe")) +
  # Remove legends for alpha and size
  scale_alpha(range = c(0.4, 1), guide = "none") +
  scale_size_continuous(range = c(1, 3), guide = "none") +
  
  theme_void() +
  theme(aspect.ratio = 1,
        text = element_text(face = "bold"),
        legend.position = "right")

p

ggplot2::ggsave(filename="C:/Users/Lympha/Documents/Repositories/enERVate/results/plots/hipposideros_larvatus_main.png", plot = p, width = 20, height = 20, dpi = 300)
