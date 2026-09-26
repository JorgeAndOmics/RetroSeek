# =============================================================================
# extract_region_fasta.R
# =============================================================================
# Bioconductor (Biostrings) replacement for `bedtools getfasta -s -nameOnly`.
# RetroSeek's range/sequence work stays inside Bioconductor rather than shelling
# out to bedtools: given a genome FASTA and a BED6 of regions, extract each
# region, reverse-complement minus-strand regions, and name each record by the
# BED name column (col 4). Used by taxonomy_classify_loci.py to cut per-locus
# marker regions before blastx.
#
# Memory note: Biostrings loads the genome into a DNAStringSet (~1 byte/base).
# That is fine for the model genomes (human ~ 3.2 GB, ample headroom here); for
# very large genomes run with bounded parallelism. `testthat` sources this file
# (the CLI block is guarded by sys.nframe()), so the pure extractor is testable
# without invoking the command line.

suppressMessages({
  library(Biostrings)
})


# Pure extractor: `genome` is a named DNAStringSet (names = seqnames), `bed` a
# data.frame with columns chrom, start0 (0-based, half-open), end, name, strand.
# Returns a DNAStringSet named by `bed$name`, minus-strand regions revcomp'd.
extract_region_fasta <- function(genome, bed) {
  if (nrow(bed) == 0L) return(Biostrings::DNAStringSet())
  pieces <- vapply(seq_len(nrow(bed)), function(i) {
    chr <- genome[[bed$chrom[i]]]  # errors on an absent seqname, mirroring bedtools
    # BED starts are 0-based; subseq() takes 1-based inclusive positions.
    sub <- Biostrings::subseq(chr,
                              start = bed$start0[i] + 1L,
                              end = bed$end[i])
    if (identical(bed$strand[i], "-")) sub <- Biostrings::reverseComplement(sub)
    as.character(sub)
  }, character(1))
  out <- Biostrings::DNAStringSet(pieces)
  names(out) <- bed$name
  out
}


main <- function() {
  suppressMessages(library(argparse))
  # This file's folder, to find utils/log.R beside it (Rscript passes --file=).
  script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  source(file.path(dirname(script_file), "..", "utils", "log.R"))
  parser <- ArgumentParser(
    description =
      "Strand-aware region FASTA extraction (Biostrings; replaces bedtools getfasta)"
  )
  parser$add_argument("--genome", required = TRUE, help = "genome FASTA")
  parser$add_argument("--bed", required = TRUE,
                      help = "regions BED6 (chrom, start0, end, name, score, strand)")
  parser$add_argument("--out", required = TRUE, help = "output FASTA")
  parser$add_argument(
    "--log", default = NULL,
    help = "job log file; the Snakemake log: path of the calling rule"
  )
  args <- parser$parse_args()
  log_job(args$log, "extract_region_fasta")
  run_main(function() extract(args))
}


# Read the genome and the BED, write the regions. A helper step inside other jobs,
# so it reports at INFO: the calling job prints the one OK headline.
extract <- function(args) {
  genome <- Biostrings::readDNAStringSet(args$genome)
  # Key on the first whitespace token of each header, as bedtools does, so a
  # header like "chr1 description" is addressed as "chr1".
  names(genome) <- sub("\\s.*$", "", names(genome))

  bed <- utils::read.table(
    args$bed, sep = "\t", stringsAsFactors = FALSE,
    col.names = c("chrom", "start0", "end", "name", "score", "strand"),
    colClasses = c("character", "integer", "integer", "character", "character",
                   "character")

  )
  regions <- extract_region_fasta(genome, bed)
  Biostrings::writeXStringSet(regions, args$out)
  log_info("%s regions extracted", format(length(regions), big.mark = ","))
}


if (sys.nframe() == 0L) main()
