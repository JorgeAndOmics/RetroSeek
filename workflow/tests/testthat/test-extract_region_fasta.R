# testthat coverage for the Biostrings region extractor (bedtools getfasta
# replacement). Sources the script (CLI guarded by sys.nframe()) and exercises
# the pure extractor: 0-based BED -> 1-based inclusive, strand-aware revcomp.

suppressMessages({
  library(testthat)
  library(Biostrings)
})

.script_dir <- file.path("..", "..", "scripts")
source(file.path(.script_dir, "extract_region_fasta.R"))


test_that("extract_region_fasta cuts the right span and respects strand", {
  genome <- Biostrings::DNAStringSet(c(chr1 = "AAAACCCCGGGGTTTT"))
  bed <- data.frame(
    chrom  = c("chr1", "chr1"),
    start0 = c(0L, 4L),     # 0-based half-open
    end    = c(4L, 8L),
    name   = c("plus", "minus"),
    strand = c("+", "-"),
    stringsAsFactors = FALSE
  )
  out <- extract_region_fasta(genome, bed)
  expect_equal(names(out), c("plus", "minus"))
  expect_equal(as.character(out[["plus"]]), "AAAA")             # [1,4] forward
  expect_equal(as.character(out[["minus"]]), "GGGG")            # [5,8]=CCCC -> revcomp GGGG
})

test_that("extract_region_fasta is empty-safe", {
  genome <- Biostrings::DNAStringSet(c(chr1 = "ACGT"))
  bed <- data.frame(chrom = character(), start0 = integer(), end = integer(),
                    name = character(), strand = character())
  expect_length(extract_region_fasta(genome, bed), 0L)
})

test_that("extract_region_fasta errors on an absent seqname (mirrors bedtools)", {
  genome <- Biostrings::DNAStringSet(c(chr1 = "ACGT"))
  bed <- data.frame(chrom = "chrX", start0 = 0L, end = 2L, name = "x", strand = "+",
                    stringsAsFactors = FALSE)
  expect_error(extract_region_fasta(genome, bed))
})
