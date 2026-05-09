# testthat tests for hotspot_analysis/masking.R
#
# Run via: make test-r

suppressMessages({
  library(testthat)
  library(Biostrings)
  library(GenomicRanges)
  library(IRanges)
})

source("../../scripts/hotspot_analysis/masking.R")


# ─────────────────────────── build_n_mask ───────────────────────────

test_that("build_n_mask returns an empty GRanges when mask_size == 0 (no masking)", {
  seqs <- Biostrings::DNAStringSet(c(chr1 = "ACGTACGTACGT"))
  out <- build_n_mask(seqs, mask_size = 0L, mask_mismatch = 0L)
  expect_s4_class(out, "GRanges")
  expect_length(out, 0L)
})

test_that("build_n_mask flags an N-run that exactly matches the motif length", {
  seqs <- Biostrings::DNAStringSet(c(chr1 = "ACGTNNNNNACGT"))
  out <- build_n_mask(seqs, mask_size = 5L, mask_mismatch = 0L)
  expect_length(out, 1L)
  expect_equal(as.character(GenomicRanges::seqnames(out)), "chr1")
})

test_that("build_n_mask reduces overlapping matches into a single interval", {
  # 8 Ns in a row with mask_size=5 produce multiple sliding-window matches;
  # they must collapse to one merged range so per-window masked-bp doesn't
  # double-count.
  seqs <- Biostrings::DNAStringSet(c(chr1 = "ACGTNNNNNNNNACGT"))
  out <- build_n_mask(seqs, mask_size = 5L, mask_mismatch = 0L)
  expect_length(out, 1L)
})


# ─────────────────────────── effective_bp_per_window ───────────────────────

test_that("effective_bp_per_window equals window width when mask is empty", {
  windows <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges   = IRanges::IRanges(start = c(1, 1001), end = c(1000, 2000))
  )
  out <- effective_bp_per_window(windows, GenomicRanges::GRanges())
  expect_equal(out, c(1000L, 1000L))
})

test_that("effective_bp_per_window subtracts only the in-window masked bp", {
  windows <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges   = IRanges::IRanges(start = c(1, 1001), end = c(1000, 2000))
  )
  # 100bp mask spanning the boundary 950-1050: 51bp in window 1, 50bp in window 2
  mask <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(start = 950, end = 1050)
  )
  out <- effective_bp_per_window(windows, mask)
  expect_equal(out, c(1000L - 51L, 1000L - 50L))
})

test_that("effective_bp_per_window never returns negative values", {
  windows <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(start = 1, end = 100)
  )
  # Mask exceeding the window width
  mask <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(start = 1, end = 100)
  )
  out <- effective_bp_per_window(windows, mask)
  expect_equal(out, 0L)
})

test_that("effective_bp_per_window handles the short last-tile correctly", {
  # Mimic cut.last.tile.in.chrom=TRUE: chr1 has total length 2500 with
  # window_size 1000 → tiles 1-1000, 1001-2000, 2001-2500 (last is 500bp).
  windows <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1"),
    ranges   = IRanges::IRanges(start = c(1, 1001, 2001),
                                end   = c(1000, 2000, 2500))
  )
  out <- effective_bp_per_window(windows, GenomicRanges::GRanges())
  expect_equal(out, c(1000L, 1000L, 500L))
})


# ─────────────────────────── pool_small_scaffolds ───────────────────────

test_that("pool_small_scaffolds keeps long chromosomes and pools short ones", {
  seqlengths <- c(chr1 = 1e8, chr2 = 5e7, scaffold_a = 5000, scaffold_b = 2000)
  out <- pool_small_scaffolds(seqlengths,
                               window_size = 10000L,
                               min_factor = 10L)
  # Threshold = 10000 * 10 = 100000. chr1, chr2 keep their names.
  # scaffolds < 100000 → "Unplaced".
  expect_equal(unname(out), c("chr1", "chr2", "Unplaced", "Unplaced"))
  expect_equal(names(out), names(seqlengths))
})

test_that("pool_small_scaffolds with min_factor=1 keeps all chroms >= window_size", {
  seqlengths <- c(chr1 = 10000, chr2 = 9999)
  out <- pool_small_scaffolds(seqlengths,
                               window_size = 10000L,
                               min_factor = 1L)
  expect_equal(unname(out), c("chr1", "Unplaced"))
})
