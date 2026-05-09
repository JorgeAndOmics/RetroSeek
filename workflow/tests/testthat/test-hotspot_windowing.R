# testthat tests for hotspot_analysis/windowing.R
#
# Run via: make test-r

suppressMessages({
  library(testthat)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})

source("../../scripts/hotspot_analysis/windowing.R")


# ─────────────────────────── tile_genome_for_hotspot ───────────────────────

test_that("tile_genome_for_hotspot produces full-width tiles for a genome divisible by window_size", {
  out <- tile_genome_for_hotspot(c(chr1 = 30000L), window_size = 10000L)
  expect_equal(length(out), 3L)
  expect_equal(BiocGenerics::width(out), c(10000L, 10000L, 10000L))
})

test_that("tile_genome_for_hotspot's last tile is shorter when chrom is not divisible", {
  out <- tile_genome_for_hotspot(c(chr1 = 25000L), window_size = 10000L)
  expect_equal(length(out), 3L)
  # last tile shorter: 25000 - 20000 = 5000
  expect_equal(BiocGenerics::width(out), c(10000L, 10000L, 5000L))
})

test_that("tile_genome_for_hotspot tiles each chromosome independently", {
  out <- tile_genome_for_hotspot(c(chr1 = 15000L, chr2 = 20000L), window_size = 10000L)
  # chr1: 10000 + 5000; chr2: 10000 + 10000 → 4 tiles total
  expect_equal(length(out), 4L)
})


# ─────────────────────────── count_hits_per_window ───────────────────────────

test_that("count_hits_per_window counts overlapping hits exactly", {
  windows <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges   = IRanges::IRanges(start = c(1, 1001), end = c(1000, 2000))
  )
  hits <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(start = c(100, 500, 1500),
                                end   = c(110, 510, 1510))
  )
  expect_equal(count_hits_per_window(windows, hits), c(2L, 1L))
})

test_that("count_hits_per_window returns zeros for empty hits", {
  windows <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(start = c(1, 1001), end = c(1000, 2000))
  )
  out <- count_hits_per_window(windows, GenomicRanges::GRanges())
  expect_equal(out, c(0L, 0L))
})


# ─────────────────────────── assemble_window_table ──────────────────────────

test_that("assemble_window_table builds a tibble with the expected columns", {
  windows <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges   = IRanges::IRanges(start = c(1, 1001), end = c(1000, 2000))
  )
  out <- assemble_window_table(
    windows,
    counts        = c(3L, 0L),
    effective_bp  = c(1000L, 950L),
    chrom_stratum = c(chr1 = "chr1"),
    label         = "Gammaretrovirus"
  )
  expect_s3_class(out, "tbl_df")
  expect_named(out, c("chrom", "chrom_stratum", "start", "end",
                      "count", "effective_bp", "label"))
  expect_equal(out$chrom, c("chr1", "chr1"))
  expect_equal(out$chrom_stratum, c("chr1", "chr1"))
  expect_equal(out$count, c(3L, 0L))
  expect_equal(out$effective_bp, c(1000L, 950L))
  expect_equal(unique(out$label), "Gammaretrovirus")
})

test_that("assemble_window_table falls back to the raw chrom name when stratum map is missing one", {
  windows <- GenomicRanges::GRanges(
    seqnames = "scaffold_99",
    ranges   = IRanges::IRanges(start = 1, end = 1000)
  )
  out <- assemble_window_table(
    windows,
    counts        = 0L,
    effective_bp  = 1000L,
    # stratum map keyed only by chr1 — scaffold_99 unmapped
    chrom_stratum = c(chr1 = "chr1"),
    label         = "Ungrouped"
  )
  expect_equal(out$chrom_stratum, "scaffold_99")
})
