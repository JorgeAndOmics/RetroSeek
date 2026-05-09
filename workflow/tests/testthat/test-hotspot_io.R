# testthat tests for hotspot_analysis/io.R + utils/chrom_names.R
#
# Run via: make test-r

suppressMessages({
  library(testthat)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})

source("../../scripts/utils/chrom_names.R")
source("../../scripts/hotspot_analysis/io.R")


# ──────────────────────── normalise_chrom_names ─────────────────────────

test_that("normalise_chrom_names extracts NCBI accession tokens from full headers", {
  headers <- c(
    "CM034567.1 Genus species chromosome 1, GRCh38",
    "NC_000001.11 Homo sapiens chromosome 1",
    "JH791234.1 unplaced scaffold"
  )
  expect_equal(
    suppressMessages(normalise_chrom_names(headers)),
    c("CM034567.1", "NC_000001.11", "JH791234.1")
  )
})

test_that("normalise_chrom_names returns NA for non-matching headers and emits a message", {
  headers <- c("CM034567.1 ok", "chr1 custom-format-no-match")
  expect_message(
    out <- normalise_chrom_names(headers),
    regexp = "did not match pattern"
  )
  expect_equal(out, c("CM034567.1", NA_character_))
})

test_that("normalise_chrom_names is silent when all headers match", {
  headers <- c("CM034567.1 a", "CM034568.1 b")
  expect_silent(normalise_chrom_names(headers))
})


# ─────────────────────────── read_hotspot_options ───────────────────────────

test_that("read_hotspot_options returns documented defaults when params block is empty", {
  config <- list(parameters = list())
  opts <- read_hotspot_options(config)
  expect_equal(opts$seed, 67L)
  expect_equal(opts$input, "valid")
  expect_false(opts$group_split)
  expect_equal(opts$window_size, 10000L)
  expect_equal(opts$mask_size, 20L)
  expect_equal(opts$pvalue_threshold, 0.05)
  expect_false(opts$validate_permutation)
  expect_true(opts$merge_adjacent)
  expect_equal(opts$merge_gap, 0L)
  expect_true(opts$strata_by_chromosome)
  expect_equal(opts$unplaced_min_factor, 10L)
})

test_that("read_hotspot_options propagates user overrides", {
  config <- list(parameters = list(
    seed                          = 1234,
    hotspot_input                 = "original",
    hotspot_window_size           = 5000,
    hotspot_validate_permutation  = TRUE,
    hotspot_merge_gap             = 100,
    hotspot_strata_by_chromosome  = FALSE
  ))
  opts <- read_hotspot_options(config)
  expect_equal(opts$seed, 1234L)
  expect_equal(opts$input, "original")
  expect_equal(opts$window_size, 5000L)
  expect_true(opts$validate_permutation)
  expect_equal(opts$merge_gap, 100L)
  expect_false(opts$strata_by_chromosome)
})


# ─────────────────────────── load_hits_gff (label assertion) ───────────────

test_that("load_hits_gff aborts when the GFF lacks an mcols$label column", {
  # Build a minimal GFF3 file via rtracklayer::export with no label column
  gr <- GenomicRanges::GRanges(
    seqnames = "CM000001.1",
    ranges   = IRanges::IRanges(start = 100, end = 200)
  )
  S4Vectors::mcols(gr)$type <- "feature"  # no `label`
  tmp <- tempfile(fileext = ".gff3")
  on.exit(unlink(tmp), add = TRUE)
  rtracklayer::export(gr, tmp, format = "gff3")
  expect_error(load_hits_gff(tmp), regexp = "label")
})

test_that("load_hits_gff returns the imported GRanges when label is present", {
  gr <- GenomicRanges::GRanges(
    seqnames = "CM000001.1",
    ranges   = IRanges::IRanges(start = 100, end = 200)
  )
  S4Vectors::mcols(gr)$label <- "Gammaretrovirus"
  tmp <- tempfile(fileext = ".gff3")
  on.exit(unlink(tmp), add = TRUE)
  rtracklayer::export(gr, tmp, format = "gff3")
  out <- load_hits_gff(tmp)
  expect_s4_class(out, "GRanges")
  expect_true("label" %in% colnames(S4Vectors::mcols(out)))
})
