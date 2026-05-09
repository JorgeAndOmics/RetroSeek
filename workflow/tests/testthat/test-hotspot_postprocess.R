# testthat tests for hotspot_analysis/postprocess.R
#
# Run via: make test-r

suppressMessages({
  library(testthat)
  library(tibble)
  library(dplyr)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})

source("../../scripts/hotspot_analysis/models.R")     # for score_windows_nb (used by recompute)
source("../../scripts/hotspot_analysis/postprocess.R")


# Helper: significant-windows tibble matching the schema produced upstream.
.fake_significant_df <- function(rows) {
  tibble::tibble(
    chrom         = rows$chrom,
    chrom_stratum = rows$chrom,
    start         = rows$start,
    end           = rows$end,
    count         = rows$count,
    effective_bp  = rep(10000L, length(rows$start)),
    label         = "Gammaretrovirus",
    qval_nb       = rep(0.001, length(rows$start))
  )
}


# ───────────────────────── select_significant_windows ──────────────────

test_that("select_significant_windows filters by q-value threshold", {
  df <- tibble::tibble(
    chrom = "chr1", chrom_stratum = "chr1",
    start = c(1L, 1001L, 2001L), end = c(1000L, 2000L, 3000L),
    count = c(10L, 1L, 5L), effective_bp = 1000L, label = "Ungrouped",
    qval_nb = c(0.001, 0.5, 0.04)
  )
  out <- select_significant_windows(df, threshold = 0.05)
  expect_equal(nrow(out), 2L)
  expect_true(all(out$qval_nb < 0.05))
})

test_that("select_significant_windows skips NA q-values silently", {
  df <- tibble::tibble(
    chrom = "chr1", chrom_stratum = "chr1",
    start = 1L, end = 1000L, count = 0L, effective_bp = 0L,
    label = "Ungrouped", qval_nb = NA_real_
  )
  out <- select_significant_windows(df, threshold = 0.05)
  expect_equal(nrow(out), 0L)
})


# ───────────────────────── merge_adjacent_hotspots ─────────────────────

test_that("merge_adjacent_hotspots merges strictly adjacent windows when gap=0", {
  rows <- list(
    chrom = c("chr1", "chr1"),
    start = c(1L, 1001L),
    end   = c(1000L, 2000L),
    count = c(5L, 7L)
  )
  merged <- merge_adjacent_hotspots(.fake_significant_df(rows), gap = 0L)
  expect_length(merged, 1L)
  expect_equal(BiocGenerics::start(merged), 1L)
  expect_equal(BiocGenerics::end(merged), 2000L)
  expect_equal(S4Vectors::mcols(merged)$count, 12L)
  expect_equal(S4Vectors::mcols(merged)$n_windows, 2L)
})

test_that("merge_adjacent_hotspots keeps non-adjacent windows separate when gap=0", {
  rows <- list(
    chrom = c("chr1", "chr1"),
    start = c(1L, 5001L),     # 4000bp gap between them
    end   = c(1000L, 6000L),
    count = c(5L, 7L)
  )
  merged <- merge_adjacent_hotspots(.fake_significant_df(rows), gap = 0L)
  expect_length(merged, 2L)
})

test_that("merge_adjacent_hotspots respects gap=large to bridge non-adjacent windows", {
  rows <- list(
    chrom = c("chr1", "chr1"),
    start = c(1L, 5001L),
    end   = c(1000L, 6000L),
    count = c(5L, 7L)
  )
  merged <- merge_adjacent_hotspots(.fake_significant_df(rows), gap = 5000L)
  expect_length(merged, 1L)
})

test_that("merge_adjacent_hotspots returns canonical empty schema for empty input", {
  out <- merge_adjacent_hotspots(.fake_significant_df(list(
    chrom = character(0), start = integer(0),
    end = integer(0), count = integer(0)
  )))
  expect_length(out, 0L)
  expect_true(all(c("label", "count", "effective_bp", "n_windows",
                    "chrom_stratum", "mu_nb_region", "pval_nb_region")
                  %in% colnames(S4Vectors::mcols(out))))
})


# ───────────────────────── apply_min_hits_filter ─────────────────────────

test_that("apply_min_hits_filter drops regions below the threshold", {
  rows <- list(
    chrom = c("chr1", "chr1"),
    start = c(1L, 5001L),
    end   = c(1000L, 6000L),
    count = c(5L, 1L)
  )
  merged <- merge_adjacent_hotspots(.fake_significant_df(rows))
  out <- apply_min_hits_filter(merged, min_hits = 3L)
  expect_length(out, 1L)
  expect_equal(S4Vectors::mcols(out)$count, 5L)
})

test_that("apply_min_hits_filter is a pass-through when min_hits <= 0", {
  rows <- list(
    chrom = c("chr1"), start = c(1L), end = c(1000L), count = c(2L)
  )
  merged <- merge_adjacent_hotspots(.fake_significant_df(rows))
  expect_equal(length(apply_min_hits_filter(merged, min_hits = 0L)),
               length(merged))
})


# ───────────────────────── assign_hotspot_ids + attach ──────────────────

test_that("assign_hotspot_ids stamps zero-padded IDs prefixed by species", {
  rows <- list(
    chrom = c("chr1", "chr1"), start = c(1L, 5001L),
    end = c(1000L, 6000L), count = c(5L, 7L)
  )
  merged <- merge_adjacent_hotspots(.fake_significant_df(rows))
  out <- assign_hotspot_ids(merged, species = "Antrozous_pallidus")
  expect_equal(S4Vectors::mcols(out)$hotspot_id,
               c("Antrozous_pallidus_HS_00001",
                 "Antrozous_pallidus_HS_00002"))
})

test_that("attach_hotspot_id_to_windows propagates IDs to overlapping windows and NA elsewhere", {
  win_df <- tibble::tibble(
    chrom = c("chr1", "chr1", "chr1"),
    chrom_stratum = "chr1",
    start = c(1L, 1001L, 8001L),
    end   = c(1000L, 2000L, 9000L),
    count = c(5L, 3L, 0L),
    effective_bp = 1000L,
    label = "Ungrouped",
    qval_nb = c(0.001, 0.005, 0.5),
    pval_nb = c(0.0005, 0.001, 0.4),
    mu_nb   = c(1, 1, 1),
    pval_perm = NA_real_, qval_perm = NA_real_
  )
  rows <- list(
    chrom = c("chr1", "chr1"), start = c(1L, 1001L),
    end = c(1000L, 2000L), count = c(5L, 3L)
  )
  merged <- merge_adjacent_hotspots(.fake_significant_df(rows))
  merged <- assign_hotspot_ids(merged, species = "X")
  attached <- attach_hotspot_id_to_windows(win_df, merged)
  expect_equal(attached$hotspot_id, c("X_HS_00001", "X_HS_00001", NA_character_))
})
