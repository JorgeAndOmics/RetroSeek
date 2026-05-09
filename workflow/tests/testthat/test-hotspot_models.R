# testthat tests for hotspot_analysis/models.R
#
# Run via: make test-r
#
# These tests focus on the deterministic, in-process branches:
#   * NB GLM fits when we plant a clear hotspot in synthetic counts.
#   * `score_windows_nb` produces well-formed p / q columns and respects NA
#     rules for windows with effective_bp == 0.
#   * The convergence ladder degrades cleanly when the data is too sparse.
#
# Permutation tests (`score_windows_perm`) are integration-level; we verify
# only that the call runs deterministically given a fixed seed and stashes
# the per-perm totals attribute used by the histogram/density plots.

suppressMessages({
  library(testthat)
  library(tibble)
  library(dplyr)
  library(GenomicRanges)
  library(IRanges)
})

source("../../scripts/hotspot_analysis/models.R")


# Helper: synthetic per-window tibble with a Poisson-ish baseline plus a
# planted hotspot at the requested rows. Makes the tests deterministic.
.fake_window_df <- function(n = 200L,
                            baseline_rate = 0.5,
                            hotspot_rows = NULL,
                            hotspot_count = 50L,
                            window_size = 10000L,
                            seed = 1L) {
  set.seed(seed)
  counts <- stats::rpois(n, baseline_rate)
  if (!is.null(hotspot_rows)) {
    counts[hotspot_rows] <- hotspot_count
  }
  tibble::tibble(
    chrom         = rep(c("chr1", "chr2"), length.out = n),
    chrom_stratum = rep(c("chr1", "chr2"), length.out = n),
    start         = seq.int(1L, by = window_size, length.out = n),
    end           = seq.int(1L, by = window_size, length.out = n) + window_size - 1L,
    count         = as.integer(counts),
    effective_bp  = rep(window_size, n),
    label         = "Ungrouped"
  )
}


# ───────────────────────── fit_nb_model ─────────────────────────────────

test_that("fit_nb_model returns insufficient_data when all windows are zero", {
  df <- .fake_window_df(n = 100L, baseline_rate = 0)  # all zeros
  fit <- fit_nb_model(df, window_size = 10000L,
                      strata_by_chromosome = TRUE)
  expect_equal(fit$status, "insufficient_data")
  expect_null(fit$model)
  expect_true(is.na(fit$theta))
})

test_that("fit_nb_model fits successfully on a dataset with planted enrichment", {
  df <- .fake_window_df(n = 400L, baseline_rate = 1,
                        hotspot_rows = c(50, 51, 200, 201),
                        hotspot_count = 50L)
  fit <- fit_nb_model(df, window_size = 10000L,
                      strata_by_chromosome = TRUE)
  expect_true(fit$status %in% c("ok", "poisson_fallback"))
  expect_false(is.null(fit$model))
})

test_that("fit_nb_model honours strata_by_chromosome=FALSE", {
  df <- .fake_window_df(n = 400L, baseline_rate = 1,
                        hotspot_rows = 100,
                        hotspot_count = 30L)
  fit <- fit_nb_model(df, window_size = 10000L,
                      strata_by_chromosome = FALSE)
  expect_true(fit$status %in% c("ok", "poisson_fallback"))
  # No chromosome term means a single intercept
  expect_false(fit$strata)
})


# ───────────────────────── score_windows_nb ────────────────────────────

test_that("score_windows_nb adds mu_nb / pval_nb / qval_nb columns", {
  df <- .fake_window_df(n = 400L, baseline_rate = 1,
                        hotspot_rows = c(50, 51, 200, 201),
                        hotspot_count = 50L)
  fit <- fit_nb_model(df, window_size = 10000L,
                      strata_by_chromosome = TRUE)
  scored <- score_windows_nb(df, fit)
  expect_true(all(c("mu_nb", "pval_nb", "qval_nb") %in% colnames(scored)))
  expect_equal(nrow(scored), nrow(df))
})

test_that("score_windows_nb's planted hotspot rows have markedly smaller p-values than the median", {
  df <- .fake_window_df(n = 400L, baseline_rate = 1,
                        hotspot_rows = c(50, 51, 200, 201),
                        hotspot_count = 50L)
  fit <- fit_nb_model(df, window_size = 10000L,
                      strata_by_chromosome = TRUE)
  scored <- score_windows_nb(df, fit)
  # Skip if model didn't converge meaningfully (rare); the planted hotspot
  # signal should otherwise dominate.
  skip_if(is.null(fit$model), "model did not fit; cannot score")
  hotspot_p <- scored$pval_nb[c(50, 51, 200, 201)]
  median_p  <- stats::median(scored$pval_nb, na.rm = TRUE)
  expect_true(all(hotspot_p < median_p, na.rm = TRUE))
  expect_true(all(scored$qval_nb[c(50, 51, 200, 201)] < 0.05, na.rm = TRUE))
})

test_that("score_windows_nb returns NA p-values for windows with effective_bp == 0", {
  df <- .fake_window_df(n = 200L, baseline_rate = 1,
                        hotspot_rows = 100, hotspot_count = 30L)
  df$effective_bp[10] <- 0L
  fit <- fit_nb_model(df, window_size = 10000L,
                      strata_by_chromosome = TRUE)
  scored <- score_windows_nb(df, fit)
  expect_true(is.na(scored$pval_nb[10]))
  expect_true(is.na(scored$qval_nb[10]))
})

test_that("score_windows_nb produces NA p-values when the fit failed entirely", {
  df <- .fake_window_df(n = 100L, baseline_rate = 0)
  fit <- fit_nb_model(df, window_size = 10000L,
                      strata_by_chromosome = TRUE)
  expect_equal(fit$status, "insufficient_data")
  scored <- score_windows_nb(df, fit)
  expect_true(all(is.na(scored$pval_nb)))
})


# ───────────────────────── score_windows_perm ──────────────────────────

test_that("score_windows_perm runs deterministically given the same seed", {
  skip_if_not_installed("regioneR")
  windows <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 50L),
    ranges   = IRanges::IRanges(
      start = seq.int(1L, by = 1000L, length.out = 50L),
      end   = seq.int(1L, by = 1000L, length.out = 50L) + 999L
    )
  )
  hits <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(start = c(1500, 12500, 23500),
                                end   = c(1510, 12510, 23510))
  )
  genome_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(start = 1L, end = 50000L)
  )
  win_df <- tibble::tibble(
    chrom = "chr1", chrom_stratum = "chr1",
    start = BiocGenerics::start(windows),
    end   = BiocGenerics::end(windows),
    count = as.integer(GenomicRanges::countOverlaps(windows, hits)),
    effective_bp = 1000L,
    label = "Ungrouped"
  )
  out_a <- score_windows_perm(win_df, windows, hits, genome_gr,
                              GenomicRanges::GRanges(),
                              n_perm = 5L, seed = 67L)
  out_b <- score_windows_perm(win_df, windows, hits, genome_gr,
                              GenomicRanges::GRanges(),
                              n_perm = 5L, seed = 67L)
  expect_equal(out_a$pval_perm, out_b$pval_perm)
  expect_equal(out_a$qval_perm, out_b$qval_perm)
  expect_equal(length(attr(out_a, "perm_totals")), 5L)
})
