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

source("../../scripts/hotspot/models.R")     # for score_windows_nb (used by recompute)
source("../../scripts/hotspot/postprocess.R")


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


# ------------------------- select_significant_windows ------------------

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


# ------------------------- merge_adjacent_hotspots ---------------------

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


# ------------------------- apply_min_hits_filter -------------------------

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


# ------------------------- assign_hotspot_ids + attach ------------------

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
    mu_nb   = c(1, 1, 1)
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


# ---------------- annotate_hotspot_composition (ADR-012) ----------------
# Detection runs on every locus in the tier; each called region is then
# DESCRIBED by what it contains. Annotation must never add, drop or re-score a
# region - only attach columns.

.loci_gr <- function() {
  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1", "chr2"),
    ranges = IRanges::IRanges(start = c(10L, 20L, 30L, 10L),
                              end   = c(15L, 25L, 35L, 15L))
  )
  S4Vectors::mcols(gr)$structure_class <- c("full", "gene", "gene", "full")
  S4Vectors::mcols(gr)$source <- c("ltr-flanked", "orphan", "ltr-flanked",
                                   "ltr-flanked")
  S4Vectors::mcols(gr)$segment <- c("Gammaretrovirus", "Betaretrovirus",
                                    "Gammaretrovirus", "Betaretrovirus")
  S4Vectors::mcols(gr)$confidence <- c("1.0", "0.5", "0.9", "0.8")
  gr
}

.regions_gr <- function() {
  GenomicRanges::GRanges(
    seqnames = c("chr1", "chr2"),
    ranges = IRanges::IRanges(start = c(1L, 1L), end = c(100L, 100L))
  )
}

test_that("composition counts structural classes and tiers per region", {
  out <- annotate_hotspot_composition(.regions_gr(), .loci_gr(), "segment")
  mc <- S4Vectors::mcols(out)
  expect_equal(mc$n_loci, c(3L, 1L))
  expect_equal(mc$n_full, c(1L, 1L))
  expect_equal(mc$n_gene, c(2L, 0L))
  expect_equal(mc$n_partial, c(0L, 0L))
  expect_equal(mc$n_ltr_flanked, c(2L, 1L))
  expect_equal(mc$n_orphan, c(1L, 0L))
})

test_that("composition reports the dominant lineage and mean confidence", {
  out <- annotate_hotspot_composition(.regions_gr(), .loci_gr(), "segment")
  mc <- S4Vectors::mcols(out)
  expect_equal(mc$dominant_taxon[1], "Gammaretrovirus")   # 2 of 3 on chr1
  expect_equal(mc$dominant_taxon[2], "Betaretrovirus")
  expect_equal(mc$mean_confidence[1], mean(c(1.0, 0.5, 0.9)))
})

test_that("composition never changes the number of regions", {
  regions <- .regions_gr()
  out <- annotate_hotspot_composition(regions, .loci_gr(), "segment")
  expect_equal(length(out), length(regions))
  expect_equal(BiocGenerics::start(out), BiocGenerics::start(regions))
})

test_that("composition is empty-safe for no regions and for no loci", {
  empty_loci <- GenomicRanges::GRanges()
  out <- annotate_hotspot_composition(.regions_gr(), empty_loci, "segment")
  expect_equal(S4Vectors::mcols(out)$n_loci, c(0L, 0L))
  expect_equal(length(annotate_hotspot_composition(GenomicRanges::GRanges(),
                                                   .loci_gr(), "segment")), 0L)
})

test_that("composition tolerates a raw-hit input carrying none of the columns", {
  bare <- GenomicRanges::GRanges(
    seqnames = "chr1", ranges = IRanges::IRanges(start = 10L, end = 15L)
  )
  out <- annotate_hotspot_composition(.regions_gr(), bare, "segment")
  mc <- S4Vectors::mcols(out)
  expect_equal(mc$n_loci, c(1L, 0L))        # still counts events
  expect_equal(mc$n_full, c(0L, 0L))        # but no structural breakdown
  expect_true(all(is.na(mc$dominant_taxon)))
})
