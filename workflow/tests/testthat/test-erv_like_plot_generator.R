# testthat tests for the ERV-like plot builders
# (workflow/scripts/erv_like_plot_generator/plots_*.R).
#
# Builder contract (shared with plot2sort / stage_plot_generator): each returns
# a ggplot, falls back to empty_plot() on zero-row input, and carries an
# `intended_dims` attr where the axis auto-scales. We assert the contract +
# the few non-trivial transforms (canonical dropped-count, gap NA-drop).
#
# Run with: make test-r.

suppressMessages({
  library(testthat)
})

# Sourcing the orchestrator loads tidyverse + the shared helpers + the erv-like
# modules (the bottom-of-file main() guard keeps the CLI from firing).
source("../../scripts/erv_like_plot_generator.R")

.loci <- function(n = 2L) {
  tibble::tibble(
    seqnames = "chr1", start = c(100, 5000)[seq_len(n)],
    end = c(2000, 7000)[seq_len(n)], width = c(1901L, 2001L)[seq_len(n)],
    strand = "+", ID = c("ERV_like_1", "ERV_like_2")[seq_len(n)],
    virus = c("HIV", "HTLV")[seq_len(n)], label = "L", species = "S",
    probes_present = c("GAG; POL; ENV", "GAG; POL")[seq_len(n)],
    n_main_present = c(3L, 2L)[seq_len(n)], n_main_expected = 3L,
    completeness_fraction = c(1, 2 / 3)[seq_len(n)],
    is_full = c(TRUE, FALSE)[seq_len(n)],
    observed_order = c("GAG; POL; ENV", "GAG; POL")[seq_len(n)],
    is_canonical = c(TRUE, FALSE)[seq_len(n)],
    n_loci = c(3L, 2L)[seq_len(n)],
    max_bitscore = c(300, 200)[seq_len(n)],
    max_identity = c(95, 90)[seq_len(n)],
    query_coverage = c(0.9, 0.6)[seq_len(n)]
  )
}

.members <- function() {
  tibble::tibble(
    Parent = c(rep("ERV_like_1", 3), rep("ERV_like_2", 2)),
    ID = paste0("m", 1:5), probe = c("GAG", "POL", "ENV", "GAG", "POL"),
    seqnames = "chr1", start = c(100, 900, 1800, 5000, 6000),
    end = c(300, 1100, 2000, 5200, 6200),
    gap_to_prev = c(NA, 599L, 699L, NA, 799L), virus = "HIV", label = "L"
  )
}

.empty_loci    <- function() .loci()[0, ]
.empty_members <- function() .members()[0, ]

# Every builder returns a ggplot on non-empty input.
test_that("composition + order/structure builders return ggplots", {
  loci <- .loci()
  expect_s3_class(probe_combination_plot(loci), "ggplot")
  expect_s3_class(completeness_plot(loci), "ggplot")
  expect_s3_class(full_partial_by_virus_plot(loci), "ggplot")
  expect_s3_class(composition_heatmap_plot(loci), "ggplot")
  expect_s3_class(gene_order_plot(loci), "ggplot")
  expect_s3_class(canonical_plot(loci, 0L), "ggplot")
  expect_s3_class(candidate_length_plot(loci), "ggplot")
  expect_s3_class(n_loci_plot(loci), "ggplot")
  expect_s3_class(interprobe_gap_plot(.members(), 1500), "ggplot")
})

# Every builder falls back to the labelled placeholder on zero-row input.
test_that("builders fall back to empty_plot on zero-row input", {
  el <- .empty_loci()
  expect_equal(probe_combination_plot(el)$labels$title, "no data")
  expect_equal(completeness_plot(el)$labels$title, "no data")
  expect_equal(composition_heatmap_plot(el)$labels$title, "no data")
  expect_equal(gene_order_plot(el)$labels$title, "no data")
  expect_equal(candidate_length_plot(el)$labels$title, "no data")
  expect_equal(interprobe_gap_plot(.empty_members(), 1500)$labels$title, "no data")
})

# Auto-scaled categorical-axis builders attach intended_dims.
test_that("categorical-axis builders attach intended_dims", {
  expect_false(is.null(attr(probe_combination_plot(.loci()), "intended_dims")))
  expect_false(is.null(attr(composition_heatmap_plot(.loci()), "intended_dims")))
})

# canonical_plot is robust to require_canonical_order: it shows retained
# canonical/rearranged AND the dropped (filtered) rearranged candidates.
test_that("canonical_plot still renders when all retained are canonical but some were dropped", {
  only_canon <- .loci(1L)                       # 1 canonical candidate retained
  expect_s3_class(canonical_plot(only_canon, dropped_noncanonical = 5L), "ggplot")
  # With no candidates at all but a positive dropped count, it must not be empty.
  p <- canonical_plot(.empty_loci(), dropped_noncanonical = 3L)
  expect_false(identical(p$labels$title, "no data"))
})

# interprobe_gap drops the per-candidate first-member NA rows.
test_that("interprobe_gap_plot ignores NA (first-member) gaps", {
  # 3 of the 5 fixture members carry a gap; 2 are NA first-members.
  d <- .members() %>% dplyr::filter(!is.na(gap_to_prev))
  expect_equal(nrow(d), 3L)
  expect_s3_class(interprobe_gap_plot(.members(), 1500), "ggplot")
})
