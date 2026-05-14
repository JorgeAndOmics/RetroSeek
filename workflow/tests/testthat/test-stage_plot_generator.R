# testthat smoke tests for workflow/scripts/stage_plot_generator/*.R
#
# The 8 middle-stage plot builders. Each must return a ggplot for valid input
# and fall through to empty_plot() (still a ggplot) for empty input. Shared
# helpers (theme, add_titles, auto_dims, ...) are reused from plot2sort.

suppressMessages({
  library(testthat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(scales)
  library(ggsci)
})

.scripts <- file.path("..", "..", "scripts")
source(file.path(.scripts, "plot2sort", "helpers.R"))
source(file.path(.scripts, "stage_plot_generator", "plots_concordance.R"))
source(file.path(.scripts, "stage_plot_generator", "plots_structure.R"))
source(file.path(.scripts, "stage_plot_generator", "plots_funnel.R"))
source(file.path(.scripts, "stage_plot_generator", "plots_multiplicity.R"))


.fake_hits_df <- function() {
  tibble::tibble(
    genome         = "G1",
    probe          = c("POL", "POL", "GAG", "ENV", "ENV"),
    virus          = "HIV",
    label          = "Lentivirus",
    species        = "G1",
    n_hits         = c(5L, 1L, 3L, 8L, 2L),
    max_bitscore   = c(300, 90, 150, 250, 110),
    max_identity   = c(95, 70, 82, 90, 75),
    query_coverage = c(0.9, 0.4, 0.6, 0.85, 0.5),
    concordance    = c("inside", "flanking", "inside", "disjoint", "inside"),
    is_candidate   = c(TRUE, TRUE, TRUE, FALSE, TRUE),
    is_valid       = c(TRUE, FALSE, TRUE, FALSE, FALSE)
  )
}

.fake_ltr_df <- function() {
  tibble::tibble(
    genome             = "G1",
    ID                 = paste0("retro_", 1:4),
    n_flanking_ltrs    = c(2L, 1L, 2L, 0L),
    has_both_ltrs      = c(TRUE, FALSE, TRUE, FALSE),
    n_domains          = c(3L, 0L, 1L, 0L),
    domain_probes      = c("GAG; POL", NA, "ENV", NA),
    has_tsd            = c(TRUE, FALSE, FALSE, FALSE),
    n_overlapping_hits = c(4L, 0L, 2L, 1L)
  )
}

.fake_reduced_df <- function() {
  tibble::tibble(
    genome = "G1", ID = paste0("POL_", 1:3),
    probe = "POL", virus = "HIV", label = "Lentivirus",
    n_hits = c(10L, 4L, 2L), n_loci = c(4L, 2L, 1L),
    max_bitscore = c(300, 200, 110)
  )
}

.fake_counts_df <- function() {
  lvls <- c("homology hits", "first-reduced", "candidate", "valid")
  tibble::tibble(
    genome = rep(c("G1", "G2"), each = 4),
    stage  = factor(rep(lvls, 2), levels = lvls),
    count  = c(10000, 4000, 1200, 600, 8000, 3000, 900, 300)
  )
}

is_gg <- function(x) inherits(x, "ggplot")


test_that("concordance + probe-yield builders return ggplots for valid input", {
  expect_true(is_gg(concordance_plot(.fake_hits_df())))
  expect_true(is_gg(probe_yield_plot(.fake_hits_df())))
})

test_that("LTR structure builders return ggplots for valid input", {
  expect_true(is_gg(ltr_completeness_plot(.fake_ltr_df())))
  expect_true(is_gg(domain_composition_plot(.fake_ltr_df())))
})

test_that("funnel builders return ggplots for valid input", {
  expect_true(is_gg(refinement_funnel_plot(.fake_counts_df())))
  expect_true(is_gg(aggregate_funnel_plot(.fake_counts_df())))
})

test_that("multiplicity builders return ggplots for valid input", {
  expect_true(is_gg(multiplicity_m1_plot(.fake_hits_df())))
  expect_true(is_gg(multiplicity_m2_plot(.fake_reduced_df())))
})

test_that("every stage builder falls through to a ggplot for empty input", {
  e_hits    <- .fake_hits_df()[0, ]
  e_ltr     <- .fake_ltr_df()[0, ]
  e_reduced <- .fake_reduced_df()[0, ]
  e_counts  <- .fake_counts_df()[0, ]
  expect_true(is_gg(concordance_plot(e_hits)))
  expect_true(is_gg(probe_yield_plot(e_hits)))
  expect_true(is_gg(ltr_completeness_plot(e_ltr)))
  expect_true(is_gg(domain_composition_plot(e_ltr)))
  expect_true(is_gg(refinement_funnel_plot(e_counts)))
  expect_true(is_gg(aggregate_funnel_plot(e_counts)))
  expect_true(is_gg(multiplicity_m1_plot(e_hits)))
  expect_true(is_gg(multiplicity_m2_plot(e_reduced)))
})

test_that("auto-scaled stage builders set an intended_dims attribute", {
  expect_false(is.null(attr(concordance_plot(.fake_hits_df()), "intended_dims")))
  expect_false(is.null(attr(probe_yield_plot(.fake_hits_df()), "intended_dims")))
  expect_false(is.null(attr(domain_composition_plot(.fake_ltr_df()), "intended_dims")))
})

test_that("warning_caption is accepted and still yields a ggplot", {
  cap <- "⚠ test caption"
  expect_true(is_gg(concordance_plot(.fake_hits_df(), warning_caption = cap)))
  expect_true(is_gg(multiplicity_m1_plot(.fake_hits_df(), warning_caption = cap)))
  expect_true(is_gg(refinement_funnel_plot(.fake_counts_df(), warning_caption = cap)))
})
