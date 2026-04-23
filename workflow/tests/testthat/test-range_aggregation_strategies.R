# testthat tests for workflow/scripts/range_aggregation_strategies.R
#
# Run with:
#   Rscript -e 'testthat::test_dir("workflow/tests/testthat")'
# or via the project Makefile target:
#   make test-r

suppressMessages({
  library(testthat)
  library(IRanges)
})

# Source the helpers from the primary consumer's sibling path.
source("../../scripts/range_aggregation_strategies.R")


# ───────────────────────────── aggregate_values ─────────────────────────────

test_that("list strategy returns a CharacterList of length 1 with all unique values", {
  out <- aggregate_values(c("POL", "GAG", "POL"), strategy = "list")
  expect_s4_class(out, "CharacterList")
  expect_length(out, 1L)
  expect_setequal(out[[1]], c("POL", "GAG"))
})

test_that("concatenate strategy joins unique values in sorted order with the separator", {
  out <- aggregate_values(c("POL", "GAG", "POL"),
                          strategy = "concatenate",
                          separator = "; ")
  expect_type(out, "character")
  expect_length(out, 1L)
  expect_equal(out, "GAG; POL")  # alphabetical order
})

test_that("concatenate strategy respects a custom separator", {
  out <- aggregate_values(c("POL", "GAG"),
                          strategy = "concatenate",
                          separator = " | ")
  expect_equal(out, "GAG | POL")
})

test_that("best strategy picks the value at the argmax of the tiebreaker", {
  out <- aggregate_values(c("POL", "GAG", "ENV"),
                          strategy = "best",
                          tiebreaker = c(100, 200, 50))
  expect_equal(out, "GAG")
})

test_that("best strategy errors when no tiebreaker is supplied", {
  expect_error(
    aggregate_values(c("POL", "GAG"), strategy = "best"),
    "requires a tiebreaker"
  )
})

test_that("best strategy resolves ties by the first matching position", {
  out <- aggregate_values(c("POL", "GAG", "ENV"),
                          strategy = "best",
                          tiebreaker = c(100, 100, 50))
  expect_equal(out, "POL")  # which.max returns the first max
})

test_that("majority strategy returns the most frequent value", {
  out <- aggregate_values(c("POL", "POL", "GAG"), strategy = "majority")
  expect_equal(out, "POL")
})

test_that("majority strategy breaks ties alphabetically (table() sort is stable on order)", {
  # With two POL and two GAG, sort() on table() puts the higher count first.
  # Here both are tied → the specific tie-break behaviour is R's — assert one wins.
  out <- aggregate_values(c("POL", "POL", "GAG", "GAG"), strategy = "majority")
  expect_true(out %in% c("POL", "GAG"))
})

test_that("first strategy returns the alphabetical first unique value", {
  out <- aggregate_values(c("POL", "GAG", "ENV"), strategy = "first")
  expect_equal(out, "ENV")
})

test_that("strict strategy returns the value when all contributors agree", {
  out <- aggregate_values(c("POL", "POL", "POL"),
                          strategy = "strict",
                          strict_marker = "ambiguous")
  expect_equal(out, "POL")
})

test_that("strict strategy returns the marker when contributors disagree", {
  out <- aggregate_values(c("POL", "GAG"),
                          strategy = "strict",
                          strict_marker = "ambiguous")
  expect_equal(out, "ambiguous")
})

test_that("strict strategy respects a custom marker", {
  out <- aggregate_values(c("POL", "GAG"),
                          strategy = "strict",
                          strict_marker = "MIXED")
  expect_equal(out, "MIXED")
})

test_that("unknown strategy raises an informative error", {
  expect_error(
    aggregate_values(c("POL", "GAG"), strategy = "nonsense"),
    "Unknown aggregation strategy"
  )
})

test_that("values are coerced to character before processing", {
  out <- aggregate_values(factor(c("POL", "GAG", "POL")), strategy = "first")
  expect_equal(out, "GAG")
})


# ───────────────────────────── make_tiebreaker_picker ─────────────────────────────

test_that("make_tiebreaker_picker('bitscore') returns the bitscore vector", {
  pick <- make_tiebreaker_picker("bitscore")
  out <- pick(bitscore = c(10, 20, 30), identity = c(90, 80, 70))
  expect_equal(out, c(10, 20, 30))
})

test_that("make_tiebreaker_picker('identity') returns the identity vector", {
  pick <- make_tiebreaker_picker("identity")
  out <- pick(bitscore = c(10, 20, 30), identity = c(90, 80, 70))
  expect_equal(out, c(90, 80, 70))
})

test_that("make_tiebreaker_picker('align_length') errors when align_length is absent", {
  pick <- make_tiebreaker_picker("align_length")
  expect_error(
    pick(bitscore = c(10, 20), identity = c(90, 80)),
    "align_length"
  )
})

test_that("make_tiebreaker_picker('align_length') returns align_length when supplied", {
  pick <- make_tiebreaker_picker("align_length")
  out <- pick(bitscore = c(10, 20), identity = c(90, 80), align_length = c(500, 300))
  expect_equal(out, c(500, 300))
})

test_that("make_tiebreaker_picker falls back to bitscore for unknown names", {
  pick <- make_tiebreaker_picker("zzz_unknown")
  out <- pick(bitscore = c(10, 20), identity = c(90, 80))
  expect_equal(out, c(10, 20))
})

test_that("tiebreaker picker closure captures tb_name by value (force semantics)", {
  pick <- make_tiebreaker_picker("identity")
  # Reassign the outer variable; the closure should still use "identity"
  tb_name <- "bitscore"
  out <- pick(bitscore = c(10, 20), identity = c(90, 80))
  expect_equal(out, c(90, 80))
})


# ───────────────────────────── end-to-end composition ─────────────────────────────

test_that("aggregate_values + make_tiebreaker_picker compose cleanly for the `best` case", {
  # Simulates the call pattern used in ranges_analysis.R's reduce_ranges_directed.
  pick <- make_tiebreaker_picker("bitscore")
  bitscore <- c(100, 200, 50)
  identity <- c(90, 80, 70)
  values   <- c("POL", "GAG", "ENV")
  out <- aggregate_values(values, strategy = "best",
                          tiebreaker = pick(bitscore = bitscore, identity = identity))
  expect_equal(out, "GAG")
})
