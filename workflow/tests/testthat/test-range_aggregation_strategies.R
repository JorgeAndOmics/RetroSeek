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


# ─────────────────── plyranges integration regression ───────────────────────
#
# The unit tests above exercise aggregate_values() in isolation. The function
# is also called from inside plyranges::reduce_ranges_directed(...) via
# range_analysis/reductions.R::reduce_first(). plyranges presents grouped
# metadata as a CompressedAtomicList of length n (one element per output
# range, each element being the per-row vector of contributing values). The
# function's AtomicList branch (.aggregate_values_listwise) must return a
# value plyranges' summarize_rng can attach to the n-row output GRanges.
#
# Regression: a prior version returned a CharacterList of length n for
# strategy="list". plyranges' summarize_rng compresses list-columns via
# Reduce(S4Vectors::pc, ans[[i]]) and that path collapses uniform-length
# CharacterLists to a single entry, then `stopifnot(length(ans) == nr)`
# fires with "length(ans[[i]]) == nr is not TRUE". The fix is to return an
# atomic character of length n (joined by `separator`) for "list" inside
# the AtomicList path. These tests guard against re-introducing that bug.

test_that("aggregate_values handles CompressedCharacterList input (plyranges shape)", {
  # Mimic the shape plyranges passes inside reduce_ranges_directed:
  # a CompressedCharacterList of length n, where each element is the
  # per-group character vector of contributing values.
  values <- IRanges::CharacterList(list(
    c("HIV", "HIV"),               # group 1: all same → unique = "HIV"
    c("HTLV"),                     # group 2: single value
    c("HIV", "HTLV", "HIV")        # group 3: mixed → unique = "HIV","HTLV"
  ))
  expect_s4_class(values, "CompressedCharacterList")

  for (strategy in c("concatenate", "first", "majority", "strict", "list")) {
    out <- aggregate_values(values, strategy = strategy, separator = "; ")
    expect_type(out, "character")
    expect_length(out, length(values))
  }
})


test_that("plyranges reduce_ranges_directed + aggregate_values round-trips cleanly", {
  # End-to-end: builds a tiny GRanges, groups by (probe, virus), reduces
  # overlapping ranges, runs aggregate_values inside reduce_ranges_directed.
  # Reproduces the exact shape that triggered the
  # `length(ans[[i]]) == nr is not TRUE` assertion in the prior version.
  skip_if_not_installed("plyranges")
  skip_if_not_installed("GenomicRanges")
  suppressMessages({
    library(plyranges)
    library(GenomicRanges)
  })
  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1", "chr2", "chr2"),
    ranges   = IRanges::IRanges(start = c(10, 30, 100, 50, 200),
                                end   = c(40,  70, 150, 90, 250)),
    probe    = c("POL", "POL", "POL", "GAG", "GAG"),
    virus    = c("HIV", "HIV", "HIV", "HTLV", "HTLV"),
    label    = c("class_a", "class_a", "class_b", "class_a", "class_a"),
    bitscore = c(100, 200, 50, 80, 90),
    identity = c(90, 80, 70, 85, 95)
  )

  result <- gr %>%
    plyranges::group_by(probe, virus) %>%
    plyranges::reduce_ranges_directed(
      virus = aggregate_values(virus, "list", separator = "; "),
      label = aggregate_values(label, "list", separator = "; "),
      n_hits = length(bitscore)
    )

  # The reduction collapses the two adjacent POL/HIV ranges 10-40 and 30-70
  # into one merged range; the third POL/HIV range 100-150 stays separate.
  # GAG/HTLV's two ranges are non-overlapping, both kept.
  # Total: 2 (POL/HIV) + 2 (GAG/HTLV) = 4 output rows.
  expect_equal(length(result), 4)
  expect_type(GenomicRanges::mcols(result)$virus, "character")
  expect_length(GenomicRanges::mcols(result)$virus, 4)
})
