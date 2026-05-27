# testthat coverage for workflow/scripts/pair_detector/pairing.R
#
# Pure probe-pairing helpers extracted from pair_detector.R. Sources only the
# leaf module; builds synthetic GRanges so no GFF3 I/O is needed.

suppressMessages({
  library(testthat)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(dplyr)
  library(tibble)
})

.script_dir <- file.path("..", "..", "scripts", "pair_detector")
source(file.path(.script_dir, "pairing.R"))


# Synthetic validated-hit GRanges (all + strand to keep findOverlaps strand-clean):
#   ENV  @ 1000-1100  (the self probe)
#   GAG  @ 1200-1300  (gap 99 from ENV  -> within a 200 bp max_gap)
#   POL  @ 1100-1150  (abuts ENV        -> within any max_gap)
#   GAG  @ 5000-5100  (gap ~3.9 kb      -> outside a 200 bp max_gap)
.fake_ranges <- function() {
  GRanges(
    seqnames = "chr1",
    ranges = IRanges(
      start = c(1000, 1200, 1100, 5000),
      end   = c(1100, 1300, 1150, 5100)
    ),
    strand = "+",
    probe  = c("ENV", "GAG", "POL", "GAG"),
    label  = "Lentivirus",
    virus  = "HIV"
  )
}


test_that("find_pairs returns one row per other-probe within max_gap", {
  pairs <- find_pairs(.fake_ranges(), probe_to_pair = "ENV", max_gap = 200L)

  expect_equal(nrow(pairs), 2L)
  expect_setequal(pairs$other.probe, c("GAG", "POL"))
  expect_true(all(pairs$self.probe == "ENV"))
})

test_that("find_pairs excludes other probes beyond max_gap", {
  pairs <- find_pairs(.fake_ranges(), probe_to_pair = "ENV", max_gap = 200L)
  # the distant GAG @ 5000 must not appear
  expect_false(5000 %in% pairs$other.start)
})

test_that("find_pairs widening max_gap captures the distant pair", {
  pairs <- find_pairs(.fake_ranges(), probe_to_pair = "ENV", max_gap = 5000L)

  expect_equal(nrow(pairs), 3L)
  expect_equal(sum(pairs$other.start == 5000), 1L)
})

test_that("find_pairs computes coord / bio / span distances", {
  pairs <- find_pairs(.fake_ranges(), probe_to_pair = "ENV", max_gap = 200L)
  expect_true(all(c("coord_distance", "bio_distance", "span_distance") %in% names(pairs)))

  # GAG @ 1200 vs ENV @ 1000: coord_distance = |1000 - 1200| = 200
  gag <- pairs[pairs$other.probe == "GAG", ]
  expect_equal(gag$coord_distance, 200)
})

test_that("get_5prime is strand-aware", {
  expect_equal(get_5prime(100, 200, "+"), 100)
  expect_equal(get_5prime(100, 200, "-"), 200)
  expect_equal(get_5prime(100, 200, "*"), 100)
})
