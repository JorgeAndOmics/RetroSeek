# testthat coverage for workflow/scripts/species_segmenter/segment.R
#
# Pure partition + species-list helpers extracted from species_segmenter.R.
# Sources only the leaf module so we avoid arrow / yaml I/O machinery.

suppressMessages({
  library(testthat)
  library(tibble)
})

.script_dir <- file.path("..", "..", "scripts", "species_segmenter")
source(file.path(.script_dir, "segment.R"))


# 5-row synthetic hit table: 3 main-probe rows (POL/GAG/ENV), 2 accessory.
.fake_hits <- function() {
  tibble::tibble(
    species = c("S1", "S1", "S2", "S2", "S3"),
    probe   = c("POL", "GAG", "ENV", "VIF", "TAT"),
    virus   = c("HIV", "HIV", "HTLV", "HTLV", "FFV")
  )
}


test_that("segment_by_probe splits main vs accessory by probe membership", {
  out <- segment_by_probe(.fake_hits(), main_probe_names = c("POL", "GAG", "ENV"))

  expect_equal(nrow(out$main), 3L)
  expect_equal(nrow(out$accessory), 2L)
  expect_setequal(out$main$probe, c("POL", "GAG", "ENV"))
  expect_setequal(out$accessory$probe, c("VIF", "TAT"))
})

test_that("segment_by_probe is exhaustive and order-preserving", {
  df <- .fake_hits()
  out <- segment_by_probe(df, main_probe_names = c("POL", "GAG", "ENV"))

  # every row lands in exactly one partition
  expect_equal(nrow(out$main) + nrow(out$accessory), nrow(df))
  # accessory rows keep their original relative order
  expect_equal(out$accessory$probe, c("VIF", "TAT"))
})

test_that("segment_by_probe sends NA and unknown probes to accessory", {
  df <- tibble::tibble(species = c("S1", "S1"), probe = c(NA_character_, "POL"))
  out <- segment_by_probe(df, main_probe_names = "POL")

  expect_equal(out$main$probe, "POL")
  expect_equal(nrow(out$accessory), 1L)
  expect_true(is.na(out$accessory$probe))
})

test_that("parse_species_list trims whitespace and drops empty tokens", {
  expect_equal(parse_species_list("S1, S2 ,, S3 "), c("S1", "S2", "S3"))
  expect_equal(parse_species_list(""), character(0))
})

test_that("species_to_backfill returns configured species with no hits", {
  expect_equal(
    species_to_backfill(all_species = c("S1", "S2", "S3"), written = c("S1", "S3")),
    "S2"
  )
  expect_length(
    species_to_backfill(all_species = c("S1", "S2"), written = c("S1", "S2")),
    0L
  )
})
