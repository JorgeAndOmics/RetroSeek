# testthat tests for emit_hotspot_manifest() in workflow/scripts/hotspot/io.R
#
# Why this file exists: the emitter used to live inline in hotspot_detector.R,
# closing over four output paths and reading `args$gff`. ADR-012 renamed that
# CLI argument to `--hits` and updated the two load sites but not the emitter,
# so `file_md5(NULL)` returned NA and EVERY hotspot run wrote
# `path: ~ / md5: ~` for its input track - provenance silently lost, with no
# error. The orchestrator's own test is an unconditional skip(), so nothing
# caught it. These tests pin the contract that input paths and checksums are
# actually recorded.
#
# Run via: make test-r

suppressMessages({
  library(testthat)
  library(yaml)
})

source("../../scripts/ranges/exporters.R")   # file_md5()
source("../../scripts/utils/chrom_names.R")
source("../../scripts/hotspot/io.R")

# Minimal on-disk inputs so md5s are real rather than mocked.
.with_manifest <- function(code) {
  dir <- withr::local_tempdir()
  fasta <- file.path(dir, "genome.fa")
  hits <- file.path(dir, "hits.gff3")
  config <- file.path(dir, "config.yaml")
  writeLines(c(">chr1", "ACGT"), fasta)
  writeLines("##gff-version 3", hits)
  writeLines("hotspot: {}", config)
  path <- file.path(dir, "Test_species.manifest.yaml")
  emit_hotspot_manifest(
    inputs  = list(fasta = fasta, hits = hits, config = config),
    outputs = list(csv = "a.csv", parquet = "a.parquet", gff3 = "a.gff3", bed = "a.bed"),
    opts = list(window_size = 500000), species = "Test_species",
    species_name = "Test species", fit_diagnostics = list(),
    counts = list(total_hits = 3L), generator_version = "RetroSeek/test",
    path = path
  )
  code(yaml::read_yaml(path), list(fasta = fasta, hits = hits, config = config))
}

test_that("every input path is recorded, none left null", {
  .with_manifest(function(m, paths) {
    expect_setequal(names(m$inputs), c("fasta", "hits", "config"))
    for (key in names(m$inputs)) {
      expect_equal(m$inputs[[key]]$path, paths[[key]])
      expect_false(is.null(m$inputs[[key]]$path))
    }
  })
})

test_that("every input checksum is a real md5, not NA", {
  .with_manifest(function(m, paths) {
    for (key in names(m$inputs)) {
      md5 <- m$inputs[[key]]$md5
      expect_false(is.null(md5))
      expect_false(is.na(md5))
      expect_match(md5, "^[0-9a-f]{32}$")
    }
    # the regression itself: the hits track must carry the genome's sibling md5,
    # not the NA that `file_md5(args$gff)` produced
    expect_equal(m$inputs$hits$md5, unname(tools::md5sum(paths$hits)))
  })
})

test_that("outputs and provenance fields survive the round trip", {
  .with_manifest(function(m, paths) {
    expect_equal(m$generator, "RetroSeek/test")
    expect_equal(m$species, "Test_species")
    expect_setequal(names(m$outputs), c("csv", "parquet", "gff3", "bed"))
    expect_equal(m$counts$total_hits, 3L)
  })
})
