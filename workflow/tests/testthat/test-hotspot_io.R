# testthat tests for hotspot_analysis/io.R + utils/chrom_names.R
#
# Run via: make test-r

suppressMessages({
  library(testthat)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})

source("../../scripts/utils/chrom_names.R")
source("../../scripts/hotspot/io.R")


# ------------------------ normalise_chrom_names -------------------------

test_that("normalise_chrom_names extracts NCBI accession tokens from full headers", {
  headers <- c(
    "CM034567.1 Genus species chromosome 1, GRCh38",
    "NC_000001.11 Homo sapiens chromosome 1",
    "JH791234.1 unplaced scaffold"
  )
  expect_equal(
    suppressMessages(normalise_chrom_names(headers)),
    c("CM034567.1", "NC_000001.11", "JH791234.1")
  )
})

test_that("normalise_chrom_names returns NA for non-matching headers and emits a message", {
  headers <- c("CM034567.1 ok", "chr1 custom-format-no-match")
  expect_message(
    out <- normalise_chrom_names(headers),
    regexp = "did not match pattern"
  )
  expect_equal(out, c("CM034567.1", NA_character_))
})

test_that("normalise_chrom_names is silent when all headers match", {
  headers <- c("CM034567.1 a", "CM034568.1 b")
  expect_silent(normalise_chrom_names(headers))
})


# --------------------------- read_hotspot_options ---------------------------

test_that("read_hotspot_options returns documented defaults when the config is empty", {
  config <- list(parameters = list(), hotspot = list())
  opts <- read_hotspot_options(config)
  expect_equal(opts$seed, 67L)
  # ADR-012: the catalog tier (per integration event) is the default, grouped by
  # the calibrated lineage call; `valid` and the group_split bool are retired.
  expect_equal(opts$input, "catalog")
  expect_equal(opts$group_by, "segment")
  expect_equal(opts$source, "ltr-flanked")
  expect_equal(opts$window_size, 500000L)
  expect_equal(opts$mask_size, 20L)
  expect_equal(opts$mask_mismatch, 3L)
  expect_equal(opts$pvalue_threshold, 0.05)
  expect_equal(opts$min_hits, 2L)
  expect_equal(opts$merge_gap, 0L)
  expect_true(opts$strata_by_chromosome)
  expect_equal(opts$unplaced_min_factor, 10L)
})

test_that("read_hotspot_options reads operational knobs from the top-level hotspot section", {
  config <- list(
    parameters = list(seed = 1234),
    hotspot = list(
      input                = "original",
      window_size          = 5000,
      merge_gap            = 100,
      strata_by_chromosome = FALSE
    )
  )
  opts <- read_hotspot_options(config)
  expect_equal(opts$seed, 1234L)            # seed still comes from parameters
  expect_equal(opts$input, "original")
  expect_equal(opts$window_size, 5000L)
  expect_equal(opts$merge_gap, 100L)
  expect_false(opts$strata_by_chromosome)
})


# --------------------------- assert_hits_on_genome --------------------------

test_that("assert_hits_on_genome aborts on a GenBank-vs-RefSeq accession mismatch", {
  hits <- GenomicRanges::GRanges(
    seqnames = c("NC_071387.1", "NC_071387.1"),
    ranges   = IRanges::IRanges(start = c(100, 500), end = c(200, 600))
  )
  seqlengths <- c("CM040288.1" = 1e6, "CM040289.1" = 2e6)
  expect_error(assert_hits_on_genome(hits, seqlengths),
               regexp = "map to a genome contig")
})

test_that("assert_hits_on_genome passes (frac == 1) when hit seqnames match the genome", {
  hits <- GenomicRanges::GRanges(
    seqnames = c("CM040288.1", "CM040289.1"),
    ranges   = IRanges::IRanges(start = c(100, 500), end = c(200, 600))
  )
  seqlengths <- c("CM040288.1" = 1e6, "CM040289.1" = 2e6)
  expect_equal(assert_hits_on_genome(hits, seqlengths), 1)
})


# --------------------------- load_hits_gff (label assertion) ---------------

test_that("load_hits_gff aborts when the GFF lacks an mcols$label column", {
  # Build a minimal GFF3 file via rtracklayer::export with no label column
  gr <- GenomicRanges::GRanges(
    seqnames = "CM000001.1",
    ranges   = IRanges::IRanges(start = 100, end = 200)
  )
  S4Vectors::mcols(gr)$type <- "feature"  # no `label`
  tmp <- tempfile(fileext = ".gff3")
  on.exit(unlink(tmp), add = TRUE)
  rtracklayer::export(gr, tmp, format = "gff3")
  expect_error(load_hits_gff(tmp), regexp = "label")
})

test_that("load_hits_gff returns the imported GRanges when label is present", {
  gr <- GenomicRanges::GRanges(
    seqnames = "CM000001.1",
    ranges   = IRanges::IRanges(start = 100, end = 200)
  )
  S4Vectors::mcols(gr)$label <- "Gammaretrovirus"
  tmp <- tempfile(fileext = ".gff3")
  on.exit(unlink(tmp), add = TRUE)
  rtracklayer::export(gr, tmp, format = "gff3")
  out <- load_hits_gff(tmp)
  expect_s4_class(out, "GRanges")
  expect_true("label" %in% colnames(S4Vectors::mcols(out)))
})


# ------------------------ load_catalog_loci (ADR-012) -------------------------
# The catalog is per-LOCUS: one row per integration event. These pin the species
# matching (stem vs display name), the tier filter, and the fail-loud contract.

.write_catalog <- function(dir) {
  path <- file.path(dir, "catalog.csv")
  readr::write_csv(tibble::tibble(
    species         = c("Mus musculus", "Mus musculus", "Mus musculus",
                        "Homo sapiens"),
    source          = c("ltr-flanked", "orphan", "ltr-flanked", "ltr-flanked"),
    seqname         = c("chr1", "chr1", "chr2", "chr1"),
    start           = c(100L, 500L, 100L, 100L),
    end             = c(200L, 600L, 200L, 200L),
    strand          = c("+", "-", "+", "+"),
    taxon_call      = c("Gammaretrovirus", "Betaretrovirus", "Gammaretrovirus",
                        "Betaretrovirus"),
    segment         = c("Gammaretrovirus", "Betaretrovirus", "Gammaretrovirus",
                        "Betaretrovirus"),
    structure_class = c("full", "gene", "partial", "gene"),
    confidence      = c("1.000", "0.400", "0.900", "0.800"),
    confidence_tag  = c("HC", "LC", "HC", "HC")
  ), path)
  path
}

test_that("load_catalog_loci matches a genome stem to its display name", {
  tmp <- tempfile(); dir.create(tmp)
  path <- .write_catalog(tmp)
  gr <- load_catalog_loci(path, "Mus_musculus",
                          list(Mus_musculus = "Mus musculus"), "both")
  expect_equal(length(gr), 3L)                      # not the Homo row
  expect_true(all(c("taxon_call", "segment", "structure_class") %in%
                    colnames(S4Vectors::mcols(gr))))
})

test_that("load_catalog_loci filters to the requested tier", {
  tmp <- tempfile(); dir.create(tmp)
  path <- .write_catalog(tmp)
  gr <- load_catalog_loci(path, "Mus_musculus",
                          list(Mus_musculus = "Mus musculus"), "ltr-flanked")
  expect_equal(length(gr), 2L)
  expect_true(all(S4Vectors::mcols(gr)$source == "ltr-flanked"))
})

test_that("load_catalog_loci counts one row per locus, not per gene", {
  # The whole point of ADR-012: a multi-gene provirus is ONE event here.
  tmp <- tempfile(); dir.create(tmp)
  gr <- load_catalog_loci(.write_catalog(tmp), "Mus_musculus",
                          list(Mus_musculus = "Mus musculus"), "both")
  expect_equal(length(gr), nrow(unique(as.data.frame(gr)[, c("seqnames", "start")])))
})

test_that("load_catalog_loci fails loud when a genome matches no rows", {
  tmp <- tempfile(); dir.create(tmp)
  path <- .write_catalog(tmp)
  # Silence here would be indistinguishable from 'this genome has no ERVs'.
  expect_error(load_catalog_loci(path, "Gallus_gallus", NULL, "both"),
               "No catalog rows")
})

test_that("load_catalog_loci fails loud when the tier is empty", {
  tmp <- tempfile(); dir.create(tmp)
  path <- .write_catalog(tmp)
  expect_error(
    load_catalog_loci(path, "Homo_sapiens", list(Homo_sapiens = "Homo sapiens"),
                      "orphan"),
    "No 'orphan' loci"
  )
})

test_that("read_hotspot_options defaults to the catalog tier and segment grouping", {
  opts <- read_hotspot_options(list(hotspot = list(), parameters = list(seed = 1L)))
  expect_equal(opts$input, "catalog")
  expect_equal(opts$group_by, "segment")
  expect_equal(opts$source, "ltr-flanked")
})
