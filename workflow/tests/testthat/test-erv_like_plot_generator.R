# testthat tests for the ERV-like structural plot builders
# (workflow/scripts/taxonomy/erv_like_plot_generator.R).
#
# The panel now reads the taxon-founded loci table produced by taxonomy_classify
# (one row per LTR-element locus). Builder contract: each returns a ggplot and
# falls back to empty_plot() on zero-row input. We assert the contract plus that
# load_taxon_loci coerces the string-typed structural columns.
#
# Run with: make test-r.

suppressMessages({
  library(testthat)
})

# Sourcing the generator loads tidyverse + shared helpers + the builders (the
# bottom-of-file main() guard keeps the CLI from firing under source()).
source("../../scripts/taxonomy/erv_like_plot_generator.R")

# A minimal taxon-founded loci frame, mirroring the classifier's LOCI_COLUMNS
# (structural columns arrive as strings from parquet; load_taxon_loci coerces).
.loci <- function(n = 2L) {
  tibble::tibble(
    id = c("L0", "L1")[seq_len(n)], seqname = "chr1",
    start = c(100, 5000)[seq_len(n)], end = c(2000, 7000)[seq_len(n)],
    strand = "+", parent = c("retro1", "retro2")[seq_len(n)],
    genes_present = c("ENV,GAG,POL", "GAG,POL")[seq_len(n)],
    n_main_genes = c(3L, 2L)[seq_len(n)],
    completeness = c(1.0, 2 / 3)[seq_len(n)],
    canonical_order = c(TRUE, FALSE)[seq_len(n)],
    structure_class = c("full", "partial")[seq_len(n)],
    domain_tier = c("domain_selected", "domain_unlisted")[seq_len(n)],
    taxon_call = c("Gammaretrovirus", "Betaretrovirus")[seq_len(n)],
    rank = "genus", resolved = "True", confidence = 1, method = "lca",
    is_mosaic = FALSE, mosaic_composition = "", erv_class = "Class I",
    species = c("Species_A", "Species_B")[seq_len(n)],
    span_bp = c(1901, 2001)[seq_len(n)]
  )
}

.empty_loci <- function() .loci()[0, ]

# Every builder returns a ggplot on non-empty input.
test_that("structural builders return ggplots", {
  loci <- .loci()
  expect_s3_class(completeness_plot(loci), "ggplot")
  expect_s3_class(canonical_order_plot(loci), "ggplot")
  expect_s3_class(gene_combinations_plot(loci), "ggplot")
  expect_s3_class(length_distribution_plot(loci), "ggplot")
  expect_s3_class(n_main_genes_plot(loci), "ggplot")
  expect_s3_class(composition_heatmap_plot(loci), "ggplot")
  expect_s3_class(structure_class_plot(loci), "ggplot")
})

# Every builder falls back to the labelled placeholder on zero-row input.
test_that("builders fall back to empty_plot on zero-row input", {
  el <- .empty_loci()
  expect_match(completeness_plot(el)$labels$title, "No loci")
  expect_match(canonical_order_plot(el)$labels$title, "No loci")
  expect_match(gene_combinations_plot(el)$labels$title, "No loci")
  expect_match(length_distribution_plot(el)$labels$title, "No loci")
  expect_match(n_main_genes_plot(el)$labels$title, "No loci")
  expect_match(composition_heatmap_plot(el)$labels$title, "No taxon-resolved loci")
  expect_match(structure_class_plot(el)$labels$title, "No loci")
})

# load_taxon_loci coerces the classifier's string-typed structural columns.
test_that("load_taxon_loci coerces structural columns + derives span_bp", {
  tmp <- tempfile(fileext = "")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))
  df <- tibble::tibble(
    id = "L0", seqname = "chr1", start = "100", end = "600", strand = "+",
    parent = "retro1", genes_present = "GAG,POL", n_main_genes = "2",
    completeness = "0.667", canonical_order = "True",
    taxon_call = "Gammaretrovirus", rank = "genus", resolved = "True",
    confidence = "1.000",
    method = "lca", per_gene = "", is_mosaic = "False",
    mosaic_composition = "", erv_class = "Class I", probe_label_set = "",
    ref_version = "abc"
  )
  arrow::write_parquet(df, file.path(tmp, "Species_A.loci.parquet"))
  loaded <- load_taxon_loci(tmp)
  expect_equal(loaded$species, "Species_A")
  expect_type(loaded$completeness, "double")
  expect_type(loaded$n_main_genes, "integer")
  expect_true(is.logical(loaded$canonical_order))
  expect_equal(loaded$span_bp, 501)            # that is 600 - 100 + 1
})


# ---------------------------------------------------------------------------
# Species labels must be the readable config `species:` names. Regression guard:
# this generator derives `species` from the parquet filename (the genome stem),
# so without the conversion every structural page showed `mus_musculus`.
# ---------------------------------------------------------------------------
test_that("erv_like species stems become readable names", {
  species_map <- list(mus_musculus = "Mus musculus",
                      HLmyoMyo6    = "Myotis myotis")
  stems <- c("mus_musculus", "HLmyoMyo6", "not_in_config")
  expect_equal(
    display_species(stems, species_map),
    c("Mus musculus", "Myotis myotis", "not in config")   # unmapped: stem, spaced
  )
})

test_that("per-host pages follow the configured species order", {
  ctx <- list(species_order = c("Species_B", "Species_A"))
  p <- structure_class_plot(.loci(), ctx)
  expect_equal(p$scales$get_scales("x")$limits, c("Species_A", "Species_B"))
  hist <- completeness_plot(.loci(), ctx)
  expect_equal(levels(hist$data$species), c("Species_B", "Species_A"))
})

test_that("rare gene combinations are pooled into one Other bar", {
  loci <- tibble::tibble(genes_present = c(rep("GAG,POL", 3), sprintf("G%02d", 1:40)))
  p <- gene_combinations_plot(loci)
  expect_equal(nrow(p$data), 31L)
  expect_true("Other (11 combinations)" %in% as.character(p$data$genes_present))
})
