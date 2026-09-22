# testthat coverage for workflow/scripts/hotspot/plots.R
#
# The hotspot pages in the house style: every builder returns a page for real
# input and a labelled placeholder for empty input, and colours carry their
# fixed meanings (structure classes, lineages).
#
# Run via: make test-r

suppressMessages({
  library(testthat)
  library(ggplot2)
  library(GenomicRanges)
})

.dir <- file.path("..", "..", "scripts")
source(file.path(.dir, "plot2sort", "style.R"))
source(file.path(.dir, "plot2sort", "helpers.R"))
source(file.path(.dir, "hotspot", "plots.R"))

.windows <- function() {
  tibble::tibble(
    chrom = rep(c("chr1", "chr2"), each = 4),
    start = rep(c(1, 500001, 1000001, 1500001), 2),
    end   = start + 499999,
    label = "Betaretrovirus",
    pval_nb = c(1e-6, 0.2, 0.5, 0.9, 0.3, 0.04, 0.6, 0.8),
    qval_nb = c(1e-4, 0.4, 0.7, 0.9, 0.5, 0.1, 0.8, 0.9)
  )
}

.hotspots <- function() {
  gr <- GRanges("chr1", IRanges(1, 500000))
  mcols(gr) <- data.frame(label = "Betaretrovirus", count = 12L, hotspot_id = "G_HS_00001",
                          n_loci = 14L, n_full = 3L, n_partial = 5L, n_gene = 6L)
  gr
}

.seqlengths <- c(chr1 = 2e6, chr2 = 2e6)

test_that("every hotspot page builds from real input", {
  expect_s3_class(plot_manhattan(.windows(), 0.05, "Mus musculus", "Betaretrovirus"), "gg")
  expect_s3_class(plot_karyotype(.seqlengths, .hotspots(), "Mus musculus"), "gg")
  expect_s3_class(plot_qq(.windows(), "Mus musculus", "Betaretrovirus"), "gg")
  expect_s3_class(plot_summary_panel(.hotspots(), .seqlengths, "Mus musculus"), "gg")
  expect_s3_class(plot_hotspot_composition(.hotspots(), "Mus musculus"), "gg")
})

test_that("empty inputs give labelled placeholders", {
  expect_match(plot_manhattan(.windows()[0, ], 0.05)$labels$title, "No windows")
  expect_match(plot_summary_panel(GRanges(), .seqlengths)$labels$title, "No hotspot")
  expect_match(plot_hotspot_composition(GRanges())$labels$title, "No hotspot")
})

test_that("hotspots in the karyotype wear their lineage colour", {
  fills <- ggplot_build(plot_karyotype(.seqlengths, .hotspots()))$data[[2]]$fill
  expect_true(.GENUS_COLOUR[["Betaretrovirus"]] %in% fills)
})

test_that("composition uses the structure-class colours", {
  p <- plot_hotspot_composition(.hotspots(), "Mus musculus")
  fills <- unique(ggplot_build(p[[1]])$data[[1]]$fill)
  expect_setequal(fills, unname(.STRUCTURE_COLOUR))
})
