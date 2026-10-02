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
  mcols(gr) <- data.frame(label = "Betaretrovirus", count = 12L,
                          hotspot_id = "G_HS_00001",
                          n_loci = 14L, n_full = 3L, n_partial = 5L, n_gene = 6L)
  gr
}

.seqlengths <- c(chr1 = 2e6, chr2 = 2e6)

test_that("every hotspot page builds from real input", {
  expect_s3_class(plot_manhattan(.windows(), 0.05, "Mus musculus", "Betaretrovirus"),
                  "gg")
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


# ---------------------------------------------------------------------------
# Small scaffolds: tested, but not drawn
# ---------------------------------------------------------------------------
# The Homo sapiens assembly has 709 sequences; drawing every one turned the
# karyotype and the Manhattan axis into an unreadable smear of accessions.
test_that("only sequences of model size are drawn, in their own order", {
  sl <- c(chr2 = 8e6, scaf1 = 1e5, chr1 = 9e6, scaf2 = 2e5)
  expect_equal(sequences_to_draw(sl, min_length = 5e6), c("chr2", "chr1"))
})

test_that("past the cap, the longest sequences are drawn", {
  sl <- c(a = 6e6, b = 9e6, c = 7e6)
  expect_equal(sequences_to_draw(sl, min_length = 5e6, max_n = 2L), c("b", "c"))
})

test_that("an assembly of small scaffolds still draws its longest ones", {
  sl <- c(a = 1e5, b = 3e5, c = 2e5)
  expect_equal(sequences_to_draw(sl, min_length = 5e6, max_n = 2L), c("b", "c"))
})

test_that("the karyotype and the Manhattan page say what they leave out", {
  sl <- c(chr1 = 2e6, chr2 = 2e6, scaf9 = 1e4)
  kary <- plot_karyotype(sl, .hotspots(), "Mus musculus", draw = c("chr1", "chr2"))
  expect_equal(levels(kary$layers[[1]]$data$chrom), c("chr1", "chr2"))
  expect_match(kary$labels$caption, "1 shorter sequence is tested but not drawn")
  manh <- plot_manhattan(.windows(), 0.05, "Mus musculus", "Betaretrovirus",
                         draw = "chr1")
  expect_setequal(unique(manh$data$chrom), "chr1")
  expect_match(manh$labels$caption, "1 shorter sequence is tested but not drawn")
})

test_that("a lineage is named in words, not as a table value", {
  p <- plot_manhattan(.windows(), 0.05, "Mus musculus", "unassigned_at_genus")
  expect_match(.subtitle_text(p), "Unassigned at genus loci per window")
  expect_no_match(.subtitle_text(p), "unassigned_at_genus")
})


# ---------------------------------------------------------------------------
# Q-Q with randomised p-values; the whole-chromosome page
# ---------------------------------------------------------------------------
source(file.path(.dir, "hotspot", "chromosomes.R"))

test_that("with the model's theta the Q-Q page draws randomised p-values", {
  w <- .windows()
  w$mu_nb <- 0.5
  w$effective_bp <- 500000
  w$count <- c(3L, 0L, 0L, 1L, 0L, 2L, 0L, 0L)
  p <- plot_qq(w, "Mus musculus", "Betaretrovirus", theta = 2)
  expect_match(.subtitle_text(p), "randomised")
  expect_equal(nrow(p$data), 8L)
})

test_that("the chromosome page shows each tested chromosome, and none when untested", {
  d <- data.frame(label = "Betaretrovirus", chrom = c("chr1", "chr2"),
                  observed = c(10L, 40L), expected = c(25, 25),
                  rate_ratio = c(0.4, 1.6), pval = c(0.9, 0.01), qval = c(0.9, 0.02))
  p <- plot_chromosome_rates(d, "Mus musculus")
  expect_s3_class(p, "ggplot")
  expect_null(plot_chromosome_rates(d[0, ], "Mus musculus"))
})

test_that("a chromosome with no loci draws without a warning", {
  d <- data.frame(label = "Betaretrovirus", chrom = c("chr1", "chr2"),
                  observed = c(0L, 50L), expected = c(25, 25),
                  rate_ratio = c(0, 2), pval = c(1, 0.001), qval = c(1, 0.002))
  expect_no_warning(ggplot2::ggplot_build(plot_chromosome_rates(d, "Mus musculus")))
})

test_that("the chromosome page keeps lineages with a call and names chromosomes", {
  d <- data.frame(label = rep(c("Betaretrovirus", "Gammaretrovirus"), each = 2),
                  chrom = c("CM1.1", "CM2.1", "CM1.1", "CM2.1"),
                  observed = c(10L, 40L, 20L, 20L), expected = 25,
                  rate_ratio = c(0.4, 1.6, 0.8, 0.8), pval = c(0.9, 0.01, 0.7, 0.7),
                  qval = c(0.9, 0.02, 0.7, 0.7))
  p <- plot_chromosome_rates(d, "Mus musculus",
                             chrom_labels = c(CM1.1 = "chr1", CM2.1 = "chrX"))
  expect_setequal(unique(as.character(p$data$label)), "Betaretrovirus")
  expect_setequal(as.character(p$data$chrom), c("chr1", "chrX"))
  quiet <- plot_chromosome_rates(transform(d, qval = 0.5), "Mus musculus")
  expect_match(quiet$labels$title, "No chromosome")
})
