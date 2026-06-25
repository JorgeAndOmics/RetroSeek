# testthat coverage for the taxonomy report aggregation (build_report). Sourcing
# the generator pulls in the shared plot helpers; the CLI block is guarded by
# the sys.nframe() check so main() does not fire during sourcing.

suppressMessages({
  library(testthat)
  library(tibble)
  library(dplyr)
})

.script_dir <- file.path("..", "..", "scripts")
source(file.path(.script_dir, "taxonomy_plot_generator.R"))


.fake_loci <- function(source) {
  tribble(
    ~species, ~genus_call,       ~rank,    ~confidence_tag, ~method,     ~is_mosaic, ~source,
    "g1",     "Gammaretrovirus", "genus",  "HC",            "placement", "False",    source,
    "g1",     "Betaretrovirus",  "genus",  "LC",            "lca",       "True",     source,
    "g1",     "UNCLASSIFIED",    "none",   "LC",            "lca",       "False",    source
  )
}


test_that("build_report counts by genus, confidence, method + mosaic/integrations", {
  rep <- build_report(.fake_loci("anchored"))

  # genus counts only over genus-rank rows
  genus <- rep %>% filter(dimension == "genus")
  expect_setequal(genus$level, c("Gammaretrovirus", "Betaretrovirus"))
  expect_true(all(genus$count == 1))

  # confidence spans all rows (HC=1, LC=2)
  conf <- rep %>% filter(dimension == "confidence")
  expect_equal(conf$count[conf$level == "HC"], 1L)
  expect_equal(conf$count[conf$level == "LC"], 2L)

  # method only over genus-rank rows
  method <- rep %>% filter(dimension == "method")
  expect_setequal(method$level, c("placement", "lca"))

  # summary totals
  summ <- rep %>% filter(dimension == "summary")
  expect_equal(summ$count[summ$level == "mosaic"], 1L)
  expect_equal(summ$count[summ$level == "integrations"], 3L)
})

test_that("build_report splits by source (anchored vs fragment)", {
  combined <- bind_rows(.fake_loci("anchored"), .fake_loci("fragment"))
  rep <- build_report(combined)
  expect_setequal(unique(rep$source), c("anchored", "fragment"))
  # each tier reports its own 3 integrations
  integ <- rep %>% filter(dimension == "summary", level == "integrations")
  expect_true(all(integ$count == 3))
})

test_that("build_report is empty-safe", {
  rep <- build_report(tibble())
  expect_equal(nrow(rep), 0L)
  expect_setequal(names(rep), c("source", "dimension", "level", "count"))
})


test_that("bucket_evidence maps hit counts to ordered buckets", {
  b <- bucket_evidence(c(0L, 1L, 2L, 5L, 6L, 25L))
  expect_equal(as.character(b), c("0", "1", "2–5", "2–5", "6+", "6+"))
  expect_true(is.factor(b))
  expect_equal(levels(b), c("0", "1", "2–5", "6+"))
})

test_that("bucket_evidence handles empty and is robust to numeric input", {
  expect_length(bucket_evidence(integer()), 0L)
  expect_equal(as.character(bucket_evidence(c(3))), "2–5")  # numeric, not integer
})


test_that("new plot builders return ggplot on data and empty_plot on empty", {
  combined <- tribble(
    ~species, ~genus_call,       ~rank,   ~confidence, ~confidence_tag, ~method,
    ~is_mosaic, ~n_blastx_hits, ~completeness, ~n_main_genes, ~source,
    "g1", "Gammaretrovirus", "genus", "0.95", "HC", "placement", "False", "8", "0.667", "2", "anchored",
    "g1", "UNCLASSIFIED",    "none",  "0.00", "LC", "lca",       "False", "0", "0.333", "1", "anchored",
    "g1", "Betaretrovirus",  "genus", "0.40", "LC", "lca",       "False", "3", "0.333", "1", "fragment"
  )
  # main() adds the numeric helper columns; mirror that here.
  combined <- dplyr::mutate(combined,
    confidence_num = as.numeric(.data$confidence),
    n_hits = as.integer(.data$n_blastx_hits),
    completeness_num = as.numeric(.data$completeness))

  builders <- list(
    evidence_depth_plot(combined),
    confidence_density_plot(combined, 0.5),
    confidence_vs_evidence_plot(combined),
    novel_structure_plot(combined),
    source_yield_plot(combined),
    genus_by_source_plot(combined)
  )
  for (p in builders) expect_s3_class(p, "ggplot")

  empty <- combined[0, ]
  expect_s3_class(evidence_depth_plot(empty), "ggplot")
  expect_s3_class(confidence_density_plot(empty, 0.5), "ggplot")
  expect_s3_class(genus_by_source_plot(empty), "ggplot")
})
