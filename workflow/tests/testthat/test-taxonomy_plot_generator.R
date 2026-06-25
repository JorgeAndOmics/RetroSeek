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
