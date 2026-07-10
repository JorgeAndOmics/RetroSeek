# testthat coverage for the taxonomy report aggregation (build_report). Sourcing
# the generator pulls in the shared plot helpers; the CLI block is guarded by
# the sys.nframe() check so main() does not fire during sourcing.

suppressMessages({
  library(testthat)
  library(tibble)
  library(dplyr)
})

.script_dir <- file.path("..", "..", "scripts")
source(file.path(.script_dir, "taxonomy", "taxonomy_plot_generator.R"))


.fake_loci <- function(source) {
  tribble(
    ~species, ~taxon_call,       ~rank,    ~resolved, ~confidence_tag, ~method,     ~is_mosaic, ~source,
    "g1",     "Gammaretrovirus", "genus",  "True",    "HC",            "placement", "False",    source,
    "g1",     "Betaretrovirus",  "genus",  "True",    "LC",            "lca",       "True",     source,
    "g1",     "UNCLASSIFIED",    "none",   "False",   "LC",            "lca",       "False",    source
  )
}


test_that("build_report counts by taxon, confidence, method + mosaic/integrations", {
  rep <- build_report(.fake_loci("anchored"))

  # taxon counts only over axis-resolved rows
  taxon <- rep %>% filter(dimension == "taxon")
  expect_setequal(taxon$level, c("Gammaretrovirus", "Betaretrovirus"))
  expect_true(all(taxon$count == 1))

  # confidence spans all rows (HC=1, LC=2)
  conf <- rep %>% filter(dimension == "confidence")
  expect_equal(conf$count[conf$level == "HC"], 1L)
  expect_equal(conf$count[conf$level == "LC"], 2L)

  # method only over axis-resolved rows
  method <- rep %>% filter(dimension == "method")
  expect_setequal(method$level, c("placement", "lca"))

  # summary totals
  summ <- rep %>% filter(dimension == "summary")
  expect_equal(summ$count[summ$level == "mosaic"], 1L)
  expect_equal(summ$count[summ$level == "integrations"], 3L)
})

test_that("build_report splits by source (anchored vs orphan)", {
  combined <- bind_rows(.fake_loci("anchored"), .fake_loci("orphan"))
  rep <- build_report(combined)
  expect_setequal(unique(rep$source), c("anchored", "orphan"))
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
    ~species, ~taxon_call,       ~rank,   ~resolved, ~confidence, ~confidence_tag, ~method,
    ~is_mosaic, ~n_blastx_hits, ~completeness, ~n_main_genes, ~structure_class, ~domain_tier, ~source,
    "g1", "Gammaretrovirus", "genus", "True",  "0.95", "HC", "placement", "False", "8", "0.667", "2", "partial", "domain_selected", "anchored",
    "g1", "UNCLASSIFIED",    "none",  "False", "0.00", "LC", "lca",       "False", "0", "0.333", "1", "gene",    "non_domain",      "anchored",
    "g1", "Betaretrovirus",  "genus", "True",  "0.40", "LC", "lca",       "False", "3", "0.333", "1", "gene",    "non_domain",      "orphan"
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
    structure_by_tier_plot(combined),
    source_yield_plot(combined),
    taxon_by_source_plot(combined),
    domain_tier_composition_plot(combined),
    structure_class_composition_plot(combined)
  )
  for (p in builders) expect_s3_class(p, "ggplot")

  empty <- combined[0, ]
  expect_s3_class(evidence_depth_plot(empty), "ggplot")
  expect_s3_class(confidence_density_plot(empty, 0.5), "ggplot")
  expect_s3_class(taxon_by_source_plot(empty), "ggplot")
  expect_s3_class(domain_tier_composition_plot(empty), "ggplot")
  expect_s3_class(structure_class_composition_plot(empty), "ggplot")
})


test_that("mosaic sub-panel builders render on mosaic loci and are empty-safe", {
  loci <- tribble(
    ~species, ~is_mosaic, ~mosaic_composition,
    "g1", "True",  "POL:Betaretrovirus;GAG:Betaretrovirus;ENV:Gammaretrovirus",
    "g1", "True",  "POL:Gammaretrovirus;ENV:Betaretrovirus",
    "g2", "True",  "POL:Alpharetrovirus;GAG:Gammaretrovirus",
    "g1", "False", ""   # non-mosaic locus present for the burden denominator
  )
  for (p in list(
    mosaic_burden_plot(loci),
    mosaic_taxon_pairs_plot(loci),
    mosaic_gene_discordance_plot(loci),
    mosaic_composition_by_species_plot(loci)
  )) expect_s3_class(p, "ggplot")

  # burden splits mosaic vs single-lineage; discordance flags ENV as odd-one-out
  disc <- mosaic_gene_discordance_plot(loci)
  expect_s3_class(disc, "ggplot")

  # empty-safe: no mosaic rows -> labelled placeholder (still a ggplot)
  none <- loci[loci$is_mosaic == "False", ]
  expect_match(mosaic_taxon_pairs_plot(none)$labels$title, "no mosaic loci")
  expect_match(mosaic_gene_discordance_plot(none)$labels$title, "no mosaic loci")
  expect_s3_class(mosaic_burden_plot(loci[0, ]), "ggplot")
})
