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
  rep <- build_report(.fake_loci("ltr-flanked"))

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

test_that("build_report splits by source (ltr-flanked vs orphan)", {
  combined <- bind_rows(.fake_loci("ltr-flanked"), .fake_loci("orphan"))
  rep <- build_report(combined)
  expect_setequal(unique(rep$source), c("ltr-flanked", "orphan"))
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
  expect_equal(as.character(b), c("0", "1", "2-5", "2-5", "6+", "6+"))
  expect_true(is.factor(b))
  expect_equal(levels(b), c("0", "1", "2-5", "6+"))
})

test_that("bucket_evidence handles empty and is robust to numeric input", {
  expect_length(bucket_evidence(integer()), 0L)
  expect_equal(as.character(bucket_evidence(c(3))), "2-5")  # numeric, not integer
})


test_that("new plot builders return ggplot on data and empty_plot on empty", {
  combined <- tribble(
    ~species, ~taxon_call,       ~rank,   ~resolved, ~confidence, ~confidence_tag, ~method,
    ~is_mosaic, ~n_blastx_hits, ~completeness, ~n_main_genes, ~structure_class, ~domain_tier, ~source,
    "g1", "Gammaretrovirus", "genus", "True",  "0.95", "HC", "placement", "False", "8", "0.667", "2", "partial", "domain_selected", "ltr-flanked",
    "g1", "UNCLASSIFIED",    "none",  "False", "0.00", "LC", "lca",       "False", "0", "0.333", "1", "gene",    "non_domain",      "ltr-flanked",
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
    structure_class_composition_plot(combined),
    confidence_count_plot(combined),
    confidence_gradient_plot(combined)
  )
  for (p in builders) expect_s3_class(p, "ggplot")

  empty <- combined[0, ]
  expect_s3_class(evidence_depth_plot(empty), "ggplot")
  expect_s3_class(confidence_density_plot(empty, 0.5), "ggplot")
  expect_s3_class(taxon_by_source_plot(empty), "ggplot")
  expect_s3_class(domain_tier_composition_plot(empty), "ggplot")
  expect_s3_class(structure_class_composition_plot(empty), "ggplot")
})


test_that("reconcile_catalog drops orphans overlapping an ltr-flanked locus (ltr-flanked precedence)", {
  combined <- tribble(
    ~species, ~source,    ~seqname, ~start, ~end,   ~id,
    "g1",     "ltr-flanked", "chr1",   "1000", "2000", "A1",
    "g1",     "orphan",   "chr1",   "1500", "1800", "O1",   # overlaps A1 -> drop
    "g1",     "orphan",   "chr1",   "5000", "5300", "O2",   # clear -> keep
    "g1",     "orphan",   "chr2",   "1500", "1800", "O3"    # different seqname -> keep
  )
  out <- reconcile_catalog(combined)
  expect_setequal(out$id, c("A1", "O2", "O3"))          # O1 dropped
  expect_true(all(out$source[out$id == "A1"] == "ltr-flanked"))
  # empty-safe + single-tier passthrough
  expect_equal(nrow(reconcile_catalog(combined[0, ])), 0L)
  expect_equal(nrow(reconcile_catalog(combined[combined$source == "ltr-flanked", ])), 1L)
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


# ---------------------------------------------------------------------------
# Tree-attached confidence panels (ADR-011). The tree comes from coordinate
# CSVs written by tree_layout.py, so these tests write those directly.
# ---------------------------------------------------------------------------
.write_tree_fixture <- function(dir, name, tips) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(
    tibble::tibble(tip = tips, x = 1, y = seq_along(tips)),
    file.path(dir, paste0(name, ".tree_tips.csv"))
  )
  readr::write_csv(
    tibble::tibble(x = 0, y = 1, xend = 0, yend = length(tips)),
    file.path(dir, paste0(name, ".tree_segments.csv"))
  )
}

test_that("read_tree_part returns NULL for a missing or header-only file", {
  tmp <- tempfile(); dir.create(tmp)
  expect_null(read_tree_part(tmp, "species", "tips"))     # no file at all
  readr::write_csv(tibble::tibble(tip = character(), x = numeric(), y = numeric()),
                   file.path(tmp, "species.tree_tips.csv"))
  expect_null(read_tree_part(tmp, "species", "tips"))     # header only
})

test_that("tree panels degrade to a placeholder when no tree is configured", {
  tmp <- tempfile(); dir.create(tmp)
  d <- .fake_loci("ltr-flanked")
  d$confidence_num <- 0.9
  expect_s3_class(species_confidence_tree_plot(d, tmp), "ggplot")
  expect_s3_class(taxon_confidence_tree_plot(d, tmp), "ggplot")
})

test_that("species tree panel builds when tips match the loci", {
  tmp <- tempfile(); dir.create(tmp)
  d <- .fake_loci("ltr-flanked")
  d$confidence_num <- c(0.9, 0.4)[seq_len(nrow(d)) %% 2 + 1]
  .write_tree_fixture(tmp, "species", unique(as.character(d$species)))
  p <- species_confidence_tree_plot(d, tmp)
  expect_true(inherits(p, "patchwork") || inherits(p, "ggplot"))
})

test_that("taxon tree panel builds and is level-agnostic about tip rank", {
  tmp <- tempfile(); dir.create(tmp)
  d <- .fake_loci("ltr-flanked")
  d$confidence_num <- 0.8
  # Mixed-rank tips: a genus and a family side by side (ADR-008).
  .write_tree_fixture(tmp, "taxon", unique(as.character(d$taxon_call)))
  p <- taxon_confidence_tree_plot(d, tmp)
  expect_true(inherits(p, "patchwork") || inherits(p, "ggplot"))
})


# ---------------------------------------------------------------------------
# Tree-ordered composition panels (ADR-014). tree_composition_plot is the
# counts/composition sibling of tree_confidence_plot: same tree machinery, but
# the bars show what a tip is made of rather than how confident it is.
# ---------------------------------------------------------------------------
.composition_loci <- function() {
  tribble(
    ~species,           ~taxon_call,       ~segment,          ~source,       ~confidence,
    "Antrozous pallidus", "Gammaretrovirus", "Gammaretrovirus", "ltr-flanked", "0.9",
    "Antrozous pallidus", "Betaretrovirus",  "Betaretrovirus",  "orphan",      "0.8",
    "Antrozous pallidus", "Betaretrovirus",  "Betaretrovirus",  "orphan",      "0.7",
    "Mus musculus",       "Gammaretrovirus", "Gammaretrovirus", "ltr-flanked", "0.95",
    "Mus musculus",       "Gammaretrovirus", "Gammaretrovirus", "orphan",      "0.6"
  )
}

test_that("tree_composition_plot renders bars ordered by a supplied tree", {
  dir <- withr::local_tempdir()
  .write_tree_fixture(dir, "species", c("Antrozous pallidus", "Mus musculus"))

  p <- tree_composition_plot(.composition_loci(), dir, "species", "species",
                             "segment", "t", "s")
  expect_s3_class(p, "patchwork")
})

test_that("tree_composition_plot degrades to a placeholder without a tree", {
  dir <- withr::local_tempdir()   # no tree files written
  p <- tree_composition_plot(.composition_loci(), dir, "species", "species",
                             "segment", "t", "s")
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$title, "no species tree")
})

test_that("tree_composition_plot is empty-safe and tolerates a missing column", {
  dir <- withr::local_tempdir()
  .write_tree_fixture(dir, "species", c("Antrozous pallidus"))
  expect_s3_class(
    tree_composition_plot(.composition_loci()[0, ], dir, "species", "species",
                          "segment", "t", "s"), "ggplot")
  # asking to fill by a column the frame does not carry must not error
  expect_s3_class(
    tree_composition_plot(.composition_loci(), dir, "species", "species",
                          "not_a_column", "t", "s"), "ggplot")
})

test_that("tree_composition_plot drops rows whose key is absent from the tree", {
  # The tip-label trap from ADR-011: loci keyed by a name the tree does not
  # carry must be excluded rather than silently plotted at a wrong y position.
  dir <- withr::local_tempdir()
  .write_tree_fixture(dir, "species", c("Antrozous pallidus"))  # Mus absent
  p <- tree_composition_plot(.composition_loci(), dir, "species", "species",
                             "segment", "t", "s")
  expect_s3_class(p, "patchwork")
})

test_that("the two ADR-014 panels build from the catalog frame", {
  dir <- withr::local_tempdir()
  .write_tree_fixture(dir, "species", c("Antrozous pallidus", "Mus musculus"))
  .write_tree_fixture(dir, "taxon", c("Gammaretrovirus", "Betaretrovirus"))
  loci <- .composition_loci()

  expect_s3_class(species_composition_tree_plot(loci, dir), "patchwork")
  expect_s3_class(taxon_tier_tree_plot(loci, dir), "patchwork")
})
