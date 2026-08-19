# testthat coverage for taxonomy_segments.R's catalog loading.
#
# The per-segment panel REUSES builders from taxonomy_plot_generator.R, so it
# must hand them the frame that script's main() prepares. Two properties matter
# and neither is free:
#
#   1. every column character - the catalog writes booleans as the strings
#      "True"/"False", and readr's default inference turns `resolved` into a
#      logical, after which the builders' `filter(resolved == "True")` matches
#      nothing (TRUE coerces to "TRUE", which is not "True");
#   2. the numeric companions (confidence_num, n_hits, completeness_num), which
#      main() derives and this path never ran.
#
# Both failures are SILENT: the builders return their empty placeholder rather
# than erroring, so every segment gets the SAME picture and nothing complains.
# Measured 2026-08-18 on the model 5: all 12 genera had byte-identical
# taxon_composition.png (md5 8c0cd59b) and confidence_gradient.png (63823bbd),
# while structure_class_composition.png - which needs neither property - was
# correctly distinct in all 12.

suppressMessages({
  library(testthat)
  library(tibble)
  library(dplyr)
  library(readr)
})

.script_dir <- file.path("..", "..", "scripts")
source(file.path(.script_dir, "taxonomy", "taxonomy_plot_generator.R"))
# Sourced before taxonomy_segments.R: load_catalog() applies both companion
# sets, and the registry tests below need structure_panel_registry().
source(file.path(.script_dir, "taxonomy", "erv_like_plot_generator.R"))
source(file.path(.script_dir, "taxonomy", "taxonomy_segments.R"))


# A catalog row set written the way taxonomy_classify_loci.py writes it: every
# value a string, booleans as "True"/"False". Carries ALL 26 catalog columns,
# because the per-segment panel now drives the whole builder set and a builder
# missing its column silently degrades to empty_plot() rather than erroring -
# the exact failure this file already guards against.
.write_catalog <- function(path) {
  tribble(
    ~species,  ~source,       ~seqname, ~start,  ~end,    ~strand, ~taxon_call,       ~rank,   ~segment,          ~segment_rank, ~resolved, ~confidence, ~confidence_tag, ~erv_class, ~structure_class, ~domain_tier,      ~oversized, ~canonical_order, ~completeness, ~n_main_genes, ~genes_present, ~is_mosaic, ~mosaic_composition,               ~n_blastx_hits, ~method,     ~id,
    "Sp one",  "ltr-flanked", "chr1",   "1000",  "9000",  "+",     "Gammaretrovirus", "genus", "Gammaretrovirus", "genus",       "True",    "0.980",     "HC",            "I",        "full",           "domain_selected", "False",    "True",           "1.000",       "3",           "GAG,POL,ENV",  "False",    "",                                "12",           "placement", "L0",
    "Sp one",  "orphan",      "chr1",   "20000", "20800", "-",     "Gammaretrovirus", "genus", "Gammaretrovirus", "genus",       "True",    "0.400",     "LC",            "I",        "gene",           "non_domain",      "False",    "False",          "0.333",       "1",           "POL",          "False",    "",                                "3",            "lca",       "L1",
    "Sp two",  "ltr-flanked", "chr2",   "5000",  "14000", "+",     "Gammaretrovirus", "genus", "Gammaretrovirus", "genus",       "True",    "0.910",     "HC",            "I",        "partial",        "domain_unlisted", "False",    "True",           "0.667",       "2",           "GAG,POL",      "True",     "GAG:Betaretrovirus;POL:Gammaretrovirus", "9",      "placement", "L2",
    "Sp two",  "ltr-flanked", "chr2",   "40000", "48000", "-",     "Betaretrovirus",  "genus", "Betaretrovirus",  "genus",       "True",    "0.750",     "HC",            "II",       "partial",        "domain_selected", "False",    "False",          "0.667",       "2",           "POL,ENV",      "False",    "",                                "7",            "lca",       "L3",
    "Sp one",  "orphan",      "chr3",   "100",   "900",   "+",     "Betaretrovirus",  "genus", "Betaretrovirus",  "genus",       "True",    "0.550",     "HC",            "II",       "gene",           "non_domain",      "False",    "False",          "0.333",       "1",           "ENV",          "False",    "",                                "0",            "lca",       "L4"
  ) %>% write_csv(path)
  path
}


test_that("resolved stays character, so the builders' 'True' filter still matches", {
  # The regression. With readr's default inference this column comes back
  # logical and every downstream `resolved == "True"` silently matches 0 rows.
  catalog <- load_catalog(.write_catalog(tempfile(fileext = ".csv")))

  expect_type(catalog$resolved, "character")
  expect_equal(nrow(dplyr::filter(catalog, resolved == "True")), 5L)
})


test_that("the numeric companions the reused builders expect are present", {
  catalog <- load_catalog(.write_catalog(tempfile(fileext = ".csv")))

  expect_true(all(c("confidence_num", "n_hits", "completeness_num") %in% names(catalog)))
  expect_equal(catalog$confidence_num, c(0.98, 0.4, 0.91, 0.75, 0.55))
  expect_equal(catalog$n_hits, c(12L, 3L, 9L, 7L, 0L))
})


test_that("taxon_composition_plot on a loaded catalog is not the empty placeholder", {
  # Guards the actual symptom rather than the mechanism: if this regressed,
  # every segment's PNG would be the same "no confident taxon calls" image.
  catalog <- load_catalog(.write_catalog(tempfile(fileext = ".csv")))
  p <- taxon_composition_plot(catalog)

  expect_equal(p$labels$title, "ERV taxon composition")
  expect_false(identical(p$labels$title, "no confident taxon calls"))
})


test_that("confidence_gradient_plot on a loaded catalog is not the empty placeholder", {
  catalog <- load_catalog(.write_catalog(tempfile(fileext = ".csv")))
  p <- confidence_gradient_plot(catalog)

  expect_equal(p$labels$title, "Confidence distribution (gradient)")
})


test_that("filtering to one segment still yields a real plot", {
  # What the per-segment loop actually does. A segment with rows must not
  # produce the same picture as any other segment.
  catalog <- load_catalog(.write_catalog(tempfile(fileext = ".csv")))
  gamma <- dplyr::filter(catalog, segment == "Gammaretrovirus")
  beta <- dplyr::filter(catalog, segment == "Betaretrovirus")

  expect_equal(nrow(gamma), 3L)
  expect_equal(nrow(beta), 2L)
  # Distinct underlying data is the thing that made the PNGs distinct.
  expect_false(identical(
    taxon_composition_plot(gamma)$data, taxon_composition_plot(beta)$data
  ))
})


test_that("add_numeric_companions tolerates an empty frame", {
  # A segment can legitimately be empty; the panel must stay well-formed.
  empty <- tibble(confidence = character(), n_blastx_hits = character(),
                  completeness = character())
  expect_equal(nrow(add_numeric_companions(empty)), 0L)
})


# ---------------------------------------------------------------------------
# Panel registry
#
# The per-segment subset used to be a hand-written list of 3 emit() calls that
# had to be kept in step with a second hand-written list in
# taxonomy_plot_generator.R's main(). It drifted, which is how the panel ended
# up at 3 plots while the full panel grew to 24. One registry, two consumers,
# so a new plot is offered to both by construction.
#
# `segment = FALSE` marks entries that are DEGENERATE for a single segment:
# erv_class is constant within a genus (measured 2026-08-18: exactly 1.00
# distinct values per segment), and the taxonomy cladogram collapses to one tip.
# ---------------------------------------------------------------------------
.full_registry <- function() c(panel_registry(), structure_panel_registry())


test_that("every registry entry is well formed and uniquely named", {
  reg <- .full_registry()
  files <- vapply(reg, function(e) e$file, character(1))

  expect_true(length(reg) > 0L)
  expect_equal(anyDuplicated(files), 0L)
  expect_true(all(grepl("\\.png$", files)))
  expect_true(all(vapply(reg, function(e) is.function(e$build), logical(1))))
  expect_true(all(vapply(reg, function(e) is.logical(e$segment), logical(1))))
})


test_that("the registry covers the whole published panel", {
  # 24 taxonomy + 7 structure. If a builder is added without a registry entry
  # it silently stops being emitted, which is the failure this pins.
  expect_equal(length(panel_registry()), 24L)
  expect_equal(length(structure_panel_registry()), 7L)
})


test_that("exactly the three degenerate plots are excluded from segments", {
  reg <- .full_registry()
  excluded <- vapply(Filter(function(e) !e$segment, reg),
                     function(e) e$file, character(1))

  expect_setequal(excluded, c("erv_class_composition.png",
                              "taxon_confidence_tree.png",
                              "taxon_tier_tree.png"))
  expect_equal(length(segment_panel(reg, "full")), 28L)
})


test_that("the species-based trees are KEPT for segments", {
  # They answer a real per-genus question - this lineage's burden across the
  # host phylogeny - unlike the taxonomy cladogram, which is one tip here.
  files <- vapply(segment_panel(.full_registry(), "full"),
                  function(e) e$file, character(1))

  expect_true("species_composition_tree.png" %in% files)
  expect_true("species_confidence_tree.png" %in% files)
})


test_that("segment_panel honours the config knob", {
  reg <- .full_registry()

  expect_equal(length(segment_panel(reg, "none")), 0L)
  curated <- vapply(segment_panel(reg, "curated"), function(e) e$file, character(1))
  expect_setequal(curated, c("taxon_composition.png", "confidence_gradient.png",
                             "structure_class_composition.png"))
})


test_that("every segment builder runs on a real catalog slice without erroring", {
  # The panel is driven from catalog.csv, so each builder must find its columns
  # there. A missing column degrades to empty_plot() rather than raising, so
  # this asserts a ggplot comes back AND that the non-tree entries are not the
  # placeholder (tree entries legitimately are, with no tree_dir supplied).
  catalog <- load_catalog(.write_catalog(tempfile(fileext = ".csv")))
  gamma <- dplyr::filter(catalog, segment == "Gammaretrovirus")
  ctx <- list(tree_dir = "", confidence_min = 0.5)

  for (e in segment_panel(.full_registry(), "full")) {
    d <- if (identical(e$data, "loci")) {
      dplyr::filter(gamma, source == "ltr-flanked")
    } else {
      gamma
    }
    p <- e$build(d, ctx)
    expect_s3_class(p, "gg")
  }
})


test_that("the loci/combined split is preserved", {
  # Composition and mosaic panels are LTR-flanked only (the taxon-founded
  # assembly); confidence and tier panels span both tiers. Losing that would
  # quietly mix orphans into plots that are meant to exclude them.
  reg <- panel_registry()
  by_file <- setNames(reg, vapply(reg, function(e) e$file, character(1)))

  expect_equal(by_file[["taxon_composition.png"]]$data, "loci")
  expect_equal(by_file[["mosaic_burden.png"]]$data, "loci")
  expect_equal(by_file[["domain_tier_composition.png"]]$data, "loci")
  expect_equal(by_file[["confidence.png"]]$data, "combined")
  expect_equal(by_file[["source_yield.png"]]$data, "combined")
  expect_equal(by_file[["structure_class_composition.png"]]$data, "combined")
})
