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
source(file.path(.script_dir, "taxonomy", "taxonomy_segments.R"))


# A catalog row set written the way taxonomy_classify_loci.py writes it: every
# value a string, booleans as "True"/"False".
.write_catalog <- function(path) {
  tribble(
    ~species,   ~source,       ~segment,          ~taxon_call,       ~rank,   ~resolved, ~confidence, ~confidence_tag, ~method,     ~is_mosaic, ~structure_class, ~n_blastx_hits, ~completeness,
    "Sp one",   "ltr-flanked", "Gammaretrovirus", "Gammaretrovirus", "genus", "True",    "0.980",     "HC",            "placement", "False",    "full",           "12",           "1.000",
    "Sp one",   "orphan",      "Gammaretrovirus", "Gammaretrovirus", "genus", "True",    "0.400",     "LC",            "lca",       "False",    "gene",           "3",            "0.333",
    "Sp two",   "ltr-flanked", "Betaretrovirus",  "Betaretrovirus",  "genus", "True",    "0.750",     "HC",            "lca",       "True",     "partial",        "7",            "0.667"
  ) %>% write_csv(path)
  path
}


test_that("resolved stays character, so the builders' 'True' filter still matches", {
  # The regression. With readr's default inference this column comes back
  # logical and every downstream `resolved == "True"` silently matches 0 rows.
  catalog <- load_catalog(.write_catalog(tempfile(fileext = ".csv")))

  expect_type(catalog$resolved, "character")
  expect_equal(nrow(dplyr::filter(catalog, resolved == "True")), 3L)
})


test_that("the numeric companions the reused builders expect are present", {
  catalog <- load_catalog(.write_catalog(tempfile(fileext = ".csv")))

  expect_true(all(c("confidence_num", "n_hits", "completeness_num") %in% names(catalog)))
  expect_equal(catalog$confidence_num, c(0.98, 0.4, 0.75))
  expect_equal(catalog$n_hits, c(12L, 3L, 7L))
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

  expect_equal(nrow(gamma), 2L)
  expect_equal(nrow(beta), 1L)
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
