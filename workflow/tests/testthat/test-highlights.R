# testthat coverage for workflow/scripts/highlights.R: the short, slide-shaped
# PDF of headline findings. Each page builder on a small catalog, the focal page
# when lineages are listed and when not, and the whole PDF end to end.

suppressMessages({
  library(testthat)
  library(dplyr)
})

.script_dir <- file.path("..", "..", "scripts")
source(file.path(.script_dir, "highlights.R"))
use_retroseek_style()

# A catalog as highlights.R reads it: readable species names, tier, lineage.
.catalog <- function() {
  tibble::tibble(
    species = c(rep("Homo sapiens", 4), rep("Mus musculus", 3)),
    source  = c("ltr-flanked", "orphan", "orphan", "ltr-flanked",
                "ltr-flanked", "ltr-flanked", "orphan"),
    segment = c("Gammaretrovirus", "Gammaretrovirus", "Betaretrovirus",
                "unassigned_at_genus", "Betaretrovirus", "Betaretrovirus",
                "Alpharetrovirus")
  )
}
.axis <- function() species_columns(c("Homo sapiens", "Mus musculus"))

test_that("every headline page is a titled composition", {
  for (page in list(loci_by_tier_page(.catalog(), .axis()),
                    lineage_mix_page(.catalog(), .axis()),
                    genus_matrix_page(.catalog(), .axis()))) {
    expect_s3_class(page, "patchwork")
    expect_true(nzchar(page_titles(list(page))))
  }
})

test_that("the genus matrix has one bubble per host and lineage, sized by loci", {
  cells <- genus_matrix_cells(.catalog())
  expect_equal(sum(cells$n), nrow(.catalog()))
  expect_equal(cells$n[cells$species == "Mus musculus" &
                         cells$segment == "Betaretrovirus"], 2L)
})

test_that("the focal page draws one panel per listed lineage, absent ones kept", {
  page <- focal_lineages_page(.catalog(), .axis(),
                              c("Alpharetrovirus", "Deltaretrovirus"))
  expect_s3_class(page, "patchwork")
  expect_null(focal_lineages_page(.catalog(), .axis(), character(0)))
})

test_that("the PDF has a key page and four pages, three without focal lineages", {
  expect_length(highlights_pages(.catalog(), .axis(), "Alpharetrovirus"), 4L)
  expect_length(highlights_pages(.catalog(), .axis(), NULL), 3L)
})
