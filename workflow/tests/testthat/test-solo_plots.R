# testthat coverage for workflow/scripts/solo_ltr/solo_plots.R
#
# The figure builders only, not the CLI. What can go wrong quietly: the
# cross-genome pages dropping out of the canonical species order (so a species
# sits on a different row here than on every other page), and a fate drawn in
# a colour other than its house colour.
#
# Run via: make test-r

suppressMessages({
  library(testthat)
  library(ggplot2)
  library(data.table)
})

.dir <- file.path("..", "..", "scripts")
source(file.path(.dir, "plot2sort", "style.R"))
source(file.path(.dir, "plot2sort", "helpers.R"))
source(file.path(.dir, "plot2sort", "tree_axis.R"))
source(file.path(.dir, "solo_ltr", "solo_plots.R"))

# Two genomes with data, listed in config order "Mus, Homo, Desmodus". Desmodus
# has no data and must still keep its row.
.order <- c("Mus musculus", "Homo sapiens", "Desmodus rotundus")
.report <- data.table(species = c("Homo sapiens", "Mus musculus"),
                      solo_to_intact_ratio = c(12.5, 30))

test_that("the ratio page puts every configured species on rows, in config order", {
  p <- solo_intact_ratio_plot(.report, tree = NULL, order = .order)
  # Bottom row first, so the config's first species is on top.
  expect_equal(p$scales$get_scales("x")$limits, rev(.order))
  expect_s3_class(p$coordinates, "CoordFlip")
})

test_that("the composition page colours each fate with its house colour", {
  candidates <- data.table(
    species = c("Homo sapiens", "Homo sapiens", "Mus musculus"),
    fate = c("solo", "intact_flank", "mono_ltr_at_orphan")
  )
  p <- class_composition_plot(candidates, tree = NULL, order = .order)
  built <- ggplot_build(p)$data[[1]]
  expect_setequal(unique(built$fill), unname(.FATE_COLOUR))
})

test_that("with a host tree the ratio page is composed beside it", {
  tree <- list(
    tips = data.frame(tip = rev(.order), x = 1, y = 1:3),
    segments = data.frame(x = 0, y = 1, xend = 0, yend = 3)
  )
  expect_s3_class(solo_intact_ratio_plot(.report, tree = tree, order = .order),
                  "patchwork")
})

test_that("funnel labels are words, never dash separators", {
  funnel <- data.table(
    stage = c("raw_hits", "accepted_hits", "merged_candidates",
              "intact_flank", "mono_ltr_at_orphan", "solo"),
    count = c(1000, 800, 500, 100, 50, 350)
  )
  p <- funnel_plot(funnel, "Homo sapiens")
  labels <- levels(p$data$label)
  expect_false(any(grepl("^- | - |->", labels)))
  expect_match(p$labels$subtitle, "^Homo sapiens\\. ")
})

test_that("empty inputs give a placeholder page rather than an error", {
  expect_s3_class(solo_intact_ratio_plot(.report[0]), "ggplot")
  expect_s3_class(class_composition_plot(data.table(species = character(),
                                                    fate = character())), "ggplot")
})
