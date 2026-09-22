# testthat coverage for workflow/scripts/plot2sort/tree_axis.R
#
# The tree axis is a coordinate bridge, not a tree library: tree_layout.py writes
# tip and segment CSVs and this module draws them beside a data panel. The things
# that can silently go wrong are alignment (rows not matching tips), ordering (the
# panel in alphabetical order while the tree is in ladderised order), and silent
# dropping (a species absent from the tree vanishing with no explanation). One
# test each.
#
# Run via: make test-r

suppressMessages({
  library(testthat)
  library(ggplot2)
  library(readr)
})

.script_dir <- file.path("..", "..", "scripts", "plot2sort")
source(file.path(.script_dir, "helpers.R"))
source(file.path(.script_dir, "tree_axis.R"))

# A three-tip tree, deliberately NOT in alphabetical order, so a test that passes
# by accident under alphabetical sorting is impossible.
.tips <- data.frame(
  tip = c("Mus musculus", "Homo sapiens", "Desmodus rotundus"),
  x = c(1, 1, 1),
  y = c(1, 2, 3),
  stringsAsFactors = FALSE
)
.segs <- data.frame(x = 0, y = 1, xend = 0, yend = 3)

.panel <- function(species = .tips$tip) {
  d <- data.frame(species = species, n = seq_along(species))
  ggplot(d, aes(x = .data$species, y = .data$n)) + geom_col()
}


test_that("read_tree_part returns NULL when no tree is configured", {
  dir <- withr::local_tempdir()
  expect_null(read_tree_part(dir, "species", "tips"))
})


test_that("read_tree_part returns NULL for a header-only file", {
  # tree_layout.py always writes the CSVs, empty when there is no tree, so the
  # DAG stays stable. An empty file must read as "no tree", not as zero tips.
  dir <- withr::local_tempdir()
  write_csv(data.frame(tip = character(), x = numeric(), y = numeric()),
            file.path(dir, "species.tree_tips.csv"))
  expect_null(read_tree_part(dir, "species", "tips"))
})


test_that("attach_tree_axis returns NULL when there is no tree", {
  expect_null(attach_tree_axis(.panel(), NULL, NULL, .tips$tip, "t", "s"))
})


test_that("attach_tree_axis returns NULL when nothing overlaps the tree", {
  # Better to emit no variant than a figure with every row dropped.
  expect_null(attach_tree_axis(.panel(), .tips, .segs,
                               c("Gallus gallus"), "t", "s"))
})


test_that("the variant orders rows by tip order, not alphabetically", {
  p <- attach_tree_axis(.panel(), .tips, .segs, .tips$tip, "t", "s")
  expect_false(is.null(p))
  # The data panel is the second element of the patchwork; its discrete scale
  # carries the tree's tip order as its limits.
  limits <- p[[2]]$scales$get_scales("x")$limits
  expect_equal(limits, c("Mus musculus", "Homo sapiens", "Desmodus rotundus"))
})


test_that("a species missing from the tree is named in the subtitle", {
  # tree_layout.py warns in the log, but a reader sees only the figure, so the
  # omission has to be visible there too.
  p <- attach_tree_axis(.panel(c(.tips$tip, "Gallus gallus")), .tips, .segs,
                        c(.tips$tip, "Gallus gallus"), "t", "base subtitle")
  expect_match(p$patches$annotation$subtitle, "Gallus gallus")
  expect_match(p$patches$annotation$subtitle, "base subtitle")
})


test_that("the tree column and the data panel share a row range", {
  # Alignment is the whole point: a mismatch silently shifts every bar by a row.
  # A discrete scale expands by 0.6 either side, which is exactly the range
  # tree_column() is given, so the two line up by construction.
  p <- attach_tree_axis(.panel(), .tips, .segs, .tips$tip, "t", "s")
  tree_y <- p[[1]]$scales$get_scales("y")$limits
  expect_equal(tree_y, c(0.4, nrow(.tips) + 0.6))
})


test_that("the data panel does not print its own row labels", {
  # The tree supplies the tip labels; printing them twice is the failure mode
  # that scale_categorical_axis() would otherwise reintroduce.
  p <- attach_tree_axis(.panel(), .tips, .segs, .tips$tip, "t", "s")
  expect_true(inherits(p[[2]]$theme$axis.text.y, "element_blank"))
})
