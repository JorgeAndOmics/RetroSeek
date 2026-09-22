# testthat coverage for workflow/scripts/plot2sort/tree_axis.R
#
# Species go on rows in one canonical order everywhere. What can go wrong quietly:
# a plot ordering species its own way, a species with no data losing its row (so
# rows shift between pages), and the tree drawn out of line with the rows.
#
# Run via: make test-r

suppressMessages({
  library(testthat)
  library(ggplot2)
})

.dir <- file.path("..", "..", "scripts", "plot2sort")
source(file.path(.dir, "style.R"))
source(file.path(.dir, "helpers.R"))
source(file.path(.dir, "tree_axis.R"))

# A three-tip tree deliberately NOT in alphabetical order, so a pass by accident
# under alphabetical sorting is impossible. y = 1 is the bottom row.
.tree <- list(
  tips = data.frame(tip = c("Mus musculus", "Homo sapiens", "Desmodus rotundus"),
                    x = c(1, 1, 1), y = c(1, 2, 3), stringsAsFactors = FALSE),
  segments = data.frame(x = 0, y = 1, xend = 0, yend = 3)
)

.panel <- function(species) {
  ggplot(data.frame(species = species, n = seq_along(species)),
         aes(x = .data$species, y = .data$n)) + geom_col() +
    labs(title = "Title", subtitle = "Subtitle")
}


test_that("read_tree_part returns NULL when no tree is configured", {
  expect_null(read_tree_part(withr::local_tempdir(), "species", "tips"))
})


test_that("read_species_tree returns NULL for a header-only tips file", {
  dir <- withr::local_tempdir()
  readr::write_csv(data.frame(tip = character(), x = numeric(), y = numeric()),
                   file.path(dir, "species.tree_tips.csv"))
  expect_null(read_species_tree(dir))
})


test_that("with a tree, rows follow the tree, bottom first", {
  expect_equal(species_order(c("Homo sapiens", "Desmodus rotundus"), .tree),
               c("Mus musculus", "Homo sapiens", "Desmodus rotundus"))
})


test_that("a species with no data keeps its row", {
  # Mus musculus is absent from the data yet still has a row, so rows never
  # shift between pages.
  expect_true("Mus musculus" %in% species_order("Homo sapiens", .tree))
})


test_that("without a tree, the config order runs top to bottom", {
  expect_equal(species_order("b", NULL, fallback_order = c("a", "b", "c")),
               c("c", "b", "a"))
})


test_that("a species outside the canonical list is kept, not dropped", {
  expect_true("Gallus gallus" %in% species_order(c("Gallus gallus"), .tree))
})


test_that("without a tree the panel is flipped with italic species labels", {
  p <- species_rows(.panel(c("Homo sapiens", "Mus musculus")), c("Homo sapiens", "Mus musculus"),
                    fallback_order = c("Homo sapiens", "Mus musculus"))
  expect_true(inherits(p$coordinates, "CoordFlip"))
  expect_equal(p$theme$axis.text.y$face, "italic")
})


test_that("with a tree the rows line up with the tips", {
  p <- species_rows(.panel(.tree$tips$tip), .tree$tips$tip, tree = .tree)
  expect_equal(p[[2]]$scales$get_scales("x")$limits, .tree$tips$tip[order(.tree$tips$y)])
  expect_equal(p[[1]]$scales$get_scales("y")$limits, c(0.4, 3.6))
})


test_that("the title sits above the tree and panel, not inside the panel", {
  p <- species_rows(.panel(.tree$tips$tip), .tree$tips$tip, tree = .tree)
  expect_equal(p$patches$annotation$title, "Title")
  expect_null(p[[2]]$labels$title)
})


test_that("a species the tree cannot place drops the tree and says why", {
  p <- species_rows(.panel(c("Homo sapiens", "Gallus gallus")),
                    c("Homo sapiens", "Gallus gallus"), tree = .tree)
  expect_true(inherits(p, "ggplot"))
  expect_match(p$labels$caption, "Gallus gallus")
})
