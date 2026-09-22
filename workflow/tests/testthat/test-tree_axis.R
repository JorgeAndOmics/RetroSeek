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


test_that("an unset tree directory means no tree, not an error", {
  expect_null(read_tree_part(NULL, "species", "tips"))
  expect_null(read_tree_part("", "species", "tips"))
  expect_null(read_species_tree(NULL))
})


test_that("a panel already on rows keeps its own x axis text", {
  p <- ggplot(data.frame(g = c("a", "b"), species = c("Homo sapiens", "Mus musculus")),
              aes(x = .data$g, y = .data$species)) + geom_tile() +
    theme(axis.text.x = element_text(angle = 30))
  out <- species_rows(p, c("Homo sapiens", "Mus musculus"), axis = "y",
                      fallback_order = c("Homo sapiens", "Mus musculus"))
  expect_equal(out$theme$axis.text.x$angle, 30)
})

test_that("on_rows takes the tree and order from a panel ctx, and works without one", {
  ctx <- list(species_tree = NULL, species_order = c("Mus musculus", "Homo sapiens"))
  p <- on_rows(.panel(c("Homo sapiens", "Mus musculus")), c("Homo sapiens", "Mus musculus"), ctx)
  expect_equal(p$scales$get_scales("x")$limits, c("Homo sapiens", "Mus musculus"))
  expect_s3_class(on_rows(.panel("Homo sapiens"), "Homo sapiens"), "ggplot")
  expect_s3_class(on_rows(.panel(c("Homo sapiens", "Mus musculus")),
                          c("Homo sapiens", "Mus musculus"),
                          list(species_tree = .tree)), "patchwork")
})

test_that("panel_ctx carries the readable config order and the host tree", {
  dir <- withr::local_tempdir()
  write.csv(data.frame(tip = c("House mouse", "Homo sapiens"), x = 1, y = 1:2),
            file.path(dir, "species.tree_tips.csv"), row.names = FALSE)
  write.csv(data.frame(x = 0, y = 1, xend = 0, yend = 2),
            file.path(dir, "species.tree_segments.csv"), row.names = FALSE)
  cfg <- list(species = list(Mus_musculus = "House mouse", Homo_sapiens = ""),
              classification = list(confidence_min = 0.7))
  ctx <- panel_ctx(cfg, species_tree_dir = dir)
  expect_equal(ctx$species_order, c("House mouse", "Homo sapiens"))
  expect_equal(ctx$species_tree$tips$tip, c("House mouse", "Homo sapiens"))
  expect_equal(ctx$confidence_min, 0.7)
  expect_null(panel_ctx(list())$species_tree)
})

test_that("render_panel feeds each entry the tier scope it declares", {
  reg <- list(
    list(name = "a", data = "loci", build = function(d, ctx) nrow(d)),
    list(name = "b", data = "combined", build = function(d, ctx) nrow(d)))
  expect_equal(render_panel(reg, data.frame(x = 1), data.frame(x = 1:3), list()),
               list(1L, 3L))
})

test_that("species_facets stacks one row per species, first configured on top", {
  d <- data.frame(species = c("Homo sapiens", "Mus musculus", "Homo sapiens"),
                  value = c(1, 2, 3))
  p <- ggplot(d, aes(x = .data$value)) + geom_histogram(bins = 3)
  ctx <- list(species_order = c("Mus musculus", "Homo sapiens"))
  faceted <- species_facets(p, ctx)
  expect_equal(levels(faceted$data$species), c("Mus musculus", "Homo sapiens"))
  expect_s3_class(faceted$facet, "FacetGrid")
})
