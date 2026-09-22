# testthat coverage for workflow/scripts/plot2sort/style.R
#
# The style module is the single source of every colour and word a figure shows,
# so its promises are tested here rather than trusted: the palette stays
# colour-blind safe, a concept keeps one colour, a genus keeps its colour no
# matter which plot asks, and no file name or raw identifier reaches a reader.
#
# Run via: make test-r

suppressMessages({
  library(testthat)
  library(ggplot2)
})

source(file.path("..", "..", "scripts", "plot2sort", "style.R"))

# Smallest CIELAB distance between any two colours after a colour-vision
# simulation. Below about 10 two colours are hard to tell apart.
.min_distance <- function(cols) {
  lab <- methods::as(colorspace::hex2RGB(cols), "LAB")@coords
  d <- as.matrix(stats::dist(lab))
  diag(d) <- Inf
  min(d)
}


test_that("the categorical palette stays distinguishable under colour blindness", {
  skip_if_not_installed("colorspace")
  pal <- unname(.TOL_MUTED)
  expect_gt(.min_distance(colorspace::deutan(pal)), 15)
  expect_gt(.min_distance(colorspace::protan(pal)), 15)
  expect_gt(.min_distance(colorspace::tritan(pal)), 15)
})


test_that("tiers and fates share colours, one per situation", {
  expect_equal(unname(.FATE_COLOUR[["intact_flank"]]), unname(.TIER_COLOUR[["ltr-flanked"]]))
  expect_equal(unname(.FATE_COLOUR[["mono_ltr_at_orphan"]]), unname(.TIER_COLOUR[["orphan"]]))
  expect_equal(unname(.FATE_COLOUR[["solo"]]), unname(.TIER_COLOUR[["solo-ltr"]]))
  expect_length(unique(.TIER_COLOUR), 3)
})


test_that("genus colours never reuse a tier colour", {
  # One colour, one meaning: indigo is LTR-flanked, so it cannot also be a genus.
  expect_length(intersect(unname(.GENUS_COLOUR), unname(.TIER_COLOUR)), 0)
})


test_that("a genus has the same colour whichever plot asks", {
  alone <- taxon_colours("Betaretrovirus")
  crowded <- taxon_colours(c("Gammaretrovirus", "Lentivirus", "Betaretrovirus", "Other"))
  expect_equal(alone[["Betaretrovirus"]], crowded[["Betaretrovirus"]])
})


test_that("spumavirus genera share a shade family and never take a genus colour", {
  cols <- taxon_colours(c("Simiispumavirus", "Equispumavirus", "Felispumavirus"))
  expect_length(unique(cols), 3)
  expect_length(intersect(unname(cols), unname(.GENUS_COLOUR)), 0)
})


test_that("unresolved and higher-rank calls are grey", {
  cols <- taxon_colours(c("unassigned_at_genus", "Orthoretrovirinae", "Other"))
  expect_true(all(cols %in% c(.GREY_OTHER, .GREY_MID)))
})


test_that("ordinal structure classes run dark to light", {
  lum <- function(hex) sum(grDevices::col2rgb(hex) * c(0.299, 0.587, 0.114))
  expect_lt(lum(.STRUCTURE_COLOUR[["full"]]), lum(.STRUCTURE_COLOUR[["partial"]]))
  expect_lt(lum(.STRUCTURE_COLOUR[["partial"]]), lum(.STRUCTURE_COLOUR[["gene"]]))
})


test_that("species show their readable config name", {
  map <- list(Desmodus_rotundus = "Desmodus rotundus")
  expect_equal(display_species("Desmodus_rotundus", map), "Desmodus rotundus")
})


test_that("a species missing from the config never shows its file name", {
  expect_equal(display_species("Myotis_myotis", list()), "Myotis myotis")
  expect_equal(display_species("Myotis_myotis", NULL), "Myotis myotis")
})


test_that("identifiers become plain words, listed or not", {
  expect_equal(display_label("mono_ltr_at_orphan"), "MonoLTR at an orphan")
  expect_equal(display_label("unassigned_at_genus"), "Unassigned at genus")
  expect_equal(display_label("some_new_value"), "Some new value")
})


test_that("taxa are italic and non-taxa are not", {
  labels <- italic_labels(c("Betaretrovirus", "Other", "Unassigned at genus"))
  expect_match(deparse(labels[[1]]), "italic")
  expect_false(grepl("italic", deparse(labels[[2]])))
  expect_false(grepl("italic", deparse(labels[[3]])))
})


test_that("no label or word in the module uses a dash as punctuation", {
  words <- c(.LABELS, names(.LABELS))
  expect_false(any(grepl(" - |->", words)))
})


test_that("a stage PDF embeds the house font", {
  skip_if_not(nzchar(Sys.which("pdffonts")), "pdffonts not installed")
  path <- withr::local_tempfile(fileext = ".pdf")
  p <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) + geom_point() + theme_retroseek()
  save_stage_pdf(list(key_page("Stage", "What it shows.", .TIER_COLOUR, "A page"), p), path)
  fonts <- system2("pdffonts", path, stdout = TRUE)
  expect_true(any(grepl("IBMPlexSans", fonts)))
})
