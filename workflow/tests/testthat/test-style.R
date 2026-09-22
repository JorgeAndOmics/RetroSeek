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


test_that("taxon legends run most abundant first, leftovers last", {
  lv <- taxon_levels(c("Betaretrovirus", "Unassigned at genus", "Gammaretrovirus", "Other"),
                     c(10, 99, 20, 50))
  expect_equal(lv, c("Gammaretrovirus", "Betaretrovirus", "Unassigned at genus", "Other"))
})

test_that("the blank theme really removes grid lines and axes", {
  p <- ggplot2::ggplot(data.frame(x = 1:2, y = 1:2), ggplot2::aes(x, y)) +
    ggplot2::geom_point() + theme_retroseek_blank()
  th <- ggplot2::calc_element("panel.grid.major.x", ggplot2::theme_get() + p$theme)
  expect_s3_class(th, "element_blank")
  expect_s3_class(ggplot2::calc_element("axis.text.x", p$theme), "element_blank")
  # Titles keep the house style.
  expect_equal(p$theme$plot.title$face, "bold")
})

test_that("ERV classes wear the colour of their defining genus", {
  # Class I is gamma-like, Class II beta-like, Class III spumaviral (Jern/Blomberg).
  expect_equal(.ERV_CLASS_COLOUR[["Class I"]], .GENUS_COLOUR[["Gammaretrovirus"]])
  expect_equal(.ERV_CLASS_COLOUR[["Class II"]], .GENUS_COLOUR[["Betaretrovirus"]])
  expect_true(.ERV_CLASS_COLOUR[["Class III"]] %in% .SPUMA_SHADES)
})

test_that("ordinal maps get darker as evidence gets stronger, and 'none' is grey", {
  lum <- function(h) sum(grDevices::col2rgb(h) * c(0.299, 0.587, 0.114))
  expect_lt(lum(.CONFIDENCE_COLOUR[["HC"]]), lum(.CONFIDENCE_COLOUR[["LC"]]))
  expect_lt(lum(.RANK_COLOUR[["genus"]]), lum(.RANK_COLOUR[["subfamily"]]))
  expect_lt(lum(.RANK_COLOUR[["subfamily"]]), lum(.RANK_COLOUR[["family"]]))
  expect_equal(.RANK_COLOUR[["none"]], .GREY_OTHER)
  expect_lt(lum(.DOMAIN_TIER_COLOUR[["domain_selected"]]),
            lum(.DOMAIN_TIER_COLOUR[["domain_unlisted"]]))
  expect_equal(.DOMAIN_TIER_COLOUR[["non_domain"]], .GREY_OTHER)
})

test_that("a small set of plain categories never borrows a tier colour", {
  cols <- category_colours(c("a", "b", "c", "d", "e", "f"))
  expect_length(intersect(cols, .TIER_COLOUR), 0L)
  expect_equal(unname(category_colours(c("a", "Other (3)"))[2]), .GREY_OTHER)
})

test_that("taxon legends italicise taxa and keep leftovers upright", {
  labs <- taxon_labels(c("Betaretrovirus", "unassigned_at_genus", "Other (4)"))
  expect_identical(labs[[1]], quote(italic("Betaretrovirus")))
  expect_identical(labs[[2]], "Unassigned at genus")
  expect_identical(labs[[3]], "Other (4)")
})

test_that("the taxon fill scale keeps fixed colours and lists taxa by abundance", {
  sc <- scale_fill_taxon(c("Gammaretrovirus", "Betaretrovirus", "Retroviridae"),
                         weights = c(5, 9, 100))
  expect_equal(sc$palette(3)[["Gammaretrovirus"]], .GENUS_COLOUR[["Gammaretrovirus"]])
  # Retroviridae is not resolved to a genus, but it is still a taxon: by abundance.
  expect_equal(sc$breaks, c("Retroviridae", "Betaretrovirus", "Gammaretrovirus"))
})

test_that("heatmap text is light on the dark end of the ramp and dark elsewhere", {
  expect_equal(ink_on_ramp(c(1, 50, 100)), c(.INK, .INK, .PAPER))
  expect_equal(ink_on_ramp(5), .INK)   # one value: no ramp, default ink
})

test_that("the LCA method has a readable label", {
  expect_equal(display_label("lca"), "Weighted LCA")
})

test_that("page_titles reads plain and tree-composed pages alike", {
  plain <- ggplot2::ggplot() + ggplot2::labs(title = "Plain page")
  composed <- patchwork::wrap_plots(ggplot2::ggplot(), ggplot2::ggplot()) +
    patchwork::plot_annotation(title = "Composed page")
  untitled <- ggplot2::ggplot()
  expect_equal(page_titles(list(plain, composed, untitled)),
               c("Plain page", "Composed page", ""))
})

test_that("taxon_factor orders lineages like the legend, so stacks match it", {
  f <- taxon_factor(c("Gammaretrovirus", "Betaretrovirus", "unassigned_at_genus"),
                    weights = c(5, 9, 100))
  expect_equal(levels(f), c("Betaretrovirus", "Gammaretrovirus", "unassigned_at_genus"))
})

test_that("the ramp fill scale runs from the lightest to the darkest ramp step", {
  sc <- scale_fill_ramp(name = "Loci")
  expect_equal(toupper(sc$palette(0)), toupper(seq_colours(9)[1]))
  expect_equal(toupper(sc$palette(1)), toupper(seq_colours(9)[9]))
})

test_that("a spumavirus genus keeps its shade whichever others share the plot", {
  alone <- taxon_colours("Felispumavirus")[["Felispumavirus"]]
  crowded <- taxon_colours(c("Bovispumavirus", "Equispumavirus", "Felispumavirus"))
  expect_equal(crowded[["Felispumavirus"]], alone)
})

test_that("the probeset's foamy-virus label takes the spumavirus shades", {
  expect_true(taxon_colours("Spumaretrovirus")[["Spumaretrovirus"]] %in% .SPUMA_SHADES)
  # The subfamily itself stays a higher rank.
  expect_equal(taxon_colours("Spumaretrovirinae")[["Spumaretrovirinae"]], .GREY_MID)
})

test_that("viruses wear shades of their genus colour, the first the genus colour itself", {
  cols <- virus_colours(c("Mouse mammary tumor virus", "Jaagsiekte sheep retrovirus",
                          "Feline leukemia virus", "Other (4)"),
                        c("Betaretrovirus", "Betaretrovirus", "Gammaretrovirus", "Other"))
  expect_equal(cols[["Jaagsiekte sheep retrovirus"]], .GENUS_COLOUR[["Betaretrovirus"]])
  expect_false(cols[["Mouse mammary tumor virus"]] == .GENUS_COLOUR[["Betaretrovirus"]])
  expect_equal(cols[["Feline leukemia virus"]], .GENUS_COLOUR[["Gammaretrovirus"]])
  expect_equal(cols[["Other (4)"]], .GREY_OTHER)
  expect_equal(anyDuplicated(unname(cols)), 0L)
})

test_that("probe colours depend only on the probe set, not the order of the data", {
  expect_equal(probe_colours(c("POL", "GAG", "ENV", "POL")),
               probe_colours(c("ENV", "GAG", "POL")))
  expect_length(probe_colours(c("POL", "GAG")), 2L)
})

test_that("a virus with no lineage label is grey rather than an error", {
  cols <- virus_colours(c("HIV", "FFV"), c(NA, "Spumaretrovirus"))
  expect_equal(cols[["HIV"]], .GREY_OTHER)
  expect_true(cols[["FFV"]] %in% .SPUMA_COLOUR)
})

test_that("past nine levels the palette keeps its colours and adds lighter tints", {
  cols <- category_colours(sprintf("p%02d", 1:12))
  expect_equal(unname(cols[1:9]), unname(.TOL_MUTED[.CATEGORY_ORDER]))
  expect_equal(anyDuplicated(unname(cols)), 0L)
  expect_length(intersect(cols[10:12], .TOL_MUTED), 0L)
})

test_that("within a genus the most abundant virus takes the genus colour", {
  cols <- virus_colours(c("Mouse mammary tumor virus", "Human endogenous retrovirus K"),
                        c("Betaretrovirus", "Betaretrovirus"), weights = c(500, 20))
  expect_equal(cols[["Mouse mammary tumor virus"]], .GENUS_COLOUR[["Betaretrovirus"]])
})
