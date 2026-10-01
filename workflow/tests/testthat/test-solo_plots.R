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


# ---- solos by the repeat class of their bait family (backlog 16) ----

.ratios <- data.table(
  species = c(rep("Homo sapiens", 4), "Mus musculus"),
  genome = c(rep("Homo_sapiens", 4), "Mus_musculus"),
  ltr_family = c("Hsap_F001", "Hsap_F002", "Hsap_F003", "Hsap_F004", "Mmus_F001"),
  ltr_family_class = c("LTR/ERV1", "LINE/L1", "SINE/Alu", "", "LTR/ERVK"),
  intact_elements = c(10L, 2L, 3L, 1L, 50L),
  solos = c(100L, 900L, 500L, 7L, 200L)
)

test_that("Dfam classes fold into a few repeat groups, blank as no label", {
  # fread reads an all-empty column (Dfam off) as logical NA.
  expect_equal(repeat_group(c(NA, NA)), c("No Dfam label", "No Dfam label"))
  # Only ERV1, ERVK and ERVL are ERV LTRs; MaLR, Gypsy and a bare LTR are not.
  expect_equal(
    repeat_group(c("LTR/ERVK", "LTR/ERV1", "LTR/ERVL", "LTR/ERVL-MaLR", "LTR/Gypsy",
                   "LTR")),
    c("ERV LTR", "ERV LTR", "ERV LTR", "Other LTR", "Other LTR", "Other LTR")
  )
  expect_equal(
    repeat_group(c("LINE/L1", "SINE/Alu", "DNA/hAT", "Unknown", "")),
    c("LINE", "SINE", "Other repeat", "Other repeat", "No Dfam label")
  )
})

test_that("the class breakdown sums families, intact elements and solos per group", {
  out <- solos_by_class(.ratios)
  homo_ltr <- out[species == "Homo sapiens" & repeat_group == "ERV LTR"]
  expect_equal(homo_ltr$families, 1L)
  expect_equal(homo_ltr$solos, 100L)
  expect_equal(homo_ltr$solos_per_intact, 10)
  expect_equal(out[species == "Homo sapiens", sum(solos)], 1507L)
  # Every group a genome has gets a row, the unlabelled family included.
  expect_setequal(out[species == "Homo sapiens", repeat_group],
                  c("ERV LTR", "LINE", "SINE", "No Dfam label"))
})

test_that("the by-class page keeps the species order and never uses fate colours", {
  p <- solos_by_class_plot(solos_by_class(.ratios), tree = NULL, order = .order)
  expect_equal(p$scales$get_scales("x")$limits, rev(.order))
  fills <- unique(ggplot_build(p)$data[[1]]$fill)
  expect_false(any(fills %in% .FATE_COLOUR))
})

test_that("the by-class page says so when no family has a Dfam label", {
  unlabelled <- copy(.ratios)[, ltr_family_class := ""]
  p <- solos_by_class_plot(solos_by_class(unlabelled), tree = NULL, order = .order)
  expect_match(p$labels$subtitle, "No family has a Dfam label")
  # Labelling may be on and simply have matched nothing: no advice to switch it on.
  expect_false(grepl("switch on", p$labels$subtitle))
})

test_that("no ratio tables at all give an empty breakdown, not an error", {
  out <- solos_by_class(data.table())
  expect_equal(nrow(out), 0L)
  expect_s3_class(solos_by_class_plot(out, tree = NULL, order = .order), "ggplot")
})

test_that("a genome without a ratio table is skipped, not fatal", {
  expect_null(read_family_ratios(tempfile(fileext = ".csv"), "Mus_musculus", NULL))
})

test_that("the subtrees words stay generic when no known kind is drawn", {
  words <- .subtrees_titles(character(0))
  expect_false(grepl("no family", words$subtitle))
})


# ---- review fixes (better_plots' read of the M6 PDFs, 2026-10-01) ----

test_that("the ratio page claims nothing about ERVs and shades 10 to 100", {
  p <- solo_intact_ratio_plot(.report, tree = NULL, order = .order)
  expect_false(grepl("ERV", p$labels$title))
  band <- ggplot_build(p)$data[[1]]
  # Log scale: the band is drawn in log10 units.
  expect_equal(c(band$ymin, band$ymax), c(1, 2))
})

test_that("the enrichment page draws Observed in a neutral colour", {
  summary_dt <- data.table(
    metric = c("same_class_sister_observed", "same_class_sister_null_mean",
               "same_class_sister_null_sd"),
    value = c("0.6", "0.3", "0.05")
  )
  built <- ggplot_build(tree_enrichment_plot(summary_dt, "Homo sapiens"))
  fills <- unique(built$data[[1]]$fill)
  expect_false(any(fills %in% .FATE_COLOUR))
})

test_that("the solo tree never gives a family a fate colour", {
  tips <- data.table(family = paste0("F", rep(1:10, each = 2)), x = 1:20, y = 1:20)
  built <- ggplot_build(solo_tree_plot(tips, NULL, "Homo sapiens"))
  fills <- unique(built$data[[1]]$colour)
  expect_false(any(fills %in% .FATE_COLOUR))
})

test_that("the family subtrees page names only the kinds it draws", {
  families <- data.table(family = c("F1", "F2"), kind = c("with_intact", "no_intact"),
                         n_solo = c(4L, 1L), n_mono = 0L, n_flank = c(2L, 0L))
  tips <- data.table(family = "F1", class = c("SOLO", "SOLO", "FLANK"),
                     x = 1:3, y = 1:3)
  segs <- data.table(family = "F1", x = 0, y = 1, xend = 1, yend = 1)
  p <- family_subtrees_plot(tips, segs, families, "Homo sapiens")
  expect_false(grepl("without", p$labels$title))
  expect_match(p$labels$subtitle, "no family without an intact member")
})
