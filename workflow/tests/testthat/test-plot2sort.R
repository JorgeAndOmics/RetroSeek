# testthat tests for workflow/scripts/plot2sort.R
#
# Run with:
#   Rscript -e 'testthat::test_dir("workflow/tests/testthat")'
# or via the project Makefile target:
#   make test-r
#
# The script defines pure helpers + plot builders at top level. testthat sources
# the file via the path below; the bottom-of-file `if (sys.nframe() == 0L) main()`
# guard prevents the CLI block from firing during sourcing.

suppressMessages({
  library(testthat)
  library(tibble)
  library(dplyr)
})

source("../../scripts/plot2sort.R")


# Synthetic plot dataframe shaped like ranges_analysis output. Single source of
# truth so each test can mutate / filter as needed without redefining columns.
.fake_plot_df <- function() {
  tibble::tribble(
    ~species, ~probe, ~label,        ~virus, ~abbreviation, ~max_bitscore, ~query_coverage, ~probe_type,
    "S1",     "POL",  "Lentivirus",  "HIV",  "HIV",         300,           0.85,            "main",
    "S1",     "POL",  "Lentivirus",  "HIV",  "HIV",         310,           0.90,            "main",
    "S1",     "POL",  "Lentivirus",  "HIV",  "HIV",         290,           0.80,            "main",
    "S1",     "GAG",  "Lentivirus",  "HIV",  "HIV",         150,           0.55,            "main",
    "S2",     "POL",  "Deltaretro",  "HTLV", "HTLV",        200,           0.65,            "main",
    "S2",     "ENV",  "Deltaretro",  "HTLV", "HTLV",        180,           0.60,            "main",
    "S3",     "VIF",  "Spumavirus",  "FFV",  "FFV",         100,           0.45,            "accessory"
  )
}


# ───────────────────────────── order_by_count ─────────────────────────────

test_that("order_by_count returns levels by count descending, ties alphabetic", {
  df <- tibble::tibble(species = c("Z", "A", "Z", "M", "M", "A", "Z"))
  out <- order_by_count(df, "species")
  expect_equal(out, c("Z", "A", "M"))   # 3, 2, 2 with A < M
})

test_that("order_by_count uses weight column when supplied", {
  df <- tibble::tibble(probe = c("POL", "GAG", "ENV"),
                       count = c(100, 50, 200))
  out <- order_by_count(df, "probe", weight = "count")
  expect_equal(out, c("ENV", "POL", "GAG"))
})

test_that("order_by_count returns empty character on zero-row input", {
  df <- tibble::tibble(species = character(0))
  expect_equal(order_by_count(df, "species"), character(0))
})


# ───────────────────────── collapse_long_tail ─────────────────────────────

test_that("collapse_long_tail is a no-op when top_n is NULL", {
  df <- tibble::tibble(probe = c("POL", "GAG", "ENV", "VIF"))
  expect_identical(collapse_long_tail(df, "probe", top_n = NULL), df)
})

test_that("collapse_long_tail is a no-op when top_n is NA", {
  df <- tibble::tibble(probe = c("POL", "GAG", "ENV", "VIF"))
  expect_identical(collapse_long_tail(df, "probe", top_n = NA), df)
})

test_that("collapse_long_tail is a no-op when top_n >= number of strata", {
  df <- tibble::tibble(probe = c("POL", "GAG", "ENV"))
  expect_identical(collapse_long_tail(df, "probe", top_n = 5), df)
})

test_that("collapse_long_tail keeps top N strata by count, folds the rest", {
  df <- tibble::tibble(
    probe = c("POL", "POL", "POL", "GAG", "GAG", "ENV", "VIF")
  )
  out <- collapse_long_tail(df, "probe", top_n = 2)
  expect_setequal(unique(out$probe), c("POL", "GAG", "Other (2)"))
  expect_equal(sum(out$probe == "Other (2)"), 2)
  expect_equal(sum(out$probe == "POL"),       3)
  expect_equal(sum(out$probe == "GAG"),       2)
})

test_that("collapse_long_tail breaks ties alphabetically (deterministic)", {
  df <- tibble::tibble(probe = c("POL", "GAG", "ENV"))   # all count 1
  out <- collapse_long_tail(df, "probe", top_n = 1)
  # alphabetical first wins → ENV kept; POL + GAG fold
  expect_setequal(unique(out$probe), c("ENV", "Other (2)"))
})

test_that("collapse_long_tail labels the bundle with collapsed count", {
  df <- tibble::tibble(probe = c("A", "B", "C", "D", "E"))
  out <- collapse_long_tail(df, "probe", top_n = 2)
  expect_true("Other (3)" %in% out$probe)
})

test_that("collapse_long_tail honours custom other_label", {
  df <- tibble::tibble(probe = c("A", "B", "C"))
  out <- collapse_long_tail(df, "probe", top_n = 1, other_label = "MISC")
  expect_true(any(grepl("^MISC \\(\\d+\\)$", out$probe)))
})


# ───────────────────────────── empty_plot ─────────────────────────────────

test_that("empty_plot returns a ggplot with the requested title", {
  p <- empty_plot("nothing here")
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "nothing here")
})

test_that("empty_plot forces a white plot.background", {
  # Same reason as add_titles: theme_void()'s default transparent background
  # composes as black in some viewers and hides the diagnostic title text.
  p <- empty_plot("placeholder")
  bg <- p$theme$plot.background
  expect_s3_class(bg, "element_rect")
  expect_equal(bg$fill, "white")
})


# ───────────────────────────── bar_plot ───────────────────────────────────

test_that("bar_plot returns empty placeholder on zero-row input", {
  p <- bar_plot(tibble::tibble(species = character(0),
                               label   = character(0),
                               count   = integer(0)))
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "no data")
})

test_that("bar_plot reorders species levels by total count descending", {
  df <- tibble::tibble(
    species = c("S1", "S2", "S3", "S2"),
    label   = c("L1", "L1", "L1", "L2"),
    count   = c(  5,    3,    1,    4)
  )
  # Totals: S2 = 7, S1 = 5, S3 = 1
  p <- bar_plot(df)
  expect_equal(levels(p$data$species), c("S2", "S1", "S3"))
})


# ─────────────────────── balloon_virus_species_plot ───────────────────────

test_that("balloon plot returns empty placeholder on zero-row input", {
  p <- balloon_virus_species_plot(tibble::tibble(
    species = character(0), virus = character(0), probe = character(0),
    label = character(0), abbreviation = character(0), count = integer(0)))
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "no data")
})

test_that("balloon plot species y-axis ordered by total contribution", {
  df <- tibble::tibble(
    species      = c("S1", "S2", "S3"),
    virus        = c("HIV","HTLV","FFV"),
    probe        = c("POL","POL","POL"),
    label        = c("Lentivirus","Deltaretro","Spumavirus"),
    abbreviation = c("HIV","HTLV","FFV"),
    count        = c(  10,    50,    1)
  )
  # Counts descending: S2 (50), S1 (10), S3 (1).
  # Y axis uses rev(...) so the largest sits at the top → levels c(S3, S1, S2).
  p <- balloon_virus_species_plot(df)
  expect_equal(levels(p$data$species), c("S3", "S1", "S2"))
})


# ─────────────────────────── density / raincloud ──────────────────────────

test_that("density_bitscore_plot empty-input guard fires on 0 rows", {
  p <- density_bitscore_plot(.fake_plot_df()[0, ], q1 = 0, median = 0, q3 = 0)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "no data")
})

test_that("density_bitscore_plot accepts x_scale = 'log10' without erroring", {
  df <- .fake_plot_df()
  expect_silent(density_bitscore_plot(df, q1 = 200, median = 220, q3 = 280,
                                      x_scale = "log10"))
})

test_that("raincloud_bitscore_plot empty-input guard fires on 0 rows", {
  p <- raincloud_bitscore_plot(.fake_plot_df()[0, ])
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "no data")
})

test_that("raincloud_bitscore_plot accepts x_scale = 'log10' without erroring", {
  df <- .fake_plot_df()
  expect_silent(raincloud_bitscore_plot(df, x_scale = "log10"))
})


# ─────────────────────── sankey two-axis ordering + collapse ──────────────

.sankey_input <- function() {
  tibble::tibble(
    species = c("S1", "S2", "S3", "S4", "S5"),
    probe   = c("POL","POL","GAG","ENV","VIF"),
    count   = c(  10,   8,    5,    3,    1)
  )
}

test_that("sankey_species_probe_plot empty-input guard fires on 0 rows", {
  p <- sankey_species_probe_plot(.sankey_input()[0, ])
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "no data")
})

test_that("sankey_species_probe_plot orders species axis by count descending", {
  p <- sankey_species_probe_plot(.sankey_input())
  expect_equal(levels(p$data$species), c("S1", "S2", "S3", "S4", "S5"))
})

test_that("sankey_species_probe_plot orders probe axis by count descending", {
  # POL = 18, GAG = 5, ENV = 3, VIF = 1
  p <- sankey_species_probe_plot(.sankey_input())
  expect_equal(levels(p$data$probe), c("POL", "GAG", "ENV", "VIF"))
})

test_that("sankey_species_probe_plot collapses long tail when top_n is set", {
  p <- sankey_species_probe_plot(.sankey_input(), top_n = 2)
  # Two species kept (S1, S2 by count); three folded into Other (3).
  expect_true("Other (3)" %in% as.character(p$data$species))
  expect_lte(length(unique(p$data$species)), 3L)   # S1 + S2 + Other
})

test_that("sankey_label_probe_plot orders both axes by count", {
  df <- tibble::tibble(
    label = c("L_big","L_big","L_med","L_small"),
    probe = c("POL","GAG","POL","ENV"),
    count = c( 100, 30, 50, 5)
  )
  # Label totals: L_big = 130, L_med = 50, L_small = 5
  # Probe totals: POL = 150, GAG = 30, ENV = 5
  p <- sankey_label_probe_plot(df)
  expect_equal(levels(p$data$label), c("L_big", "L_med", "L_small"))
  expect_equal(levels(p$data$probe), c("POL", "GAG", "ENV"))
})

test_that("sankey_species_label_plot orders both axes by count", {
  df <- tibble::tibble(
    species = c("S_big","S_med","S_small"),
    label   = c("L_alpha","L_alpha","L_beta"),
    count   = c( 100,    50,    1)
  )
  p <- sankey_species_label_plot(df)
  expect_equal(levels(p$data$species), c("S_big", "S_med", "S_small"))
  expect_equal(levels(p$data$label),   c("L_alpha", "L_beta"))
})


# ───────────────────────── new plot: query_coverage ───────────────────────

test_that("query_coverage_plot empty-input guard fires on 0 rows", {
  p <- query_coverage_plot(.fake_plot_df()[0, ])
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "no data")
})

test_that("query_coverage_plot returns placeholder when column missing", {
  df <- .fake_plot_df() %>% dplyr::select(-query_coverage)
  p <- query_coverage_plot(df)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "query_coverage column missing")
})

test_that("query_coverage_plot returns ggplot with x = query_coverage", {
  p <- query_coverage_plot(.fake_plot_df())
  expect_s3_class(p, "ggplot")
  expect_equal(rlang::as_label(p$mapping$x), "query_coverage")
})

test_that("query_coverage_plot returns placeholder when column is all NA", {
  df <- .fake_plot_df()
  df$query_coverage <- NA_real_
  p <- query_coverage_plot(df)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "query_coverage all NA")
})


# ───────────────────────── new plot: heatmap_probe_species ────────────────

test_that("heatmap_probe_species_plot empty-input guard fires on 0 rows", {
  p <- heatmap_probe_species_plot(.fake_plot_df()[0, ])
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "no data")
})

test_that("heatmap reorders both axes by marginal count", {
  df <- .fake_plot_df()
  p  <- heatmap_probe_species_plot(df)
  # Species totals (main+accessory): S1 = 4, S2 = 2, S3 = 1
  # Probe totals: POL = 4, GAG = 1, ENV = 1, VIF = 1 (ties → alphabetic)
  expect_equal(levels(p$data$species), c("S1", "S2", "S3"))
  expect_equal(levels(p$data$probe)[1], "POL")
})

test_that("heatmap fills missing probe×species cells with zero", {
  df <- .fake_plot_df()
  p  <- heatmap_probe_species_plot(df)
  # Original rows: 7. Cells: 3 species × 4 probes = 12.
  expect_equal(nrow(p$data), 12L)
  expect_true(0L %in% p$data$count)
})


# ───────────────────────── new plot: waffle_virus ─────────────────────────

test_that("waffle_virus_plot empty-input guard fires on 0 rows", {
  skip_if_not_installed("waffle")
  p <- waffle_virus_plot(.fake_plot_df()[0, ])
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "no data")
})

test_that("waffle_virus_plot returns ggplot at unit_hits = 1", {
  skip_if_not_installed("waffle")
  p <- waffle_virus_plot(.fake_plot_df(), unit_hits = 1L)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title,   "Hits per virus")
  expect_equal(p$labels$caption, "1 square = 1 hit")
})

test_that("waffle_virus_plot uses unit_hits to bucket squares", {
  skip_if_not_installed("waffle")
  df <- tibble::tibble(virus = rep("HIV", 100))
  p_one  <- waffle_virus_plot(df, unit_hits = 1L)
  p_ten  <- waffle_virus_plot(df, unit_hits = 10L)
  # Caption distinguishes the two configurations.
  expect_equal(p_one$labels$caption, "1 square = 1 hit")
  expect_equal(p_ten$labels$caption, "1 square = 10 hits")
})

test_that("waffle_virus_plot auto-derives unit_hits when input would exceed cap", {
  skip_if_not_installed("waffle")
  # 10,000 hits with a 400-square cap → auto-derive forces unit_hits to 25
  df <- tibble::tibble(virus = rep("HIV", 10000))
  p <- waffle_virus_plot(df, unit_hits = NULL)
  expect_match(p$labels$caption, "auto-scaled from 10000 total hits")
})


# ───────────────────────── title / subtitle injection ─────────────────────

test_that("add_titles prepends subset_label when supplied", {
  p <- ggplot2::ggplot()
  q <- add_titles(p, title = "Foo", subtitle = "Bar", subset_label = "Main")
  expect_equal(q$labels$title,    "Main — Foo")
  expect_equal(q$labels$subtitle, "Bar")
})

test_that("add_titles leaves title untouched when subset_label is NULL or empty", {
  p <- ggplot2::ggplot()
  q1 <- add_titles(p, "Foo", "Bar", subset_label = NULL)
  q2 <- add_titles(p, "Foo", "Bar", subset_label = "")
  expect_equal(q1$labels$title, "Foo")
  expect_equal(q2$labels$title, "Foo")
})

test_that("add_titles forces a white plot.background", {
  # theme_void() makes plot.background transparent, which renders as black in
  # some viewers and hides the title. add_titles must pin the background to
  # white regardless of the underlying theme.
  p <- ggplot2::ggplot() + ggplot2::theme_void()
  q <- add_titles(p, "T", "S")
  bg <- q$theme$plot.background
  expect_s3_class(bg, "element_rect")
  expect_equal(bg$fill, "white")
})

test_that("query_coverage_plot threads subset_label into title", {
  p <- query_coverage_plot(.fake_plot_df(), subset_label = "Accessory")
  expect_equal(p$labels$title, "Accessory — Probe query coverage density")
  expect_match(p$labels$subtitle, "alignment length / probe length")
})


# ───────────────────────────── auto_dims ──────────────────────────────────

test_that("auto_dims keeps the base canvas when n is at or below base_strata", {
  d <- auto_dims(5, axis = "x", base_w = 15, base_h = 12)
  expect_equal(d$w, 15)
  expect_equal(d$h, 12)
})

test_that("auto_dims grows width with n on the x axis", {
  d <- auto_dims(102, axis = "x", base_w = 15, base_h = 12,
                 per_stratum = 0.18, base_strata = 12)
  expect_equal(d$w, 15 + (102 - 12) * 0.18)
  expect_equal(d$h, 12)
})

test_that("auto_dims grows height with n on the y axis", {
  d <- auto_dims(60, axis = "y", base_w = 15, base_h = 12,
                 per_stratum = 0.20, base_strata = 12)
  expect_equal(d$w, 15)
  expect_equal(d$h, 12 + (60 - 12) * 0.20)
})

test_that("auto_dims clamps at the cap", {
  # 1000 strata at 0.18 in/stratum would request ~178 in width — must cap.
  d <- auto_dims(1000, axis = "x", base_w = 15, base_h = 12,
                 per_stratum = 0.18, cap = 60)
  expect_equal(d$w, 60)
})


# ───────────────────────────── intended_dims attribute ────────────────────

test_that("species-axis builders attach an intended_dims attribute", {
  df <- .fake_plot_df() %>% group_count()
  p_bar <- bar_plot(df)
  p_bal <- balloon_virus_species_plot(df)
  p_hm  <- heatmap_probe_species_plot(.fake_plot_df())
  expect_true(is.list(attr(p_bar, "intended_dims")))
  expect_true(is.list(attr(p_bal, "intended_dims")))
  expect_true(is.list(attr(p_hm,  "intended_dims")))
  expect_true(all(c("w", "h") %in% names(attr(p_bar, "intended_dims"))))
})


# ───────────────────────────── save_plot ──────────────────────────────────

test_that("save_plot writes the file and respects dims override", {
  tmp <- tempfile("plot2sort_test_", fileext = "")
  dir.create(tmp)
  p <- empty_plot("hello")
  out <- save_plot("test.png", p, tmp,
                   dims = list(w = 4, h = 3),
                   base_w = 15, base_h = 12, dpi = 72)
  expect_true(file.exists(out))
  expect_match(out, "test\\.png$")
  unlink(tmp, recursive = TRUE)
})

test_that("save_plot falls back to base_w / base_h when dims is NULL", {
  tmp <- tempfile("plot2sort_test_", fileext = "")
  dir.create(tmp)
  p <- empty_plot("hello2")
  out <- save_plot("test2.png", p, tmp,
                   dims = NULL,
                   base_w = 4, base_h = 3, dpi = 72)
  expect_true(file.exists(out))
  unlink(tmp, recursive = TRUE)
})


# ───────────────────────────── verify_required_columns ────────────────────

test_that("verify_required_columns is silent when columns are present", {
  expect_silent(verify_required_columns(.fake_plot_df(),
                                        c("species", "probe", "label")))
})

test_that("verify_required_columns lists every missing column at once", {
  err <- expect_error(
    verify_required_columns(.fake_plot_df(),
                            c("species", "probe", "missing_one", "missing_two"))
  )
  expect_match(conditionMessage(err), "missing_one")
  expect_match(conditionMessage(err), "missing_two")
})


# ─────────────────── multi-value aggregation warning ───────────────────────
#
# aggregation_warning() flags `list` / `concatenate` on virus/label (entry
# explosion); stamp_warning_caption() applies it to a finished plot. Both are
# defined in plot2sort/helpers.R and shared with stage_plot_generator.R.

.cfg_with_agg <- function(virus = "best", label = "best") {
  list(parameters = list(aggregation = list(virus = virus, label = label)))
}

test_that("aggregation_warning returns NULL when virus/label are singular", {
  expect_null(aggregation_warning(.cfg_with_agg("best", "best")))
  expect_null(aggregation_warning(.cfg_with_agg("first", "majority")))
})

test_that("aggregation_warning fires for list / concatenate and names the offender", {
  w_list <- aggregation_warning(.cfg_with_agg("list", "best"))
  expect_type(w_list, "character")
  expect_match(w_list, "virus=list")

  w_concat <- aggregation_warning(.cfg_with_agg("best", "concatenate"))
  expect_match(w_concat, "label=concatenate")

  w_both <- aggregation_warning(.cfg_with_agg("list", "concatenate"))
  expect_match(w_both, "virus=list")
  expect_match(w_both, "label=concatenate")
})

test_that("aggregation_warning returns NULL when there is no aggregation block", {
  expect_null(aggregation_warning(list(parameters = list())))
  expect_null(aggregation_warning(list()))
})

test_that("stamp_warning_caption is a no-op for NULL / empty captions", {
  p <- ggplot2::ggplot()
  expect_identical(stamp_warning_caption(p, NULL), p)
  expect_identical(stamp_warning_caption(p, ""), p)
})

test_that("stamp_warning_caption attaches the caption when supplied", {
  p <- stamp_warning_caption(ggplot2::ggplot(), "watch out")
  expect_true(inherits(p, "ggplot"))
  expect_equal(p$labels$caption, "watch out")
})
