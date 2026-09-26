# testthat tests for workflow/scripts/ranges/plot_dataframe.R
#
# Focus: build_plot_dataframe() must resolve Label / Abbreviation under EVERY
# aggregation strategy, not just `concatenate`.
#
# Why this matters: in the plyranges reduce path `aggregate_values(strategy =
# "list")` returns a separator-joined string, not a CharacterList (see the
# "IMPORTANT - why list is concatenated here" note in
# range_aggregation_strategies.R), and ranges/io.R DEFAULTS agg_virus to "list".
# A config that omits the aggregation block therefore lands here. Before this
# fix only the "concatenate" branch split the string, so match() against the
# probe metadata missed and every multi-virus locus silently got NA label and
# NA abbreviation.
#
# Run with: make test-r

suppressMessages({
  library(testthat)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(dplyr)
  library(tibble)
})

.script_dir <- file.path("..", "..", "scripts")
source(file.path(.script_dir, "ranges", "plot_dataframe.R"))

.probe_meta <- tibble::tibble(
  Name         = c("HIV", "HTLV"),
  Label        = c("Lentivirus", "Deltaretrovirus"),
  Abbreviation = c("HIV", "HTLV")
)

# One locus whose `virus` mcol carries two identities, in whichever shape the
# caller's aggregation strategy produced.
.fake_gr <- function(virus) {
  GenomicRanges::GRanges(
    seqnames    = "chr1",
    ranges      = IRanges::IRanges(start = 10, end = 50),
    strand      = "+",
    probe       = "POL",
    virus       = virus,
    species     = "Test_species",
    domain_tier = "domain_selected"
  )
}

test_that("concatenate splits the joined string and resolves both labels", {
  df <- build_plot_dataframe(.fake_gr("HIV; HTLV"), .probe_meta,
                             main_probes = "POL",
                             agg_virus_strategy = "concatenate")
  expect_equal(nrow(df), 2L)
  expect_setequal(df$virus, c("HIV", "HTLV"))
  expect_false(any(is.na(df$label)))
  expect_setequal(df$label, c("Lentivirus", "Deltaretrovirus"))
})

test_that("list (the ranges/io.R default) resolves labels too", {
  df <- build_plot_dataframe(.fake_gr("HIV; HTLV"), .probe_meta,
                             main_probes = "POL",
                             agg_virus_strategy = "list")
  expect_equal(nrow(df), 2L)
  expect_setequal(df$virus, c("HIV", "HTLV"))
  expect_false(any(is.na(df$label)))          # regression: was NA for both rows
  expect_false(any(is.na(df$abbreviation)))
})

test_that("a genuine list-column virus is unnested rather than stringified", {
  gr <- .fake_gr(I(list(c("HIV", "HTLV"))))
  df <- build_plot_dataframe(gr, .probe_meta, main_probes = "POL",
                             agg_virus_strategy = "list")
  expect_equal(nrow(df), 2L)
  expect_setequal(df$virus, c("HIV", "HTLV"))
  expect_false(any(is.na(df$label)))
})

test_that("a single-virus locus is untouched by either strategy", {
  for (strategy in c("best", "concatenate", "list")) {
    df <- build_plot_dataframe(.fake_gr("HIV"), .probe_meta,
                               main_probes = "POL",
                               agg_virus_strategy = strategy)
    expect_equal(nrow(df), 1L, info = strategy)
    expect_equal(df$label, "Lentivirus", info = strategy)
    expect_equal(df$probe_type, "main", info = strategy)
  }
})

# ---------------------------------------------------------------------------
# attach_probe_category: main / accessory / mixed per range, from a probe column
# that is either separator-joined text or a CharacterList.
# ---------------------------------------------------------------------------
.probe_gr <- function(probe) {
  gr <- GenomicRanges::GRanges("c1", IRanges::IRanges(seq_along(probe), width = 1L))
  S4Vectors::mcols(gr)$probe <- probe
  gr
}

test_that("attach_probe_category reads joined text", {
  probe <- c("POL", "REX", "POL; REX", "", "POL; ", NA, "GAG; POL")
  out <- attach_probe_category(.probe_gr(probe), c("POL", "GAG"))
  expect_equal(S4Vectors::mcols(out)$probe_category,
               c("main", "accessory", "mixed", NA, "main", "accessory", "main"))
})

test_that("attach_probe_category reads a CharacterList", {
  probe <- IRanges::CharacterList(list("POL", c("POL", "REX"), character(0),
                                       c("", "REX")))
  out <- attach_probe_category(.probe_gr(probe), c("POL"))
  expect_equal(S4Vectors::mcols(out)$probe_category,
               c("main", "mixed", NA, "accessory"))
})

test_that("attach_probe_category keeps an empty range set typed", {
  out <- attach_probe_category(GenomicRanges::GRanges(), c("POL"))
  expect_identical(S4Vectors::mcols(out)$probe_category, character(0))
})
