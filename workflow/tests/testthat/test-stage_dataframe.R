# testthat tests for workflow/scripts/range_analysis/stage_dataframe.R
#
# The three builders that produce the per-genome middle-stage parquets
# consumed by stage_plot_generator.R: spatial concordance + M1 (hits),
# structural completeness + domain composition (ltr), M2 (reduced).

suppressMessages({
  library(testthat)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})

.script_dir <- file.path("..", "..", "scripts", "range_analysis")
source(file.path(.script_dir, "stage_dataframe.R"))


# Synthetic gr_virus: 3 loci with distinct probe IDs and an M1 (n_hits) tally.
.fake_gr_virus <- function() {
  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1"),
    ranges   = IRanges::IRanges(start = c(100, 5000, 50000),
                                end   = c(200, 5100, 50100)),
    strand   = "+"
  )
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    probe          = c("POL", "GAG", "ENV"),
    virus          = c("HIV", "HIV", "HTLV"),
    label          = c("Lentivirus", "Lentivirus", "Deltaretrovirus"),
    species        = "Test_species",
    n_hits         = c(5L, 2L, 1L),
    max_bitscore   = c(300, 150, 90),
    max_identity   = c(95, 80, 70),
    query_coverage = c(0.9, 0.6, 0.4),
    ID             = c("POL_1", "GAG_1", "ENV_1")
  )
  gr
}

# retro_1 overlaps locus 1 (100-200); retro_2 sits ~600 bp from locus 2.
.fake_retros <- function() {
  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges   = IRanges::IRanges(start = c(80, 5700), end = c(260, 6000)),
    strand   = "+"
  )
  S4Vectors::mcols(gr)$type <- "LTR_retrotransposon"
  S4Vectors::mcols(gr)$ID   <- c("retro_1", "retro_2")
  gr
}


# ───────────────────────────── build_stage_hits_df ──────────────────────────

test_that("build_stage_hits_df classifies inside / flanking / disjoint concordance", {
  gv     <- .fake_gr_virus()
  retros <- .fake_retros()
  cand   <- gv[c(1, 2)]   # loci 1 & 2 are candidates
  valid  <- gv[1]         # only locus 1 is valid
  df <- build_stage_hits_df(gv, retros, cand, valid)

  expect_equal(nrow(df), 3L)
  expect_equal(df$concordance[df$probe == "POL"], "inside")    # overlaps retro_1
  expect_equal(df$concordance[df$probe == "GAG"], "flanking")  # ~600 bp from retro_2
  expect_equal(df$concordance[df$probe == "ENV"], "disjoint")  # far from all retros
  expect_equal(df$is_candidate, c(TRUE, TRUE, FALSE))
  expect_equal(df$is_valid,     c(TRUE, FALSE, FALSE))
  expect_equal(df$n_hits, c(5L, 2L, 1L))   # M1 carried through unchanged
})

test_that("build_stage_hits_df returns a typed empty tibble for empty input", {
  empty <- .fake_gr_virus()[FALSE]
  df <- build_stage_hits_df(empty, .fake_retros(), empty, empty)
  expect_equal(nrow(df), 0L)
  expect_true(all(c("concordance", "is_candidate", "is_valid", "n_hits")
                  %in% colnames(df)))
})


# ───────────────────────────── build_stage_ltr_df ───────────────────────────

test_that("build_stage_ltr_df counts flanking LTRs, domains and TSDs per retrotransposon", {
  retros <- .fake_retros()
  flank <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(start = c(80, 240, 5700),
                                end   = c(110, 260, 5740)),
    strand   = "+"
  )
  S4Vectors::mcols(flank)$Parent <- c("retro_1", "retro_1", "retro_2")
  doms <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(start = c(120, 150), end = c(140, 170)),
    strand   = "+"
  )
  S4Vectors::mcols(doms)$Parent <- c("retro_1", "retro_1")
  S4Vectors::mcols(doms)$probe  <- c("POL", "GAG")
  ltr_data <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(start = 78, end = 82), strand = "+"
  )
  S4Vectors::mcols(ltr_data)$type   <- "target_site_duplication"
  S4Vectors::mcols(ltr_data)$Parent <- "retro_1"

  df <- build_stage_ltr_df(retros, flank, doms, ltr_data, .fake_gr_virus())
  expect_equal(nrow(df), 2L)
  r1 <- df[df$ID == "retro_1", ]
  r2 <- df[df$ID == "retro_2", ]
  expect_equal(r1$n_flanking_ltrs, 2L)
  expect_true(r1$has_both_ltrs)
  expect_equal(r1$n_domains, 2L)
  expect_equal(r1$domain_probes, "GAG; POL")   # sorted unique set
  expect_true(r1$has_tsd)
  expect_equal(r2$n_flanking_ltrs, 1L)
  expect_false(r2$has_both_ltrs)
  expect_equal(r2$n_domains, 0L)
  expect_false(r2$has_tsd)
})

test_that("build_stage_ltr_df returns a typed empty tibble for no retrotransposons", {
  empty <- .fake_retros()[FALSE]
  df <- build_stage_ltr_df(empty, empty, empty, empty, .fake_gr_virus())
  expect_equal(nrow(df), 0L)
  expect_true(all(c("n_flanking_ltrs", "has_both_ltrs", "n_domains",
                    "domain_probes", "has_tsd") %in% colnames(df)))
})


# ──────────────────────────── build_stage_reduced_df ────────────────────────

test_that("build_stage_reduced_df carries n_loci (M2) and n_hits (M1 roll-up)", {
  gg <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(start = c(100, 5000), end = c(300, 5100)),
    strand   = "+"
  )
  S4Vectors::mcols(gg) <- S4Vectors::DataFrame(
    probe = c("POL", "GAG"), virus = c("HIV", "HIV"),
    label = c("Lentivirus", "Lentivirus"),
    n_hits = c(7L, 2L), n_loci = c(3L, 1L), max_bitscore = c(300, 150),
    ID = c("POL_1", "GAG_1")
  )
  df <- build_stage_reduced_df(gg)
  expect_equal(nrow(df), 2L)
  expect_equal(df$n_loci, c(3L, 1L))
  expect_equal(df$n_hits, c(7L, 2L))
  expect_true(all(c("n_loci", "n_hits", "probe", "virus") %in% colnames(df)))
})

test_that("build_stage_reduced_df returns a typed empty tibble for empty input", {
  gg <- GenomicRanges::GRanges()
  S4Vectors::mcols(gg) <- S4Vectors::DataFrame(
    probe = character(0), virus = character(0), label = character(0),
    n_hits = integer(0), n_loci = integer(0), max_bitscore = numeric(0),
    ID = character(0)
  )
  df <- build_stage_reduced_df(gg)
  expect_equal(nrow(df), 0L)
  expect_true("n_loci" %in% colnames(df))
})
