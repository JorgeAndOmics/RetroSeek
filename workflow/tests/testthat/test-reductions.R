# testthat tests for workflow/scripts/range_analysis/reductions.R
#
# Focus areas:
#   * attach_tiebreak_rank() — the per-row deterministic rank that makes the
#     `best` / `majority` aggregation strategies reproducible. It must survive
#     reduce_ranges_directed's internal position re-sort.
#   * reduce_first() determinism — the headline risk: the `best` pick must be
#     identical regardless of input row order.
#   * reduce_first() n_hits — guards the pre-existing `length()` vs `lengths()`
#     bug (n_hits must be the per-merged-range contributor count = M1).
#   * reduce_global() n_loci — the M2 multiplicity column.
#
# Run with: make test-r  (or testthat::test_dir on this directory).

suppressMessages({
  library(testthat)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(plyranges)
  library(dplyr)
})

.script_dir <- file.path("..", "..", "scripts")
source(file.path(.script_dir, "range_aggregation_strategies.R"))
source(file.path(.script_dir, "range_analysis", "reductions.R"))


# Filtered-BLAST-like GRanges carrying exactly the mcols reduce_first reads.
# All rows are probe = POL so reduce_global has a single probe group to merge
# across; virus varies so reduce_first keeps per-virus loci that reduce_global
# then collapses. `order_idx` permutes the rows to probe order-invariance.
.fake_blast_gr <- function(order_idx = NULL) {
  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1", "chr1"),
    ranges   = IRanges::IRanges(start = c(10, 20, 30, 500),
                                end   = c(50, 60, 70, 600)),
    strand   = "+",
    probe          = c("POL", "POL", "POL", "POL"),
    virus          = c("HIV", "HIV", "HTLV", "HIV"),
    label          = c("Lenti_A", "Lenti_B", "Delta_X", "Lenti_C"),
    species        = "Test_species",
    bitscore       = c(100, 100, 50, 80),
    identity       = c(90, 95, 70, 85),
    evalue         = c(1e-10, 1e-12, 1e-5, 1e-8),
    align_length   = c(300L, 320L, 200L, 280L),
    query_coverage = c(0.8, 0.9, 0.5, 0.7),
    min_gapwidth   = 200L
  )
  if (!is.null(order_idx)) gr <- gr[order_idx]
  gr
}

.agg_opts <- function(virus = "best", label = "best") {
  list(
    agg_virus = virus, agg_label = label, agg_species = "first",
    agg_best_tiebreaker  = "bitscore",
    agg_concat_separator = "; ", agg_strict_marker = "ambiguous"
  )
}


# ───────────────────────────── attach_tiebreak_rank ─────────────────────────

test_that("attach_tiebreak_rank adds a per-row rank following the key chain", {
  gr <- attach_tiebreak_rank(.fake_blast_gr(), "bitscore", "identity", "evalue")
  rk <- S4Vectors::mcols(gr)$tiebreak_rank
  expect_length(rk, 4L)
  # rows 1 & 2 tie on bitscore (100); row 2 (qcov 0.9) outranks row 1 (qcov 0.8).
  expect_gt(rk[2], rk[1])
  # row 2 is the chain winner (top bitscore + top query_coverage).
  expect_equal(which.max(rk), 2L)
  # row 3 (bitscore 50) is the chain loser.
  expect_equal(which.min(rk), 3L)
})

test_that("attach_tiebreak_rank gives each logical row the same rank under any input order", {
  rank_of_winner <- function(idx) {
    gr <- attach_tiebreak_rank(.fake_blast_gr(idx), "bitscore", "identity", "evalue")
    mc <- S4Vectors::mcols(gr)
    # The chain winner is the row with bitscore 100 AND query_coverage 0.9.
    mc$tiebreak_rank[mc$bitscore == 100 & mc$query_coverage == 0.9]
  }
  expect_equal(rank_of_winner(c(1, 2, 3, 4)), 4L)
  expect_equal(rank_of_winner(c(4, 3, 2, 1)), 4L)
  expect_equal(rank_of_winner(c(2, 4, 1, 3)), 4L)
})

test_that("attach_tiebreak_rank returns an empty GRanges with the column present", {
  out <- attach_tiebreak_rank(.fake_blast_gr()[FALSE], "bitscore", "identity", "evalue")
  expect_length(out, 0L)
  expect_true("tiebreak_rank" %in% colnames(S4Vectors::mcols(out)))
})


# ─────────────────── reduce_first determinism (the headline risk) ───────────

test_that("reduce_first 'best' label follows the deterministic tie-break chain", {
  gv <- reduce_first(.fake_blast_gr(), "virus", .agg_opts())
  # Rows 1 & 2 (POL/HIV, bitscore tie 100) merge into the chr1:10-60 locus.
  # The tie-break chain breaks the bitscore tie on query_coverage desc, so
  # row 2 (qcov 0.9, label Lenti_B) wins over row 1 (qcov 0.8, label Lenti_A).
  locus <- gv[GenomicRanges::start(gv) == 10L]
  expect_length(locus, 1L)
  expect_equal(S4Vectors::mcols(locus)$label, "Lenti_B")
})

test_that("reduce_first 'best' pick is invariant to input row order", {
  pick_label <- function(idx) {
    gv <- reduce_first(.fake_blast_gr(idx), "virus", .agg_opts())
    as.character(S4Vectors::mcols(gv[GenomicRanges::start(gv) == 10L])$label)
  }
  # If the rank column did not survive reduce_ranges_directed's internal
  # position re-sort, these would diverge — the regression guard for the
  # whole tie-break approach.
  expect_equal(pick_label(c(1, 2, 3, 4)), "Lenti_B")
  expect_equal(pick_label(c(4, 3, 2, 1)), "Lenti_B")
  expect_equal(pick_label(c(2, 4, 1, 3)), "Lenti_B")
})


# ───────────── reduce_first n_hits — per-merged-range contributor count ──────

test_that("reduce_first n_hits counts per-merged-range contributors (lengths, not length)", {
  gv <- reduce_first(.fake_blast_gr(), "virus", .agg_opts())
  # POL/HIV {10-60} has 2 contributing raw hits (rows 1,2); {500-600} has 1
  # (row 4); POL/HTLV {30-70} has 1 (row 3). Total raw hits = 4.
  expect_equal(sort(S4Vectors::mcols(gv)$n_hits), c(1L, 1L, 2L))
  expect_equal(sum(S4Vectors::mcols(gv)$n_hits), 4L)
})


# ───────────────────────── reduce_global n_loci (M2) ────────────────────────

test_that("reduce_global emits an n_loci column", {
  gv <- reduce_first(.fake_blast_gr(), "virus", .agg_opts())
  gg <- reduce_global(gv, .agg_opts())
  expect_true("n_loci" %in% colnames(S4Vectors::mcols(gg)))
})

test_that("n_loci sums to the number of gr_virus loci collapsed", {
  gv <- reduce_first(.fake_blast_gr(), "virus", .agg_opts())
  gg <- reduce_global(gv, .agg_opts())
  # Every gr_virus locus collapses into exactly one gr_global locus, so the
  # total of n_loci across gr_global equals length(gr_virus).
  expect_equal(sum(S4Vectors::mcols(gg)$n_loci), length(gv))
})

test_that("n_loci is distinct from n_hits — at least one global locus collapses >1 locus", {
  gv <- reduce_first(.fake_blast_gr(), "virus", .agg_opts())
  gg <- reduce_global(gv, .agg_opts())
  # The chr1:10-70 global locus merges the POL/HIV and POL/HTLV loci.
  expect_true(any(S4Vectors::mcols(gg)$n_loci > 1L))
  # n_hits rolls up raw hits; n_loci counts loci — n_hits >= n_loci per locus.
  expect_true(all(S4Vectors::mcols(gg)$n_hits >= S4Vectors::mcols(gg)$n_loci))
})
