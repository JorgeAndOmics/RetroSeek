# testthat tests for workflow/scripts/range_analysis/erv_assembly.R
#
# Focus areas:
#   * assemble_erv_like() — the chaining of >=2 distinct main-probe loci into
#     composite candidates, the group_by toggle, completeness, canonical-order
#     handling, strand awareness, and determinism.
#   * pure helpers .chain_by_gap / .check_canonical / .observed_order.
#
# Run with: make test-r  (or testthat::test_dir on this directory).

suppressMessages({
  library(testthat)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})

.script_dir <- file.path("..", "..", "scripts")
source(file.path(.script_dir, "range_analysis", "erv_assembly.R"))


# valid_hits-shaped GRanges: one row per main-or-accessory probe locus, carrying
# exactly the mcols assemble_erv_like reads.
.valid_gr <- function(seqnames, start, end, strand, probe, virus,
                      label = "Lab", species = "Test_sp",
                      max_bitscore = 100, max_identity = 90,
                      min_evalue = 1e-10, query_coverage = 0.5, ID = NULL) {
  n <- length(start)
  gr <- GenomicRanges::GRanges(
    seqnames = seqnames,
    ranges   = IRanges::IRanges(start = start, end = end),
    strand   = strand
  )
  if (is.null(ID)) ID <- paste0(probe, "_", seq_len(n))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    probe = probe, virus = virus, label = label, species = species,
    max_bitscore   = rep_len(max_bitscore, n),
    max_identity   = rep_len(max_identity, n),
    min_evalue     = rep_len(min_evalue, n),
    query_coverage = rep_len(query_coverage, n),
    ID = ID
  )
  gr
}

.MAIN <- c("GAG", "POL", "ENV")           # canonical 5'->3' order for tests
.agg  <- list(agg_concat_separator = "; ")

# Convenience wrapper with test-friendly defaults.
.run <- function(gr, main = .MAIN, group_by = "virus", max_gap = 1500,
                 canon = FALSE, thresh = 1.0) {
  assemble_erv_like(gr, main, group_by, max_gap, canon, thresh, .agg)
}


# ───────────────────────── pure helpers ─────────────────────────

test_that(".chain_by_gap chains within gap and breaks beyond it", {
  # gap(2nd) = 50-40-1 = 9 (<=20 -> same); gap(3rd) = 1000-90-1 = 909 (>20 -> new)
  expect_equal(.chain_by_gap(c(1, 50, 1000), c(40, 90, 1100), 20), c(1L, 1L, 2L))
})

test_that(".chain_by_gap bridges via the running max end, not the previous end", {
  # A long first locus (ends 100) keeps later close loci in the same cluster
  # even though the immediate predecessor ends earlier.
  expect_equal(.chain_by_gap(c(1, 30, 60), c(100, 40, 70), 5), c(1L, 1L, 1L))
})

test_that(".observed_order is strand-aware and de-duplicated", {
  expect_equal(.observed_order(c("GAG", "POL", "ENV"), c(100, 300, 500), "+"),
               c("GAG", "POL", "ENV"))
  expect_equal(.observed_order(c("GAG", "POL", "ENV"), c(100, 300, 500), "-"),
               c("ENV", "POL", "GAG"))
  expect_equal(.observed_order(c("GAG", "GAG", "POL"), c(100, 150, 300), "+"),
               c("GAG", "POL"))
})

test_that(".check_canonical follows main_probes order, strand-aware", {
  expect_true(.check_canonical(c("GAG", "POL", "ENV"), .MAIN, "+"))
  expect_false(.check_canonical(c("POL", "GAG", "ENV"), .MAIN, "+"))
  expect_true(.check_canonical(c("GAG", "POL", "ENV"), .MAIN, "-"))   # already 5'->3'
  expect_true(.check_canonical(c("ENV", "POL", "GAG"), .MAIN, "*"))   # decreasing ok on *
  expect_true(.check_canonical(c("GAG"), .MAIN, "+"))                 # trivially ordered
})


# ───────────────────────── composition ─────────────────────────

test_that("a full ERV (all main probes) yields one full candidate with children", {
  gr <- .valid_gr(
    seqnames = "chr1", strand = "+",
    start = c(100, 300, 500), end = c(200, 400, 600),
    probe = c("GAG", "POL", "ENV"), virus = "HIV"
  )
  res <- .run(gr)
  expect_length(res$parents, 1L)
  m <- S4Vectors::mcols(res$parents)
  expect_equal(m$n_main_present, 3L)
  expect_equal(m$completeness_fraction, 1)
  expect_true(m$is_full)
  expect_equal(m$observed_order, "GAG; POL; ENV")
  expect_true(m$is_canonical)
  # Children are the three constituent loci, parented to the candidate.
  expect_length(res$children, 3L)
  expect_true(all(as.character(S4Vectors::mcols(res$children)$Parent) == m$ID))
})

test_that("a partial assembly (GAG+POL, no ENV) is a fragment, not full", {
  gr <- .valid_gr("chr1", c(100, 300), c(200, 400), "+",
                  probe = c("GAG", "POL"), virus = "HIV")
  res <- .run(gr)
  expect_length(res$parents, 1L)
  m <- S4Vectors::mcols(res$parents)
  expect_equal(m$n_main_present, 2L)
  expect_equal(m$completeness_fraction, 2 / 3)
  expect_false(m$is_full)
})

test_that("an isolated single main probe is not assembled", {
  gr <- .valid_gr("chr1", 100, 200, "+", probe = "GAG", virus = "HIV")
  expect_length(.run(gr)$parents, 0L)
})

test_that("main probes farther apart than max_join_distance do not chain", {
  gr <- .valid_gr("chr1", c(100, 100000), c(200, 100100), "+",
                  probe = c("GAG", "POL"), virus = "HIV")
  expect_length(.run(gr, max_gap = 1500)$parents, 0L)
})

test_that("repeats of the same probe do not satisfy the >=2 distinct rule", {
  gr <- .valid_gr("chr1", c(100, 300, 500), c(200, 400, 600), "+",
                  probe = c("GAG", "GAG", "GAG"), virus = "HIV")
  expect_length(.run(gr)$parents, 0L)
})

test_that("accessory probes are ignored entirely (not counted, not children)", {
  # GAG + POL chain; a VIF accessory locus sits between them but must neither
  # count toward distinctness nor appear as a child.
  gr <- .valid_gr("chr1", c(100, 300, 500), c(200, 400, 600), "+",
                  probe = c("GAG", "VIF", "POL"), virus = "HIV")
  res <- .run(gr)
  expect_length(res$parents, 1L)
  m <- S4Vectors::mcols(res$parents)
  expect_equal(m$n_main_present, 2L)
  expect_false(grepl("VIF", m$probes_present))
  expect_length(res$children, 2L)
})


# ───────────────────────── grouping ─────────────────────────

test_that("group_by=virus keeps co-located different-virus candidates separate", {
  gr <- .valid_gr(
    "chr1", c(100, 110, 300, 310), c(200, 210, 400, 410), "+",
    probe = c("GAG", "GAG", "POL", "POL"),
    virus = c("HIV", "HTLV", "HIV", "HTLV")
  )
  res <- .run(gr, group_by = "virus")
  expect_length(res$parents, 2L)
  expect_setequal(as.character(S4Vectors::mcols(res$parents)$virus),
                  c("HIV", "HTLV"))
})

test_that("group_by=none chains across viruses into a single candidate", {
  gr <- .valid_gr(
    "chr1", c(100, 110, 300, 310), c(200, 210, 400, 410), "+",
    probe = c("GAG", "GAG", "POL", "POL"),
    virus = c("HIV", "HTLV", "HIV", "HTLV")
  )
  res <- .run(gr, group_by = "none")
  expect_length(res$parents, 1L)
  expect_equal(S4Vectors::mcols(res$parents)$n_loci, 4L)
})


# ───────────────────────── canonical order ─────────────────────────

test_that("is_canonical is recorded regardless of the require flag", {
  gr <- .valid_gr("chr1", c(100, 300, 500), c(200, 400, 600), "+",
                  probe = c("POL", "GAG", "ENV"), virus = "HIV")   # rearranged
  res <- .run(gr, canon = FALSE)
  expect_length(res$parents, 1L)
  expect_false(S4Vectors::mcols(res$parents)$is_canonical)
  expect_equal(res$dropped_noncanonical, 0L)
})

test_that("require_canonical_order drops rearranged candidates and tallies them", {
  gr <- .valid_gr(
    c("chr1", "chr1", "chr1", "chr2", "chr2", "chr2"),
    c(100, 300, 500, 100, 300, 500), c(200, 400, 600, 200, 400, 600), "+",
    probe = c("GAG", "POL", "ENV", "POL", "GAG", "ENV"),   # chr1 canonical, chr2 rearranged
    virus = "HIV"
  )
  res <- .run(gr, canon = TRUE)
  expect_length(res$parents, 1L)
  expect_equal(as.character(GenomicRanges::seqnames(res$parents)), "chr1")
  expect_equal(res$dropped_noncanonical, 1L)
})


# ───────────────────────── strand ─────────────────────────

test_that("a minus-strand element reads canonical when genes descend in coordinate", {
  # On '-' the 5'->3' direction runs against genomic coordinates: ENV at the
  # low coordinate, GAG at the high one, so observed 5'->3' is GAG,POL,ENV.
  gr <- .valid_gr("chr1", c(100, 300, 500), c(200, 400, 600), "-",
                  probe = c("ENV", "POL", "GAG"), virus = "HIV")
  m <- S4Vectors::mcols(.run(gr)$parents)
  expect_equal(m$observed_order, "GAG; POL; ENV")
  expect_true(m$is_canonical)
})

test_that("opposite-strand loci of the same virus do not chain", {
  gr <- .valid_gr("chr1", c(100, 300), c(200, 400), c("+", "-"),
                  probe = c("GAG", "POL"), virus = "HIV")
  expect_length(.run(gr)$parents, 0L)
})


# ───────────────────────── gap boundary ─────────────────────────

test_that("the max_join_distance boundary is inclusive", {
  # GAG ends at 200. A gap of exactly max_gap (=100) means POL starts at 301.
  chain <- .valid_gr("chr1", c(100, 301), c(200, 401), "+",
                     probe = c("GAG", "POL"), virus = "HIV")
  expect_length(.run(chain, max_gap = 100)$parents, 1L)
  # One bp farther (start 302, gap 101) breaks the chain.
  brk <- .valid_gr("chr1", c(100, 302), c(200, 402), "+",
                   probe = c("GAG", "POL"), virus = "HIV")
  expect_length(.run(brk, max_gap = 100)$parents, 0L)
})


# ───────────────────────── multi-value virus ─────────────────────────

test_that("a multi-virus locus explodes into each virus group", {
  # GAG belongs to both HIV and HTLV; POL only to HIV. Under group_by=virus the
  # HIV group sees GAG+POL (a candidate); the HTLV group sees only GAG (none).
  gr <- .valid_gr("chr1", c(100, 300), c(200, 400), "+",
                  probe = c("GAG", "POL"), virus = c("HIV; HTLV", "HIV"))
  res <- .run(gr, group_by = "virus")
  expect_length(res$parents, 1L)
  expect_equal(as.character(S4Vectors::mcols(res$parents)$virus), "HIV")
})


# ───────────────────────── degenerate ─────────────────────────

test_that("empty and single-locus inputs return typed-empty results", {
  empty <- .run(.valid_gr("chr1", 100, 200, "+", probe = "GAG", virus = "HIV")[0])
  expect_length(empty$parents, 0L)
  expect_length(empty$children, 0L)
  expect_equal(empty$dropped_noncanonical, 0L)
  expect_true("completeness_fraction" %in% colnames(S4Vectors::mcols(empty$parents)))
})


# ───────────────────────── determinism ─────────────────────────

test_that("candidate IDs, coordinates and order are invariant to input row order", {
  base <- .valid_gr(
    c("chr1", "chr1", "chr1", "chr2", "chr2"),
    c(100, 300, 500, 100, 300), c(200, 400, 600, 200, 400), "+",
    probe = c("GAG", "POL", "ENV", "GAG", "POL"),
    virus = "HIV"
  )
  df <- function(idx) build_erv_like_df(.run(base[idx])$parents)
  d1 <- df(c(1, 2, 3, 4, 5))
  d2 <- df(c(5, 4, 3, 2, 1))
  d3 <- df(c(3, 1, 5, 2, 4))
  expect_identical(d1, d2)
  expect_identical(d1, d3)
  # Two candidates: a full one on chr1, a fragment on chr2, IDs by position.
  expect_equal(d1$ID, c("ERV_like_1", "ERV_like_2"))
  expect_equal(d1$n_main_present, c(3L, 2L))
})


# ───────────────────────── build_erv_like_df ─────────────────────────

test_that("build_erv_like_df returns a typed-empty tibble for no candidates", {
  out <- build_erv_like_df(.empty_erv_parents())
  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 0L)
  expect_true(all(c("ID", "completeness_fraction", "is_full", "is_canonical")
                  %in% names(out)))
})


# ───────────────────────── build_erv_like_members_df ──────────────────────

test_that("build_erv_like_members_df computes per-candidate running-max gaps", {
  gr <- .valid_gr("chr1", c(100, 300, 700), c(200, 400, 800), "+",
                  probe = c("GAG", "POL", "ENV"), virus = "HIV")
  res <- .run(gr)
  m <- build_erv_like_members_df(res$children)
  expect_equal(nrow(m), 3L)
  # First member of the candidate has NA; subsequent gaps are start - prev end - 1.
  m <- m[order(m$start), ]
  expect_true(is.na(m$gap_to_prev[1]))
  expect_equal(m$gap_to_prev[2], 300L - 200L - 1L)   # 99
  expect_equal(m$gap_to_prev[3], 700L - 400L - 1L)   # 299 (running max end = 400)
})

test_that("build_erv_like_members_df returns typed-empty on no children", {
  out <- build_erv_like_members_df(.empty_erv_children())
  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 0L)
  expect_true(all(c("Parent", "gap_to_prev", "probe") %in% names(out)))
})
