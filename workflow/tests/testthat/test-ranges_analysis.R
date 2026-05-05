# testthat coverage for workflow/scripts/range_analysis/*.R
#
# Focused on the per-hit `query_coverage` path introduced after the broken
# probe_lengths lookup was retired. Sources only the leaf modules so we
# don't drag in plyranges / GenomicRanges machinery the orchestrator needs.

suppressMessages({
  library(testthat)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(dplyr)
  library(tibble)
})

.script_dir <- file.path("..", "..", "scripts", "range_analysis")
source(file.path(.script_dir, "granges_build.R"))


# Synthetic per-hit BLAST tibble — minimal columns build_blast_gr reads.
.fake_blast_df <- function(viruses, probes, align_len, qstart = 1L, qend = 100L) {
  n <- length(viruses)
  tibble::tibble(
    accession        = rep("chr1", n),
    hsp_sbjct_start  = seq.int(1L, by = 1000L, length.out = n),
    hsp_sbjct_end    = seq.int(1L, by = 1000L, length.out = n) + 99L,
    strand           = rep("+", n),
    label            = rep("Gammaretrovirus", n),
    virus            = viruses,
    hsp_bits         = rep(50, n),
    hsp_identity     = rep(40L, n),
    hsp_align_length = align_len,
    species          = rep("Test_species", n),
    probe            = probes,
    hsp_evalue       = rep(1e-10, n),
    hsp_query_start  = rep(qstart, n),
    hsp_query_end    = rep(qend, n)
  )
}


test_that("build_blast_gr attaches per-hit query_coverage when probe_lengths provided", {
  blast_df <- .fake_blast_df(
    viruses   = c("ALV", "RSV"),
    probes    = c("POL", "POL"),
    align_len = c(100L, 200L)
  )
  probe_lengths <- c("ALV|POL" = 500L, "RSV|POL" = 1000L)
  gr <- build_blast_gr(blast_df, probe_lengths = probe_lengths)
  qcov <- S4Vectors::mcols(gr)$query_coverage

  expect_length(qcov, 2L)
  expect_equal(qcov[1], 100 / 500)   # 0.2
  expect_equal(qcov[2], 200 / 1000)  # 0.2
})


test_that("query_coverage is clamped to [0, 1] when align_length exceeds probe length", {
  blast_df <- .fake_blast_df(
    viruses   = c("ALV"),
    probes    = c("POL"),
    align_len = c(5000L)            # absurdly large vs probe
  )
  probe_lengths <- c("ALV|POL" = 500L)
  gr <- build_blast_gr(blast_df, probe_lengths = probe_lengths)
  expect_equal(S4Vectors::mcols(gr)$query_coverage, 1)
})


test_that("query_coverage is NA for (virus, probe) keys missing from probe_lengths", {
  blast_df <- .fake_blast_df(
    viruses   = c("ALV", "UnknownVirus"),
    probes    = c("POL", "POL"),
    align_len = c(100L, 100L)
  )
  probe_lengths <- c("ALV|POL" = 500L)   # UnknownVirus|POL absent
  gr <- build_blast_gr(blast_df, probe_lengths = probe_lengths)
  qcov <- S4Vectors::mcols(gr)$query_coverage
  expect_equal(qcov[1], 0.2)
  expect_true(is.na(qcov[2]))
})


test_that("build_blast_gr returns NA query_coverage when probe_lengths is NULL", {
  blast_df <- .fake_blast_df(
    viruses = c("ALV"), probes = c("POL"), align_len = c(100L)
  )
  gr <- build_blast_gr(blast_df, probe_lengths = NULL)
  expect_true(all(is.na(S4Vectors::mcols(gr)$query_coverage)))
})


test_that("query_coverage matches the (virus, probe) key on heterogeneous input", {
  # Cross-check the join by mixing several (virus, probe) combos and probe
  # lengths and confirming each hit picks up its own normaliser.
  blast_df <- .fake_blast_df(
    viruses   = c("ALV", "ALV", "RSV", "MMTV"),
    probes    = c("POL", "GAG", "ENV", "PR160"),
    align_len = c(50L,  75L,   60L,   90L)
  )
  probe_lengths <- c(
    "ALV|POL"   = 500L,
    "ALV|GAG"   = 300L,
    "RSV|ENV"   = 400L,
    "MMTV|PR160" = 900L
  )
  gr <- build_blast_gr(blast_df, probe_lengths = probe_lengths)
  qcov <- S4Vectors::mcols(gr)$query_coverage
  expect_equal(qcov, c(50/500, 75/300, 60/400, 90/900))
})
