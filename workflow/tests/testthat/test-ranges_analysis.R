# testthat coverage for workflow/scripts/ranges/*.R
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

.script_dir <- file.path("..", "..", "scripts", "ranges")
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


# ---------------------------------------------------------------------------
# annotate_anchored_hits KEEPS every candidate (nothing discarded) and labels
# each with `Parent` (greatest-overlap element — the taxonomic classifier's
# grouping anchor), a per-provirus `domain_tier`, and a per-hit
# `domain_hit_class`. See ranges/validation.R + ADR-009.
# ---------------------------------------------------------------------------
source(file.path(.script_dir, "validation.R"))

# Three elements exercising every domain_tier:
#   retroA — carries a config-matched POL domain      -> domain_selected
#   retroB — carries a protein domain, none config    -> domain_unlisted
#   retroC — no protein domain at all                 -> non_domain
.tier_fixture <- function() {
  retros <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(c(100, 1000, 2000), c(500, 1500, 2500)),
    strand = "+", ID = c("retroA", "retroB", "retroC")
  )
  # config-matched subset (extract_domains_with_probes output): only retroA's POL
  domains_w_probes <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(120, 200), strand = "+",
    Parent = "retroA", probe = "POL"
  )
  # full protein_match superset (extract_all_domains): retroA's + retroB's
  all_domains <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(c(120, 1020), c(200, 1100)), strand = "+",
    Parent = c("retroA", "retroB")
  )
  candidates <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(c(150, 350, 1100, 2100), c(300, 450, 1200, 2200)),
    strand = "+", probe = c("POL", "GAG", "POL", "ENV")
  )
  list(retros = retros, domains_w_probes = domains_w_probes,
       all_domains = all_domains, candidates = candidates)
}

test_that("annotate_anchored_hits keeps all candidates and attaches greatest-overlap Parent", {
  f <- .tier_fixture()
  out <- annotate_anchored_hits(f$candidates, f$retros, f$domains_w_probes, f$all_domains)
  expect_equal(length(out), length(f$candidates))   # nothing discarded
  expect_equal(as.character(S4Vectors::mcols(out)$Parent),
               c("retroA", "retroA", "retroB", "retroC"))
})

test_that("domain_tier is element-wise: selected / unlisted / non_domain", {
  f <- .tier_fixture()
  out <- annotate_anchored_hits(f$candidates, f$retros, f$domains_w_probes, f$all_domains)
  expect_equal(as.character(S4Vectors::mcols(out)$domain_tier),
               c("domain_selected", "domain_selected", "domain_unlisted", "non_domain"))
})

test_that("membership domain_hit_class flags the hit's OWN gene (grain differs from tier)", {
  f <- .tier_fixture()
  out <- annotate_anchored_hits(f$candidates, f$retros, f$domains_w_probes, f$all_domains,
                                hit_domain_mode = "membership")
  # POL hit in retroA -> its gene matches the config POL domain -> substring_match.
  # GAG hit in retroA -> element is domain_selected, but GAG is not the matched
  # gene -> no_substring_match (the two grains legitimately disagree).
  expect_equal(as.character(S4Vectors::mcols(out)$domain_hit_class),
               c("substring_match", "no_substring_match",
                 "no_substring_match", "no_substring_match"))
})

test_that("positional domain_hit_class is co-localization and adds a non_domain level", {
  f <- .tier_fixture()
  out <- annotate_anchored_hits(f$candidates, f$retros, f$domains_w_probes, f$all_domains,
                                hit_domain_mode = "positional")
  # c1 POL(150-300) overlaps the POL config domain(120-200) -> substring_match
  # c2 GAG(350-450) overlaps NO domain                       -> non_domain
  # c3 POL(1100-1200) overlaps retroB's non-config domain    -> no_substring_match
  # c4 ENV(2100-2200) overlaps no domain                     -> non_domain
  expect_equal(as.character(S4Vectors::mcols(out)$domain_hit_class),
               c("substring_match", "non_domain", "no_substring_match", "non_domain"))
})

test_that("annotate_anchored_hits on empty input returns a typed-empty GRanges", {
  f <- .tier_fixture()
  out <- annotate_anchored_hits(f$candidates[FALSE], f$retros, f$domains_w_probes, f$all_domains)
  expect_equal(length(out), 0L)
  expect_true(all(c("Parent", "domain_tier", "domain_hit_class")
                  %in% names(S4Vectors::mcols(out))))
})


# ---------------------------------------------------------------------------
# find_unanchored_hits is the exact strand-aware complement of
# find_candidate_hits: the reduced BLAST hits overlapping NO retrotransposon.
# These are the non-LTR-associated hits recovered into the orphan tier.
# ---------------------------------------------------------------------------
test_that("find_unanchored_hits returns hits overlapping no retrotransposon", {
  retros <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(100, 500), strand = "+", ID = "retroA"
  )
  hits <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(c(150, 1000), c(300, 1100)), strand = "+",
    probe = c("POL", "ENV")
  )
  # hit 1 overlaps retroA (anchored); hit 2 is far away (unanchored)
  un <- find_unanchored_hits(hits, retros)
  expect_equal(length(un), 1L)
  expect_equal(IRanges::start(un), 1000L)
  # complement invariant: candidate + unanchored partition the input exactly
  cand <- find_candidate_hits(hits, retros)
  expect_equal(length(cand) + length(un), length(hits))
})

test_that("find_unanchored_hits keeps everything when there are no retrotransposons", {
  hits <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(c(1, 1000), c(99, 1100)), strand = "+",
    probe = c("POL", "GAG")
  )
  retros <- GenomicRanges::GRanges()
  expect_equal(length(find_unanchored_hits(hits, retros)), 2L)
})

test_that("find_unanchored_hits is strand-aware (opposite-strand retro does not anchor)", {
  retros <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(100, 500), strand = "-", ID = "retroA"
  )
  hits <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(150, 300), strand = "+", probe = "POL"
  )
  # + hit vs - retro: no strand-aware overlap -> the hit is unanchored
  expect_equal(length(find_unanchored_hits(hits, retros)), 1L)
})

test_that("find_unanchored_hits on empty input returns empty", {
  retros <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(100, 500), strand = "+", ID = "retroA"
  )
  expect_equal(length(find_unanchored_hits(GenomicRanges::GRanges(), retros)), 0L)
})


# ---------------------------------------------------------------------------
# cluster_orphan_hits assigns a synthetic Parent by spatial proximity so the
# classifier groups orphans into single multi-gene loci (ADR-010).
# ---------------------------------------------------------------------------
test_that("cluster_orphan_hits groups hits within merge_gap and splits far ones", {
  hits <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(c(100, 500, 20000, 20100), c(200, 600, 20050, 20200)),
    strand = "+", probe = c("POL", "GAG", "ENV", "POL")
  )
  # gap 1&2 = 300 (< 1000 -> same cluster); gap 2&3 = 19400 (> 1000 -> split);
  # gap 3&4 = 50 (< 1000 -> same cluster). Expect two clusters: {1,2} and {3,4}.
  out <- cluster_orphan_hits(hits, merge_gap = 1000L)
  parents <- as.character(S4Vectors::mcols(out)$Parent)
  expect_equal(parents[1], parents[2])
  expect_equal(parents[3], parents[4])
  expect_false(parents[1] == parents[3])
  expect_equal(length(unique(parents)), 2L)
})

test_that("cluster_orphan_hits merges everything under a large gap, splits under a tiny one", {
  hits <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(c(100, 20000), c(200, 20100)), strand = "+",
    probe = c("POL", "GAG")
  )
  expect_equal(length(unique(cluster_orphan_hits(hits, 100000L)$Parent)), 1L)  # one cluster
  expect_equal(length(unique(cluster_orphan_hits(hits, 10L)$Parent)),     2L)  # two clusters
})

test_that("cluster_orphan_hits on empty input returns a typed-empty GRanges with Parent", {
  out <- cluster_orphan_hits(GenomicRanges::GRanges(), 1000L)
  expect_equal(length(out), 0L)
  expect_true("Parent" %in% names(S4Vectors::mcols(out)))
})
