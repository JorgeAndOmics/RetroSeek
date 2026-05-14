# -----------------------------------------------------------------------------
# range_analysis / reductions.R
# -----------------------------------------------------------------------------
# Two-stage reduction of the filtered BLAST GRanges:
#   1. `reduce_first(gr, ...)`: merge overlapping ranges per (probe, virus)
#      or (probe, label) group, depending on `merge_option`. Composite
#      metrics replace the pre-refactor mean_bitscore / mean_identity.
#   2. `reduce_global(gr_first, ...)`: further merge per probe across the
#      virus / label groupings, re-aggregating the composite metrics.
#
# The aggregation strategies (concatenate / list / first / etc.) are dispatched
# via `aggregate_values()` from range_aggregation_strategies.R, which handles
# both flat-character and AtomicList input shapes that plyranges can present.

suppressMessages({
  library(GenomicRanges)
  library(plyranges)
  library(IRanges)
  library(S4Vectors)
  library(dplyr)
})


# Attach a per-row `tiebreak_rank` integer column encoding the full
# deterministic key chain as one monotonic value (higher = more preferred).
#
# Why a per-row rank instead of sorting the input GRanges: plyranges'
# `reduce_ranges_directed` re-sorts each merged range's contributing rows by
# genomic position internally, so a pre-sort of the input is undone before
# `aggregate_values` sees the rows. A per-row rank rides along through that
# reordering as ordinary mcols data — `aggregate_values(strategy = "best")`
# then picks `which.max(tiebreak_rank)`, the chain-winning contributor,
# deterministically and reproducibly across R / plyranges versions.
#
# Key chain: primary_col desc -> query_coverage desc -> identity_col desc ->
# evalue_col asc -> seqnames -> start -> end -> label name asc. `primary_col`
# is the configured best_tiebreaker; the remainder is the fixed deterministic
# tail. Metric column names differ between reduction stages, so they are
# passed by name.
attach_tiebreak_rank <- function(gr, primary_col, identity_col, evalue_col) {
  if (length(gr) == 0L) {
    S4Vectors::mcols(gr)$tiebreak_rank <- integer(0)
    return(gr)
  }
  mc <- S4Vectors::mcols(gr)
  ord <- order(
    -as.numeric(mc[[primary_col]]),
    -as.numeric(mc$query_coverage),
    -as.numeric(mc[[identity_col]]),
     as.numeric(mc[[evalue_col]]),
     as.character(GenomicRanges::seqnames(gr)),
     GenomicRanges::start(gr),
     GenomicRanges::end(gr),
     as.character(mc$label),
    method = "radix"
  )
  rank <- integer(length(ord))
  rank[ord] <- rev(seq_along(ord))   # order-first row -> highest rank value
  S4Vectors::mcols(gr)$tiebreak_rank <- rank
  gr
}


# Resolve the configured `best_tiebreaker` name to the actual mcols column at
# each reduction stage. reduce_first sees the raw per-hit columns; reduce_global
# sees only the first-reduction composites — and `align_length` is gone by
# then, so it falls back to the bitscore composite (graceful, unlike the
# pre-refactor make_tiebreaker_picker which errored).
.primary_tiebreak_col <- function(best_tiebreaker, stage = c("first", "global")) {
  stage <- match.arg(stage)
  if (identical(stage, "first")) {
    switch(best_tiebreaker,
           bitscore     = "bitscore",
           identity     = "identity",
           align_length = "align_length",
           "bitscore")
  } else {
    switch(best_tiebreaker,
           bitscore = "max_bitscore",
           identity = "max_identity",
           "max_bitscore")
  }
}


# First reduction: per (probe, virus) or (probe, label) group, merge overlapping
# ranges and compute composite metrics from contributing rows.
#
# Args:
#   gr               : BLAST GRanges (filtered), with `min_gapwidth` and
#                      per-hit `query_coverage` mcols already attached by
#                      build_blast_gr.
#   merge_option     : "virus" or "label" (selects the secondary group key).
#   agg              : aggregation-options list from read_pipeline_options().
reduce_first <- function(gr, merge_option, agg) {
  by_label <- identical(merge_option, "label")

  # Per-row deterministic tie-break rank (see attach_tiebreak_rank). `best`
  # aggregation reads it via which.max, so the pick is reproducible regardless
  # of how reduce_ranges_directed reorders contributors internally.
  primary_col <- .primary_tiebreak_col(agg$agg_best_tiebreaker, "first")
  gr <- attach_tiebreak_rank(gr, primary_col, "identity", "evalue")

  gr_list <- BiocGenerics::Reduce(
    function(acc, sub_gr) c(acc, list(sub_gr)),
    split(gr, ~ probe),
    accumulate = FALSE,
    init = list()
  )

  reduce_one <- function(sub_gr) {
    if (length(sub_gr) == 0L) return(sub_gr)
    gap_val <- unique(S4Vectors::mcols(sub_gr)$min_gapwidth)
    if (length(gap_val) != 1L) gap_val <- gap_val[1]

    if (by_label) {
      sub_gr %>%
        plyranges::group_by(probe, label) %>%
        plyranges::reduce_ranges_directed(
          min.gapwidth    = gap_val,
          virus           = aggregate_values(virus,   agg$agg_virus,
                                             tiebreaker = tiebreak_rank,
                                             separator = agg$agg_concat_separator,
                                             strict_marker = agg$agg_strict_marker),
          species         = aggregate_values(species, agg$agg_species,
                                             tiebreaker = tiebreak_rank,
                                             separator = agg$agg_concat_separator,
                                             strict_marker = agg$agg_strict_marker),
          # lengths() (plural): per-merged-range contributor count. length()
          # (singular) returns the CompressedNumericList element count — the
          # number of merged ranges, not the per-range hit tally.
          n_hits          = lengths(bitscore),
          max_bitscore    = max(bitscore),
          sum_bitscore    = sum(bitscore),
          median_bitscore = median(bitscore),
          max_identity    = max(identity),
          min_evalue      = min(evalue),
          query_coverage  = pmin(1, sum(query_coverage))
        )
    } else {
      sub_gr %>%
        plyranges::group_by(probe, virus) %>%
        plyranges::reduce_ranges_directed(
          min.gapwidth    = gap_val,
          virus           = aggregate_values(virus,   agg$agg_virus,
                                             tiebreaker = tiebreak_rank,
                                             separator = agg$agg_concat_separator,
                                             strict_marker = agg$agg_strict_marker),
          species         = aggregate_values(species, agg$agg_species,
                                             tiebreaker = tiebreak_rank,
                                             separator = agg$agg_concat_separator,
                                             strict_marker = agg$agg_strict_marker),
          label           = aggregate_values(label,   agg$agg_label,
                                             tiebreaker = tiebreak_rank,
                                             separator = agg$agg_concat_separator,
                                             strict_marker = agg$agg_strict_marker),
          n_hits          = lengths(bitscore),
          max_bitscore    = max(bitscore),
          sum_bitscore    = sum(bitscore),
          median_bitscore = median(bitscore),
          max_identity    = max(identity),
          min_evalue      = min(evalue),
          query_coverage  = pmin(1, sum(query_coverage))
        )
    }
  }

  reduced_list <- lapply(gr_list, reduce_one)
  # Drop empty results so bind_ranges doesn't choke.
  reduced_list <- reduced_list[vapply(reduced_list, length, integer(1)) > 0L]
  if (length(reduced_list) == 0L) return(gr[FALSE])
  plyranges::bind_ranges(reduced_list)
}


# Global reduction: merge per probe across the virus/label groupings of the
# first reduction. Aggregation is applied element-wise, composite metrics are
# re-aggregated. The output is intentionally not arrange()d — plyranges'
# arrange.Ranges hits an internal Reduce length-mismatch on the
# CompressedAtomicList output of reduce_ranges_directed under some shapes,
# and downstream consumers (rtracklayer::export, overlap queries) don't
# require a specific order. The group_by state is dropped via ungroup() so
# subsequent operations don't carry it implicitly.
reduce_global <- function(gr_first, agg) {
  if (length(gr_first) == 0L) return(gr_first)

  # Per-row deterministic tie-break rank, same mechanism as reduce_first. The
  # input columns here are the first-reduction composites (max_bitscore etc.).
  primary_col <- .primary_tiebreak_col(agg$agg_best_tiebreaker, "global")
  gr_first <- attach_tiebreak_rank(gr_first, primary_col,
                                   "max_identity", "min_evalue")

  reduced <- gr_first %>%
    plyranges::group_by(probe) %>%
    plyranges::reduce_ranges_directed(
      species         = aggregate_values(species, agg$agg_species,
                                         tiebreaker = tiebreak_rank,
                                         separator = agg$agg_concat_separator,
                                         strict_marker = agg$agg_strict_marker),
      virus           = aggregate_values(virus,   agg$agg_virus,
                                         tiebreaker = tiebreak_rank,
                                         separator = agg$agg_concat_separator,
                                         strict_marker = agg$agg_strict_marker),
      label           = aggregate_values(label,   agg$agg_label,
                                         tiebreaker = tiebreak_rank,
                                         separator = agg$agg_concat_separator,
                                         strict_marker = agg$agg_strict_marker),
      # n_loci = per-merged-range count of contributing gr_virus loci (M2),
      # via lengths() (plural) on the CompressedNumericList. n_hits rolls the
      # M1 tallies up via sum(). Both read the input list independently.
      n_loci          = lengths(n_hits),
      n_hits          = sum(n_hits),
      max_bitscore    = max(max_bitscore),
      sum_bitscore    = sum(sum_bitscore),
      median_bitscore = median(median_bitscore),
      max_identity    = max(max_identity),
      min_evalue      = min(min_evalue),
      query_coverage  = pmin(1, sum(query_coverage)),
      type            = "proviral_sequence"
    )

  dplyr::ungroup(reduced)
}


# Attach a probe-scoped sequential ID column ("ENV_1", "ENV_2", ...).
attach_probe_id <- function(gr) {
  if (length(gr) == 0L) {
    S4Vectors::mcols(gr)$ID <- character(0)
    return(gr)
  }
  gr %>%
    plyranges::group_by(probe) %>%
    dplyr::mutate(ID = paste0(probe, "_", seq_along(probe))) %>%
    dplyr::ungroup()
}
