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
  pick_tb <- make_tiebreaker_picker(agg$agg_best_tiebreaker)
  by_label <- identical(merge_option, "label")

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
                                             tiebreaker = pick_tb(bitscore = bitscore, identity = identity, align_length = align_length),
                                             separator = agg$agg_concat_separator,
                                             strict_marker = agg$agg_strict_marker),
          species         = aggregate_values(species, agg$agg_species,
                                             tiebreaker = pick_tb(bitscore = bitscore, identity = identity, align_length = align_length),
                                             separator = agg$agg_concat_separator,
                                             strict_marker = agg$agg_strict_marker),
          n_hits          = length(bitscore),
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
                                             tiebreaker = pick_tb(bitscore = bitscore, identity = identity, align_length = align_length),
                                             separator = agg$agg_concat_separator,
                                             strict_marker = agg$agg_strict_marker),
          species         = aggregate_values(species, agg$agg_species,
                                             tiebreaker = pick_tb(bitscore = bitscore, identity = identity, align_length = align_length),
                                             separator = agg$agg_concat_separator,
                                             strict_marker = agg$agg_strict_marker),
          label           = aggregate_values(label,   agg$agg_label,
                                             tiebreaker = pick_tb(bitscore = bitscore, identity = identity, align_length = align_length),
                                             separator = agg$agg_concat_separator,
                                             strict_marker = agg$agg_strict_marker),
          n_hits          = length(bitscore),
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
  pick_tb <- make_tiebreaker_picker(agg$agg_best_tiebreaker)
  if (length(gr_first) == 0L) return(gr_first)

  reduced <- gr_first %>%
    plyranges::group_by(probe) %>%
    plyranges::reduce_ranges_directed(
      species         = aggregate_values(species, agg$agg_species,
                                         tiebreaker = pick_tb(bitscore = max_bitscore, identity = max_identity),
                                         separator = agg$agg_concat_separator,
                                         strict_marker = agg$agg_strict_marker),
      virus           = aggregate_values(virus,   agg$agg_virus,
                                         tiebreaker = pick_tb(bitscore = max_bitscore, identity = max_identity),
                                         separator = agg$agg_concat_separator,
                                         strict_marker = agg$agg_strict_marker),
      label           = aggregate_values(label,   agg$agg_label,
                                         tiebreaker = pick_tb(bitscore = max_bitscore, identity = max_identity),
                                         separator = agg$agg_concat_separator,
                                         strict_marker = agg$agg_strict_marker),
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
