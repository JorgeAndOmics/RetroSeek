# -------------------
# AGGREGATION STRATEGIES
# -------------------
#
# Reusable helper functions for reducing multi-contributor metadata to a
# single output per merged range in ranges_analysis.R. Extracted into its
# own file so testthat can source just this module without pulling in the
# rest of the pipeline's startup (config read, FASTA load, etc.).
#
# Vocabulary (see docs/configuration.md and docs/adr/ADR-002):
#   list        - multi-value preserved as CharacterList / list<string>
#   concatenate - single string joined by a separator
#   best        - row with highest `tiebreaker` wins
#   majority    - mode
#   first       - alphabetical first unique
#   strict      - value if unanimous, else a marker
#
# The helpers are pure - they do not read from config themselves; callers
# pass strategy name and knobs explicitly. This makes them unit-testable
# without a config fixture.

suppressMessages({
  library(IRanges)
})


#' Reduce a vector of metadata values to a single per-range output.
#'
#' @param values        Character vector of contributor values.
#' @param strategy      One of list | concatenate | best | majority | first | strict.
#' @param tiebreaker    Numeric vector aligned with `values`, required by `best`.
#' @param separator     String separator for `concatenate`. Default "; ".
#' @param strict_marker Value returned by `strict` when contributors disagree.
#' @return Either a single character (for all strategies except `list`) or an
#'   IRanges::CharacterList of length 1 (for `list`).
# Byte order, the same on every machine. Plain sort() follows the locale's
# collation, so `first`, a `majority` tie or a joined list could differ between
# two machines for the same data (found 2026-09-23).
.byte_sort <- function(x) sort(x, method = "radix")

# The most frequent value; ties go to the first in byte order.
.majority <- function(x) {
  counts <- table(factor(x, levels = .byte_sort(unique(x))))
  names(counts)[which.max(counts)]
}

aggregate_values <- function(values, strategy,
                             tiebreaker    = NULL,
                             separator     = "; ",
                             strict_marker = "ambiguous") {
  # plyranges' reduce_ranges_directed presents grouped metadata columns as
  # CompressedAtomicList objects, one top-level element per merged range,
  # each element being the per-row vector of contributing values. The
  # required output is a single aggregated value per merged range - i.e.
  # a vector of length(values), or a CharacterList of length(values) for
  # the `list` strategy. Calling as.character() on a multi-row AtomicList
  # element fails ("top-level elements of length <= 1"), so we must apply
  # the strategy element-wise.
  if (inherits(values, c("AtomicList", "List"))) {
    return(.aggregate_values_listwise(values, strategy, tiebreaker,
                                      separator, strict_marker))
  }
  # Atomic-vector path (unit tests, single-group calls): apply globally.
  values_chr <- as.character(values)
  uniq <- unique(values_chr)
  switch(
    strategy,
    list        = IRanges::CharacterList(list(uniq)),
    concatenate = paste(.byte_sort(uniq), collapse = separator),
    best        = {
      if (is.null(tiebreaker)) {
        stop("aggregate_values(strategy='best') requires a tiebreaker vector.")
      }
      values_chr[which.max(tiebreaker)][1]
    },
    majority    = .majority(values_chr),
    first       = .byte_sort(uniq)[1],
    strict      = if (length(uniq) == 1L) uniq else strict_marker,
    stop("Unknown aggregation strategy: ", strategy)
  )
}


# Element-wise variant for AtomicList input. Returns an atomic character
# vector of length(values) so plyranges can attach the column to the
# per-merged-range output GRanges of reduce_ranges_directed.
#
# IMPORTANT - why "list" is concatenated here:
# plyranges (>=1.22) `summarize_rng` compresses list-columns via
# `Reduce(S4Vectors::pc, ans[[i]])` followed by `as(., "CompressedList")`.
# With a uniform-inner-length CharacterList of length n, that reducer
# collapses to a single concatenated entry instead of preserving the per-
# group structure (verified locally: a 5-element CharacterList of length-1
# entries reduces to length 1, then `stopifnot(length == nr)` fires). To
# avoid the broken path, "list" returns an atomic character with values
# joined by `separator`. Downstream consumers that need a true CharacterList
# can split on `separator` post-reduce - the GFF3 / parquet exporters in
# this pipeline already treat it as a flat string.
#
# When `aggregate_values` is called directly (unit tests, non-plyranges
# callers) the atomic-vector branch above still returns CharacterList for
# strategy="list", preserving the original semantics in that context.
.aggregate_values_listwise <- function(values, strategy, tiebreaker,
                                       separator, strict_marker) {
  n <- length(values)
  vapply(seq_len(n), function(i) {
    elem <- as.character(values[[i]])
    if (length(elem) == 0L) return(NA_character_)
    uniq <- unique(elem)
    switch(
      strategy,
      list        = paste(.byte_sort(uniq), collapse = separator),  # see comment above
      concatenate = paste(.byte_sort(uniq), collapse = separator),
      best        = {
        if (is.null(tiebreaker)) {
          stop("aggregate_values(strategy='best') requires a tiebreaker vector.")
        }
        tb_i <- if (inherits(tiebreaker, c("AtomicList", "List"))) {
          tiebreaker[[i]]
        } else {
          tiebreaker
        }
        elem[which.max(tb_i)][1]
      },
      majority    = .majority(elem),
      first       = .byte_sort(uniq)[1],
      strict      = if (length(uniq) == 1L) uniq else strict_marker,
      stop("Unknown aggregation strategy: ", strategy)
    )
  }, character(1))
}
