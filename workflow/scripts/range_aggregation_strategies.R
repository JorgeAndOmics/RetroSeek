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
#   list        — multi-value preserved as CharacterList / list<string>
#   concatenate — single string joined by a separator
#   best        — row with highest `tiebreaker` wins
#   majority    — mode
#   first       — alphabetical first unique
#   strict      — value if unanimous, else a marker
#
# The helpers are pure — they do not read from config themselves; callers
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
aggregate_values <- function(values, strategy,
                             tiebreaker    = NULL,
                             separator     = "; ",
                             strict_marker = "ambiguous") {
  values_chr <- as.character(values)
  uniq <- unique(values_chr)
  switch(
    strategy,
    list        = IRanges::CharacterList(list(uniq)),
    concatenate = paste(sort(uniq), collapse = separator),
    best        = {
      if (is.null(tiebreaker)) {
        stop("aggregate_values(strategy='best') requires a tiebreaker vector.")
      }
      values_chr[which.max(tiebreaker)][1]
    },
    majority    = names(sort(table(values_chr), decreasing = TRUE))[1],
    first       = sort(uniq)[1],
    strict      = if (length(uniq) == 1L) uniq else strict_marker,
    stop("Unknown aggregation strategy: ", strategy)
  )
}


#' Build a tiebreaker-picker closure bound to a configured name.
#'
#' Returns a function that selects one of the candidate vectors based on
#' the tiebreaker name captured at closure construction time. Used by
#' ranges_analysis.R to resolve the configured `best_tiebreaker` once at
#' startup, then call the closure at each aggregation site.
#'
#' @param tb_name One of bitscore | identity | align_length.
#' @return A function(bitscore, identity, align_length = NULL) returning
#'   whichever argument matches `tb_name`.
make_tiebreaker_picker <- function(tb_name = "bitscore") {
  force(tb_name)
  function(bitscore = NULL, identity = NULL, align_length = NULL) {
    switch(
      tb_name,
      bitscore     = bitscore,
      identity     = identity,
      align_length = if (!is.null(align_length)) align_length
                     else stop("best_tiebreaker='align_length' not available in this context yet"),
      bitscore  # default fallback
    )
  }
}
