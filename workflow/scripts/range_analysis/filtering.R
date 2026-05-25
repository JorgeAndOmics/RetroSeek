# -----------------------------------------------------------------------------
# range_analysis / filtering.R
# -----------------------------------------------------------------------------
# Threshold-based filtering of the BLAST GRanges. Per-probe minimum lengths
# are looked up by abbreviation; unknown probes get a fallback of 0 (no length
# filter), matching the pre-refactor behaviour.

suppressMessages({
  library(GenomicRanges)
  library(S4Vectors)
})


# Apply width / bitscore / identity thresholds. Width is checked against the
# probe-specific minimum length from `probe_min_length` (named integer vector).
filter_blast_gr <- function(gr, probe_min_length, bitscore_threshold, identity_threshold) {
  probe_chr <- as.character(S4Vectors::mcols(gr)$probe)
  per_row_min_w <- ifelse(
    !is.na(probe_min_length[probe_chr]),
    probe_min_length[probe_chr],
    0L
  )
  keep <- (BiocGenerics::width(gr) > per_row_min_w) &
    (S4Vectors::mcols(gr)$bitscore > bitscore_threshold) &
    (S4Vectors::mcols(gr)$identity > identity_threshold)
  gr[keep]
}


# Attach a per-row `min_gapwidth` mcols column derived from probe_min_length.
# Used as the gap parameter when reducing overlapping ranges per probe-group:
# a probe-group's reduce uses the unique min_gapwidth across its rows.
attach_min_gapwidth <- function(gr, probe_min_length) {
  probe_chr <- as.character(S4Vectors::mcols(gr)$probe)
  S4Vectors::mcols(gr)$min_gapwidth <- ifelse(
    !is.na(probe_min_length[probe_chr]),
    probe_min_length[probe_chr],
    0L
  )
  gr
}
