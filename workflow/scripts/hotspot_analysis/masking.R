# -----------------------------------------------------------------------------
# hotspot_analysis / masking.R
# -----------------------------------------------------------------------------
# Build the N-run mask for the genome, compute per-window callable bp, and
# pool small scaffolds into a single chromosome stratum so the NB GLM does not
# spend degrees of freedom on singleton-coefficient short contigs.
#
# Pure module: takes Biostrings / GenomicRanges objects in, returns
# GenomicRanges / integer / character vectors out. No I/O.

suppressMessages({
  library(Biostrings)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})


#' Build a GRanges of stretches of N bp in the genome.
#'
#' Wraps `Biostrings::vmatchPattern` over an N-of-length-`mask_size` motif,
#' allowing `mask_mismatch` mismatches per match. The result is reduced so
#' overlapping matches collapse into single intervals (otherwise per-window
#' masked-bp would double-count).
#'
#' Returns an empty GRanges (length 0) when `mask_size == 0`, signalling
#' "no masking" in a way the rest of the pipeline can consume uniformly.
#'
#' Caller contract: `seqs` must already carry normalised chromosome names
#' (set by `load_genome_for_hotspot()`). The resulting mask GRanges then has
#' seqlevels matching the windows produced from the same genome.
build_n_mask <- function(seqs, mask_size, mask_mismatch) {
  if (as.integer(mask_size) == 0L) {
    return(GenomicRanges::GRanges())
  }
  pattern <- Biostrings::DNAString(paste(rep("N", mask_size), collapse = ""))
  match_idx <- Biostrings::vmatchPattern(pattern, seqs, max.mismatch = mask_mismatch)
  mask_gr <- as(match_idx, "GRanges")
  GenomicRanges::reduce(mask_gr)
}


#' Compute the per-window callable bp = window width - masked bp in window.
#'
#' Used as the offset in the NB GLM (log(effective_bp_w)). Implemented via
#' `findOverlaps` + `pintersect` rather than per-window apply, so cost is
#' O(W + M) plus the overlap step, instead of O(W * M).
#'
#' Returns an integer vector of length(windows). Always >= 0.
effective_bp_per_window <- function(windows, mask) {
  widths <- as.integer(BiocGenerics::width(windows))
  if (length(mask) == 0L) {
    return(widths)
  }
  hits <- GenomicRanges::findOverlaps(windows, mask)
  if (length(hits) == 0L) {
    return(widths)
  }
  pieces <- GenomicRanges::pintersect(
    windows[S4Vectors::queryHits(hits)],
    mask[S4Vectors::subjectHits(hits)]
  )
  piece_widths <- as.integer(BiocGenerics::width(pieces))
  qh <- as.integer(S4Vectors::queryHits(hits))
  masked_bp <- rep(0L, length(windows))
  agg <- tapply(piece_widths, qh, sum, simplify = TRUE)
  masked_bp[as.integer(names(agg))] <- as.integer(agg)
  pmax(widths - masked_bp, 0L)
}


#' Map seqlengths -> stratum factor; pool short scaffolds into 'Unplaced'.
#'
#' Per the plan agent's recommendation: scaffolds shorter than
#' `min_factor * window_size` are not given their own intercept in the NB GLM
#' (they have too few windows for theta to be identifiable). Instead they
#' share a single 'Unplaced' stratum.
#'
#' @return named character vector with the same names as `seqlengths`; each
#'   value is either the chromosome name itself (kept as-is) or "Unplaced".
pool_small_scaffolds <- function(seqlengths, window_size, min_factor) {
  threshold <- as.integer(window_size) * as.integer(min_factor)
  stratum <- ifelse(seqlengths >= threshold, names(seqlengths), "Unplaced")
  setNames(stratum, names(seqlengths))
}
