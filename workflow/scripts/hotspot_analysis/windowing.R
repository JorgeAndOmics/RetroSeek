# -----------------------------------------------------------------------------
# hotspot_analysis / windowing.R
# -----------------------------------------------------------------------------
# Tile the genome into fixed-width windows, count per-window hits, and
# assemble the per-window tibble that the NB GLM consumes.
#
# Pure module: no I/O, no logging, no config reads. Inputs are GRanges /
# named integer vectors / character vectors; outputs are GRanges and tibbles.

suppressMessages({
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(tibble)
})


#' Tile the genome into non-overlapping windows of `window_size` bp.
#'
#' Last tile per chromosome may be shorter than `window_size` (we set
#' `cut.last.tile.in.chrom = TRUE`); the per-window offset in the NB GLM
#' accounts for that. Gotcha #14 in `.claude/memory/gotchas.md` documents the
#' short-tail behaviour.
tile_genome_for_hotspot <- function(seqlengths, window_size) {
  GenomicRanges::tileGenome(
    seqlengths             = seqlengths,
    tilewidth              = as.integer(window_size),
    cut.last.tile.in.chrom = TRUE
  )
}


#' Count per-window overlaps with a GRanges of hits.
#'
#' Thin wrapper kept here for testability and for the orchestrator's logging
#' to read naturally as a pipeline step.
count_hits_per_window <- function(windows, hits) {
  as.integer(GenomicRanges::countOverlaps(windows, hits))
}


#' Assemble the per-window tibble that drives `models.R::fit_nb_model`.
#'
#' One row per window. Carries chrom (raw chromosome name), chrom_stratum
#' (potentially pooled to "Unplaced" by `pool_small_scaffolds`), start, end,
#' count, effective_bp, and the per-call label. The label column is
#' constant within a single per-label call but is written into the tibble
#' so downstream concatenation across labels keeps everything in one frame.
#'
#' @param windows         GRanges of windows
#' @param counts          integer vector aligned with `windows`
#' @param effective_bp    integer vector aligned with `windows`
#' @param chrom_stratum   named character vector keyed by chromosome name,
#'                        produced by `pool_small_scaffolds()`
#' @param label           single character — the retrovirus genus or
#'                        "Ungrouped" when `hotspot_group_split = false`
assemble_window_table <- function(windows, counts, effective_bp,
                                  chrom_stratum, label) {
  chrom <- as.character(GenomicRanges::seqnames(windows))
  stratum <- unname(chrom_stratum[chrom])
  # Any chrom not present in the stratum map keeps its own name as a fallback
  # so the GLM doesn't drop it via NA factor coercion.
  stratum[is.na(stratum)] <- chrom[is.na(stratum)]
  tibble::tibble(
    chrom        = chrom,
    chrom_stratum = stratum,
    start        = as.integer(BiocGenerics::start(windows)),
    end          = as.integer(BiocGenerics::end(windows)),
    count        = as.integer(counts),
    effective_bp = as.integer(effective_bp),
    label        = as.character(label)
  )
}
