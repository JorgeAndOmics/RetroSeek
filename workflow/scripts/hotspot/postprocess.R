# -----------------------------------------------------------------------------
# hotspot / postprocess.R
# -----------------------------------------------------------------------------
# Convert per-window scored tibbles into hotspot regions:
#   1. Select windows whose qval is below threshold.
#   2. Merge adjacent (or near-adjacent within `gap` bp) significant windows
#      into a single region; sum counts and effective_bp across the span.
#   3. Recompute a merged-region p-value by re-evaluating the fitted NB at
#      summed (count, effective_bp). Conservative under spatial correlation
#      (sums of correlated NB counts have higher variance than the model
#      assumes), which is the right direction for a discovery filter. Do NOT
#      BH-adjust again - this is a post-hoc summary, not a new test family.
#   4. Apply min-hits filter AFTER merging - a singleton tile next to a
#      3-hit tile is a legitimate 4-hit hotspot.
#
# Pure module. Inputs: scored tibble + fitted model. Output: GRanges.

suppressMessages({
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(dplyr)
  library(tibble)
})


# Canonical empty hotspot GRanges with the full mcols schema. All
# postprocess paths return this on the no-data case so that concatenating
# per-label hotspot GRanges (`do.call(c, list)`) succeeds - `BiocGenerics::c`
# requires matching mcols columns across operands.
.empty_merged_gr <- function() {
  gr <- GenomicRanges::GRanges()
  S4Vectors::mcols(gr)$label          <- character(0)
  S4Vectors::mcols(gr)$count          <- integer(0)
  S4Vectors::mcols(gr)$effective_bp   <- integer(0)
  S4Vectors::mcols(gr)$n_windows      <- integer(0)
  S4Vectors::mcols(gr)$chrom_stratum  <- character(0)
  S4Vectors::mcols(gr)$mu_nb_region   <- numeric(0)
  S4Vectors::mcols(gr)$pval_nb_region <- numeric(0)
  gr
}


#' Select windows below the q-value threshold.
#'
#' Skips windows with NA qval (effective_bp == 0 etc).
select_significant_windows <- function(window_df, threshold) {
  dplyr::filter(window_df, !is.na(.data$qval_nb), .data$qval_nb < threshold)
}


#' Merge adjacent significant windows into hotspot regions.
#'
#' Uses `GenomicRanges::reduce(min.gapwidth = gap + 1L)`: gap == 0 -> only
#' strictly adjacent windows merge; gap == window_size -> windows separated
#' by up to one window also merge. Aggregates `count`, `effective_bp`, and
#' `n_windows` per merged region.
#'
#' Returns a GRanges with mcols: chrom_stratum, label, count, effective_bp,
#' n_windows, hotspot_id (assigned in caller order).
merge_adjacent_hotspots <- function(significant_df, gap = 0L) {
  if (nrow(significant_df) == 0L) {
    return(.empty_merged_gr())
  }
  gr <- GenomicRanges::GRanges(
    seqnames = significant_df$chrom,
    ranges   = IRanges::IRanges(start = significant_df$start,
                                end   = significant_df$end)
  )
  merged <- GenomicRanges::reduce(gr, min.gapwidth = as.integer(gap) + 1L)

  # Map each input window back to its merged region. `reduce()` produces a
  # disjoint partition over input regions, so every query has exactly one hit.
  hits <- GenomicRanges::findOverlaps(gr, merged)
  region_idx <- as.integer(S4Vectors::subjectHits(hits))
  query_idx  <- as.integer(S4Vectors::queryHits(hits))

  agg <- tibble::tibble(
    region_idx    = region_idx,
    count         = significant_df$count[query_idx],
    effective_bp  = significant_df$effective_bp[query_idx],
    chrom_stratum = significant_df$chrom_stratum[query_idx],
    label         = significant_df$label[query_idx]
  ) %>%
    dplyr::group_by(.data$region_idx) %>%
    dplyr::summarise(
      count         = sum(.data$count),
      effective_bp  = sum(.data$effective_bp),
      n_windows     = dplyr::n(),
      chrom_stratum = dplyr::first(.data$chrom_stratum),
      label         = dplyr::first(.data$label),
      .groups       = "drop"
    ) %>%
    dplyr::arrange(.data$region_idx)

  S4Vectors::mcols(merged)$label          <- agg$label
  S4Vectors::mcols(merged)$count          <- as.integer(agg$count)
  S4Vectors::mcols(merged)$effective_bp   <- as.integer(agg$effective_bp)
  S4Vectors::mcols(merged)$n_windows      <- as.integer(agg$n_windows)
  S4Vectors::mcols(merged)$chrom_stratum  <- agg$chrom_stratum
  # Placeholder columns; recompute_merged_pvalue overwrites these. Adding
  # them here keeps the mcols schema invariant across all postprocess paths.
  S4Vectors::mcols(merged)$mu_nb_region   <- rep(NA_real_, length(merged))
  S4Vectors::mcols(merged)$pval_nb_region <- rep(NA_real_, length(merged))
  merged
}


#' Recompute the per-region p-value by re-evaluating the fitted NB at the
#' summed (count, effective_bp) of each merged region.
#'
#' Per the plan agent's recommendation (point D), this is more honest than
#' Fisher / Stouffer combination of per-window p-values - those assume
#' independence, which is false for adjacent windows. The resulting p is
#' conservative because spatial correlation inflates the variance.
#'
#' Adds `mu_nb_region` and `pval_nb_region` to the GRanges mcols. Does NOT
#' add `qval_nb_region` - these are post-hoc summaries, not new tests.
recompute_merged_pvalue <- function(merged_gr, fit) {
  if (length(merged_gr) == 0L) {
    return(merged_gr)
  }
  m <- S4Vectors::mcols(merged_gr)
  newdata <- tibble::tibble(
    chrom         = as.character(GenomicRanges::seqnames(merged_gr)),
    chrom_stratum = as.character(m$chrom_stratum),
    start         = as.integer(BiocGenerics::start(merged_gr)),
    end           = as.integer(BiocGenerics::end(merged_gr)),
    count         = as.integer(m$count),
    effective_bp  = as.integer(m$effective_bp),
    label         = as.character(m$label)
  )
  scored <- score_windows_nb(newdata, fit)
  S4Vectors::mcols(merged_gr)$mu_nb_region   <- scored$mu_nb
  S4Vectors::mcols(merged_gr)$pval_nb_region <- scored$pval_nb
  merged_gr
}


#' Drop hotspot regions whose total count is below `min_hits`.
#'
#' Applied AFTER merging - see plan-agent point E. Threshold is on the
#' aggregated count across the merged span, not per-window.
apply_min_hits_filter <- function(merged_gr, min_hits) {
  if (length(merged_gr) == 0L || is.null(min_hits) || as.integer(min_hits) <= 0L) {
    return(merged_gr)
  }
  keep <- as.integer(S4Vectors::mcols(merged_gr)$count) >= as.integer(min_hits)
  merged_gr[keep]
}


#' Assign sequential hotspot IDs across a (post-merge, post-filter) GRanges.
#'
#' Caller passes a `species` prefix so cross-genome aggregation later can
#' tell IDs apart without reassignment.
assign_hotspot_ids <- function(merged_gr, species) {
  if (length(merged_gr) == 0L) {
    S4Vectors::mcols(merged_gr)$hotspot_id <- character(0)
    return(merged_gr)
  }
  ids <- sprintf("%s_HS_%05d", species, seq_along(merged_gr))
  S4Vectors::mcols(merged_gr)$hotspot_id <- ids
  merged_gr
}


#' Attach hotspot_id back to the per-window scored tibble for the CSV export.
#'
#' Each per-window row gets the id of the merged region it falls inside, or
#' NA if it isn't part of any hotspot.
attach_hotspot_id_to_windows <- function(window_df, merged_gr) {
  if (nrow(window_df) == 0L) {
    return(dplyr::mutate(window_df, hotspot_id = character(0)))
  }
  hotspot_id_col <- rep(NA_character_, nrow(window_df))
  if (length(merged_gr) > 0L) {
    win_gr <- GenomicRanges::GRanges(
      seqnames = window_df$chrom,
      ranges   = IRanges::IRanges(start = window_df$start, end = window_df$end)
    )
    hits <- GenomicRanges::findOverlaps(win_gr, merged_gr)
    if (length(hits) > 0L) {
      ids <- as.character(S4Vectors::mcols(merged_gr)$hotspot_id)
      hotspot_id_col[as.integer(S4Vectors::queryHits(hits))] <-
        ids[as.integer(S4Vectors::subjectHits(hits))]
    }
  }
  dplyr::mutate(window_df, hotspot_id = hotspot_id_col)
}


#' Annotate merged hotspot regions with the composition of the loci inside them
#' (ADR-012).
#'
#' Detection deliberately runs on EVERY locus in the chosen tier: splitting an
#' already sparse count matrix by structural class would starve the NB fit.
#' Instead each called region is described by what it is made of, so a hotspot
#' can be read as intact-provirus-driven or fragment-driven without a second
#' query. This is annotation, not filtering - no region is added or removed.
#'
#' Adds per region:
#'   n_loci                     loci overlapping the region
#'   n_full / n_partial / n_gene   structure_class breakdown
#'   n_ltr_flanked / n_orphan   tier breakdown
#'   dominant_taxon             most frequent value of `group_col` (ties -> first
#'                              alphabetically, so the output is deterministic)
#'   mean_confidence            mean numeric confidence, NA when unavailable
#'
#' Columns absent from `loci` yield zero counts / NA rather than an error, so a
#' raw-hit input (which carries none of them) still passes through unharmed.
#'
#' @param merged_gr GRanges of merged hotspot regions
#' @param loci      GRanges of the loci that were counted
#' @param group_col mcols column naming the lineage (e.g. "segment")
annotate_hotspot_composition <- function(merged_gr, loci, group_col = "segment") {
  n <- length(merged_gr)
  .rep <- function(value) rep(value, n)
  S4Vectors::mcols(merged_gr)$n_loci          <- .rep(0L)
  S4Vectors::mcols(merged_gr)$n_full          <- .rep(0L)
  S4Vectors::mcols(merged_gr)$n_partial       <- .rep(0L)
  S4Vectors::mcols(merged_gr)$n_gene          <- .rep(0L)
  S4Vectors::mcols(merged_gr)$n_ltr_flanked   <- .rep(0L)
  S4Vectors::mcols(merged_gr)$n_orphan        <- .rep(0L)
  S4Vectors::mcols(merged_gr)$dominant_taxon  <- .rep(NA_character_)
  S4Vectors::mcols(merged_gr)$mean_confidence <- .rep(NA_real_)
  if (n == 0L || length(loci) == 0L) return(merged_gr)

  ov <- GenomicRanges::findOverlaps(loci, merged_gr, ignore.strand = TRUE)
  if (length(ov) == 0L) return(merged_gr)
  locus_i  <- as.integer(S4Vectors::queryHits(ov))
  region_i <- as.integer(S4Vectors::subjectHits(ov))

  mc <- S4Vectors::mcols(loci)
  col_of <- function(name) {
    if (name %in% colnames(mc)) as.character(mc[[name]][locus_i]) else NULL
  }
  structure_class <- col_of("structure_class")
  source_col      <- col_of("source")
  taxon           <- col_of(group_col)
  conf <- if ("confidence" %in% colnames(mc)) {
    suppressWarnings(as.numeric(mc$confidence[locus_i]))
  } else {
    NULL
  }

  # tapply over the region index keeps this a single pass per statistic.
  tally <- function(mask) {
    out <- integer(n)
    if (is.null(mask)) return(out)
    t <- table(region_i[mask])
    out[as.integer(names(t))] <- as.integer(t)
    out
  }
  S4Vectors::mcols(merged_gr)$n_loci <- tally(rep(TRUE, length(region_i)))
  if (!is.null(structure_class)) {
    S4Vectors::mcols(merged_gr)$n_full    <- tally(structure_class == "full")
    S4Vectors::mcols(merged_gr)$n_partial <- tally(structure_class == "partial")
    S4Vectors::mcols(merged_gr)$n_gene    <- tally(structure_class == "gene")
  }
  if (!is.null(source_col)) {
    S4Vectors::mcols(merged_gr)$n_ltr_flanked <- tally(source_col == "ltr-flanked")
    S4Vectors::mcols(merged_gr)$n_orphan      <- tally(source_col == "orphan")
  }
  if (!is.null(taxon)) {
    dom <- vapply(split(taxon, region_i), function(v) {
      v <- v[!is.na(v) & nzchar(v)]
      if (length(v) == 0L) return(NA_character_)
      t <- table(v)
      # sort() breaks count ties alphabetically -> deterministic output
      names(sort(t, decreasing = TRUE))[1]
    }, character(1))
    S4Vectors::mcols(merged_gr)$dominant_taxon[as.integer(names(dom))] <- unname(dom)
  }
  if (!is.null(conf)) {
    mc_mean <- vapply(split(conf, region_i), function(v) {
      if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
    }, numeric(1))
    S4Vectors::mcols(merged_gr)$mean_confidence[as.integer(names(mc_mean))] <-
      unname(mc_mean)
  }
  merged_gr
}
