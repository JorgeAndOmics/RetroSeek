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


#' Count, per region, the overlaps selected by `mask` (all of them when TRUE).
#' A NULL mask means the column is absent: every region counts 0.
.tally_regions <- function(region_i, n, mask) {
  out <- integer(n)
  if (is.null(mask)) return(out)
  t <- table(region_i[mask])
  out[as.integer(names(t))] <- as.integer(t)
  out
}

#' The most frequent non-blank value among one region's loci; NA when none.
#' sort() breaks count ties in the table's level order, so the output is
#' deterministic for a given locale.
.dominant_value <- function(v) {
  v <- v[!is.na(v) & nzchar(v)]
  if (length(v) == 0L) return(NA_character_)
  names(sort(table(v), decreasing = TRUE))[1]
}

#' Mean of one region's numeric values, NA when all are missing.
.mean_or_na <- function(v) {
  if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
}

#' Fill the per-region summary `column` of `gr` from `values` split by region.
.fill_by_region <- function(gr, column, values, region_i, reduce, type) {
  if (is.null(values)) return(gr)
  per_region <- vapply(split(values, region_i), reduce, type)
  S4Vectors::mcols(gr)[[column]][as.integer(names(per_region))] <- unname(per_region)
  gr
}

#' The zero / NA composition columns every region starts with.
.empty_composition <- function(merged_gr) {
  n <- length(merged_gr)
  for (column in c("n_loci", "n_full", "n_partial", "n_gene",
                   "n_ltr_flanked", "n_orphan")) {
    S4Vectors::mcols(merged_gr)[[column]] <- rep(0L, n)
  }
  S4Vectors::mcols(merged_gr)$dominant_taxon  <- rep(NA_character_, n)
  S4Vectors::mcols(merged_gr)$mean_confidence <- rep(NA_real_, n)
  merged_gr
}

#' Which locus overlaps which region, as two parallel integer vectors; NULL when
#' there are no regions, no loci, or no overlap.
.overlap_index <- function(loci, regions) {
  if (length(regions) == 0L || length(loci) == 0L) return(NULL)
  ov <- GenomicRanges::findOverlaps(loci, regions, ignore.strand = TRUE)
  if (length(ov) == 0L) return(NULL)
  list(locus = as.integer(S4Vectors::queryHits(ov)),
       region = as.integer(S4Vectors::subjectHits(ov)))
}

#' `values == level`, or NULL when the column is absent (so it counts zero).
.equal_or_null <- function(values, level) {
  if (is.null(values)) NULL else values == level
}

#' Numeric confidence at the overlapping loci; NULL when the column is absent.
#' Orphan loci carry a blank confidence, which reads as NA.
.locus_confidence <- function(mc, locus_i) {
  if (!"confidence" %in% colnames(mc)) return(NULL)
  suppressWarnings(as.numeric(mc$confidence[locus_i]))
}

#' A loci column at the overlapping loci, as character; NULL when absent.
.locus_column <- function(mc, name, locus_i) {
  if (name %in% colnames(mc)) as.character(mc[[name]][locus_i]) else NULL
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
  merged_gr <- .empty_composition(merged_gr)
  n <- length(merged_gr)
  hits <- .overlap_index(loci, merged_gr)
  if (is.null(hits)) return(merged_gr)
  locus_i  <- hits$locus
  region_i <- hits$region

  mc <- S4Vectors::mcols(loci)
  structure_class <- .locus_column(mc, "structure_class", locus_i)
  tier <- .locus_column(mc, "source", locus_i)
  counts <- list(
    n_loci        = rep(TRUE, length(region_i)),
    n_full        = .equal_or_null(structure_class, "full"),
    n_partial     = .equal_or_null(structure_class, "partial"),
    n_gene        = .equal_or_null(structure_class, "gene"),
    n_ltr_flanked = .equal_or_null(tier, "ltr-flanked"),
    n_orphan      = .equal_or_null(tier, "orphan")
  )
  for (column in names(counts)) {
    S4Vectors::mcols(merged_gr)[[column]] <-
      .tally_regions(region_i, n, counts[[column]])
  }
  conf <- .locus_confidence(mc, locus_i)
  merged_gr <- .fill_by_region(merged_gr, "dominant_taxon",
                               .locus_column(mc, group_col, locus_i), region_i,
                               .dominant_value, character(1))
  .fill_by_region(merged_gr, "mean_confidence", conf, region_i,
                  .mean_or_na, numeric(1))
}
