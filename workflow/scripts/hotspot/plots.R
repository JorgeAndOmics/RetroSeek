# -----------------------------------------------------------------------------
# hotspot / plots.R
# -----------------------------------------------------------------------------
# The hotspot detector's pages (NB GLM only), in the house style
# (plot2sort/style.R, docs/visual_style.md; hotspot_detector.R sources style.R
# and helpers.R before this file):
#
#   * `plot_manhattan()`           - genome-wide -log10(qval_nb) against stitched
#                                    position, chromosomes in alternating shades.
#   * `plot_karyotype()`           - chromosome ideograms with hotspots overlaid,
#                                    each in its lineage colour.
#   * `plot_qq()`                  - observed against expected -log10(p) under the
#                                    uniform null (NB series).
#   * `plot_summary_panel()`       - hotspots per Mb per chromosome, and hotspot
#                                    widths.
#   * `plot_hotspot_composition()` - what each hotspot is made of, by structural
#                                    class.
#
# `species` is the readable species name; it leads each page's subtitle. The
# detector writes every page into one PDF per genome.
#
# Pure module: pass tibbles / GRanges in, get ggplot objects out.

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(scales)
  library(patchwork)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(rlang)
})


# Build a chromosome offset table giving the cumulative bp at which each
# chromosome's local coordinates begin in a stitched genome-wide x-axis.
# Used by Manhattan and per-chrom density bar.
.chrom_offsets <- function(seqlengths) {
  ord <- order(names(seqlengths))
  sl_ord <- seqlengths[ord]
  offsets <- c(0L, cumsum(as.numeric(sl_ord)))[seq_along(sl_ord)]
  tibble::tibble(
    chrom        = names(sl_ord),
    chrom_length = as.numeric(sl_ord),
    chrom_offset = as.numeric(offsets),
    chrom_centre = as.numeric(offsets) + as.numeric(sl_ord) / 2
  )
}


# Smaller titles for the panels inside a composed page, under the page title.
.panel_title <- theme(plot.title = element_text(size = 12, face = "bold"))


#' Manhattan plot: genome-wide -log10(qval_nb) against stitched genomic position,
#' chromosomes in alternating shades so neighbours stay apart. Dashed line at
#' the significance threshold.
plot_manhattan <- function(window_df, threshold, species = NULL, label = NULL) {
  if (nrow(window_df) == 0L) return(empty_plot("No windows to plot"))
  seqlengths <- tapply(window_df$end, window_df$chrom, max)
  offsets <- .chrom_offsets(seqlengths)
  df <- window_df %>%
    dplyr::filter(!is.na(.data$qval_nb)) %>%
    dplyr::left_join(offsets, by = "chrom") %>%
    dplyr::mutate(
      cum_pos = .data$chrom_offset + .data$start,
      neg_log10_q = -log10(pmax(.data$qval_nb, .Machine$double.xmin)),
      band = factor(match(.data$chrom, offsets$chrom) %% 2L)
    )
  if (nrow(df) == 0L) return(empty_plot("No callable windows"))
  p <- ggplot(df, aes(x = .data$cum_pos, y = .data$neg_log10_q, colour = .data$band)) +
    geom_point(size = 0.6, alpha = 0.8) +
    geom_hline(yintercept = -log10(threshold), colour = .INK_SOFT, linetype = "dashed") +
    scale_colour_manual(values = c(`0` = .DATA_COLOUR, `1` = .GREY_MID), guide = "none") +
    scale_x_continuous(breaks = offsets$chrom_centre, labels = offsets$chrom,
                       expand = expansion(mult = 0.01)) +
    labs(x = NULL, y = expression(-log[10](q))) +
    # Sequence names are a dense axis of accession codes: the one place a tilted
    # label is the lesser evil.
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          panel.grid.major.x = element_blank())
  add_titles(p, "Where integrations cluster",
             sprintf(paste("%s loci per window, tested against a negative binomial model.",
                           "Dashed line: q = %s."),
                     label %||% "All", format(threshold)),
             subset_label = species)
}


#' Karyotype / ideogram plot: chromosomes as horizontal bars with hotspots
#' overlaid in their lineage colour. Pure ggplot2 (faceted) rather than ggbio,
#' so the module is robust to ggbio API drift across versions.
plot_karyotype <- function(seqlengths, hotspots, species = NULL) {
  if (length(seqlengths) == 0L) return(empty_plot("No chromosomes"))
  chrom_df <- tibble::tibble(
    chrom = factor(names(seqlengths), levels = names(seqlengths)),
    start_mb = 0,
    end_mb = as.numeric(seqlengths) / 1e6
  )
  p <- ggplot() +
    geom_rect(data = chrom_df,
              aes(xmin = .data$start_mb, xmax = .data$end_mb, ymin = -0.4, ymax = 0.4),
              fill = .GREY_OTHER, colour = .GREY_MID, linewidth = 0.3)
  if (length(hotspots) > 0L) {
    hs_df <- tibble::tibble(
      chrom = factor(as.character(GenomicRanges::seqnames(hotspots)),
                     levels = names(seqlengths)),
      start_mb = as.numeric(BiocGenerics::start(hotspots)) / 1e6,
      end_mb   = as.numeric(BiocGenerics::end(hotspots)) / 1e6,
      label    = as.character(S4Vectors::mcols(hotspots)$label)
    )
    p <- p +
      geom_rect(data = hs_df,
                aes(xmin = .data$start_mb, xmax = .data$end_mb,
                    ymin = -0.4, ymax = 0.4, fill = .data$label),
                colour = NA) +
      scale_fill_manual(values = taxon_colours(hs_df$label), labels = taxon_labels,
                        name = NULL)
  }
  p <- p +
    facet_grid(rows = vars(.data$chrom), switch = "y") +
    scale_y_continuous(breaks = NULL, expand = expansion(mult = 0)) +
    scale_x_continuous(labels = scales::comma_format(suffix = " Mb")) +
    labs(x = NULL, y = NULL) +
    theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7, face = "plain"),
          panel.grid.major.y = element_blank(),
          panel.spacing.y = unit(0.05, "lines"))
  add_titles(p, "Hotspots along the chromosomes",
             sprintf("%d hotspots across %d sequences.", length(hotspots), length(seqlengths)),
             subset_label = species)
}


#' Q-Q plot of -log10(p-values) for the NB series against the uniform null.
plot_qq <- function(window_df, species = NULL, label = NULL) {
  p <- window_df$pval_nb
  p <- p[!is.na(p) & p > 0]
  if (length(p) == 0L) return(empty_plot("No p-values to plot"))
  df <- tibble::tibble(
    observed = -log10(sort(p)),
    expected = -log10((seq_along(p) - 0.5) / length(p))
  )
  plot <- ggplot(df, aes(x = .data$expected, y = .data$observed)) +
    geom_abline(slope = 1, intercept = 0, colour = .INK_SOFT, linetype = "dashed") +
    geom_point(size = 0.7, alpha = 0.7, colour = .DATA_COLOUR) +
    labs(x = expression(Expected~~-log[10](p)), y = expression(Observed~~-log[10](p)))
  add_titles(plot, "Is the model calibrated?",
             sprintf(paste("%s windows: observed p-values against the uniform null. Points",
                           "on the dashed line mean a calibrated model."),
                     label %||% "All"),
             subset_label = species)
}


#' Two panels: hotspots per Mb per chromosome, and hotspot widths.
plot_summary_panel <- function(hotspots, seqlengths, species = NULL) {
  if (length(hotspots) == 0L) return(empty_plot("No hotspot passed the threshold"))
  hs_df <- tibble::tibble(
    chrom = as.character(GenomicRanges::seqnames(hotspots)),
    width_kb = as.numeric(BiocGenerics::width(hotspots)) / 1e3
  )
  per_chrom_density <- hs_df %>%
    dplyr::count(.data$chrom, name = "n") %>%
    dplyr::mutate(
      chrom_length_mb = as.numeric(seqlengths[.data$chrom]) / 1e6,
      hotspots_per_mb = .data$n / .data$chrom_length_mb
    ) %>%
    dplyr::arrange(.data$hotspots_per_mb) %>%
    dplyr::mutate(chrom = factor(.data$chrom, levels = .data$chrom))

  p_density <- ggplot(per_chrom_density,
                      aes(x = .data$chrom, y = .data$hotspots_per_mb)) +
    geom_col(fill = .DATA_COLOUR, width = 0.7) +
    coord_flip() +
    labs(title = "Hotspots per Mb", x = NULL, y = "Hotspots per Mb") +
    theme(panel.grid.major.y = element_blank()) + .panel_title

  p_widths <- ggplot(hs_df, aes(x = .data$width_kb)) +
    geom_histogram(bins = 30L, fill = .DATA_COLOUR) +
    labs(title = "Hotspot widths", x = "Width (kb)", y = "Hotspots") + .panel_title

  (p_density / p_widths) +
    patchwork::plot_annotation(
      title = "Hotspot summary",
      subtitle = paste(c(species, sprintf("%d hotspots.", length(hotspots))), collapse = ". "),
      theme = theme_retroseek())
}


#' Per-hotspot structural composition (ADR-012).
#'
#' Detection counts every locus in the tier, so a called region can be driven by
#' intact proviruses or by single-gene fragments and the p-value looks the same.
#' This page makes that visible: hotspots as rows, largest on top, by structural
#' class as counts and as shares. It is why composition is annotated rather than
#' filtered on.
#'
#' Empty-safe, and degrades to a placeholder when the input carried no
#' structure_class (a raw-hit run).
plot_hotspot_composition <- function(hotspots, species = NULL, tier_note = NULL) {
  if (length(hotspots) == 0L) return(empty_plot("No hotspot passed the threshold"))
  mc <- S4Vectors::mcols(hotspots)
  if (!all(c("n_full", "n_partial", "n_gene") %in% colnames(mc)) ||
      sum(mc$n_full, mc$n_partial, mc$n_gene, na.rm = TRUE) == 0L) {
    return(empty_plot("No structural composition on this input tier"))
  }
  id <- if ("hotspot_id" %in% colnames(mc)) {
    as.character(mc$hotspot_id)
  } else {
    sprintf("%s:%d", as.character(GenomicRanges::seqnames(hotspots)),
            BiocGenerics::start(hotspots))
  }
  comp <- tibble::tibble(
    hotspot = rep(id, 3L),
    class = rep(c("full", "partial", "gene"), each = length(id)),
    n = c(as.integer(mc$n_full), as.integer(mc$n_partial), as.integer(mc$n_gene))
  ) %>%
    dplyr::mutate(
      # Largest on top once the axis is flipped.
      hotspot = factor(.data$hotspot, levels = id[order(as.numeric(mc$n_loci))]),
      class = factor(.data$class, levels = c("full", "partial", "gene"))
    )
  bars <- function(position) {
    ggplot(comp, aes(x = .data$hotspot, y = .data$n, fill = .data$class)) +
      geom_col(position = position, width = 0.7) +
      coord_flip() +
      scale_fill_manual(values = .STRUCTURE_COLOUR, labels = display_label, name = NULL) +
      labs(x = NULL) +
      theme(panel.grid.major.y = element_blank()) + .panel_title
  }
  p_counts <- bars(position_stack(reverse = TRUE)) +
    labs(title = "Loci per hotspot", y = "Loci")
  p_share <- bars(position_fill(reverse = TRUE)) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "The same, as shares", y = "Share of the hotspot's loci")

  (p_counts / p_share) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = "What each hotspot is made of",
      subtitle = paste(c(species, tier_note %||% "Loci by structural class."), collapse = ". "),
      theme = theme_retroseek())
}


# Local %||% (rlang's is shadowed when scripts source via `source()` rather
# than `library(rlang)` - we already library(rlang) at the top, but keeping
# this local fallback makes the module self-contained for unit testing).
`%||%` <- function(x, y) if (is.null(x) || (length(x) == 1L && is.na(x))) y else x
