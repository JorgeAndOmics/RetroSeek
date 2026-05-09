# -----------------------------------------------------------------------------
# hotspot_analysis / plots.R
# -----------------------------------------------------------------------------
# All hotspot-detector plot families:
#
#   * `plot_manhattan()`           — genome-wide -log10(qval_nb) vs cumulative
#                                    position. One series per chromosome.
#   * `plot_karyotype()`           — chromosome ideograms with hotspots overlaid.
#   * `plot_qq()`                  — observed vs expected -log10(p) under the
#                                    uniform null. Overlays NB and (when
#                                    present) permutation series.
#   * `plot_summary_panel()`       — patchwork composition: hotspots-per-Mb
#                                    bar chart per chromosome + hotspot-width
#                                    histogram.
#   * `plot_perm_histogram()`      — permutation null histogram per label
#                                    (only emitted when validation pass ran).
#   * `plot_perm_density()`        — same data as kernel density.
#   * `save_plots_pdf_pages()`     — multi-page PDF helper.
#   * `save_blank_pdf()`           — emits a 1-page placeholder for the
#                                    Snakemake static-output contract when
#                                    a plot family is conditionally suppressed.
#
# Pure module: pass tibbles / GRanges in, get ggplot objects out.

suppressMessages({
  library(ggplot2)
  library(ggsci)
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


#' Manhattan plot: genome-wide -log10(qval_nb) vs stitched genomic position.
#'
#' Points coloured by chromosome via `ggsci::scale_colour_futurama()`
#' (interpolated when chrom count > palette size). Dashed horizontal line at
#' the significance threshold.
plot_manhattan <- function(window_df, threshold,
                           title = "Hotspot Manhattan plot",
                           subtitle = NULL) {
  if (nrow(window_df) == 0L) {
    return(.empty_plot(title, "No windows to plot"))
  }
  seqlengths <- tapply(window_df$end, window_df$chrom, max)
  offsets <- .chrom_offsets(seqlengths)
  df <- window_df %>%
    dplyr::filter(!is.na(.data$qval_nb)) %>%
    dplyr::left_join(offsets, by = "chrom") %>%
    dplyr::mutate(
      cum_pos = .data$chrom_offset + .data$start,
      neg_log10_q = -log10(pmax(.data$qval_nb, .Machine$double.xmin))
    )
  if (nrow(df) == 0L) {
    return(.empty_plot(title, "No callable windows"))
  }
  pal_n <- length(unique(df$chrom))
  base_pal <- ggsci::pal_futurama()(min(pal_n, 12L))
  colours <- if (pal_n > length(base_pal)) {
    grDevices::colorRampPalette(base_pal)(pal_n)
  } else {
    base_pal
  }
  ggplot(df, aes(x = .data$cum_pos, y = .data$neg_log10_q, colour = .data$chrom)) +
    geom_point(size = 0.6, alpha = 0.7) +
    geom_hline(yintercept = -log10(threshold),
               colour = "red", linetype = "dashed") +
    scale_colour_manual(values = colours, guide = "none") +
    scale_x_continuous(
      breaks = offsets$chrom_centre,
      labels = offsets$chrom,
      expand = expansion(mult = 0.01)
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Chromosome",
      y = expression(-log[10](q[NB]))
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}


#' Karyotype / ideogram plot: chromosomes as horizontal bars with hotspots
#' overlaid. Implemented in pure ggplot2 (faceted) rather than ggbio so the
#' module is robust to ggbio API drift across versions. Each chromosome
#' gets its own facet row; hotspots appear as coloured rectangles.
plot_karyotype <- function(seqlengths, hotspots,
                           title = "Hotspot karyotype",
                           subtitle = NULL) {
  if (length(seqlengths) == 0L) {
    return(.empty_plot(title, "No chromosomes"))
  }
  chrom_df <- tibble::tibble(
    chrom = factor(names(seqlengths), levels = names(seqlengths)),
    start_mb = 0,
    end_mb = as.numeric(seqlengths) / 1e6
  )
  p <- ggplot() +
    geom_rect(
      data = chrom_df,
      aes(xmin = .data$start_mb, xmax = .data$end_mb,
          ymin = -0.4, ymax = 0.4),
      fill = "grey92", colour = "grey60", linewidth = 0.3
    )
  if (length(hotspots) > 0L) {
    hs_df <- tibble::tibble(
      chrom = factor(as.character(GenomicRanges::seqnames(hotspots)),
                     levels = names(seqlengths)),
      start_mb = as.numeric(BiocGenerics::start(hotspots)) / 1e6,
      end_mb   = as.numeric(BiocGenerics::end(hotspots)) / 1e6,
      label    = as.character(S4Vectors::mcols(hotspots)$label),
      count    = as.integer(S4Vectors::mcols(hotspots)$count)
    )
    p <- p +
      geom_rect(
        data = hs_df,
        aes(xmin = .data$start_mb, xmax = .data$end_mb,
            ymin = -0.4, ymax = 0.4, fill = .data$label),
        colour = NA, alpha = 0.85
      ) +
      ggsci::scale_fill_futurama(name = "Label")
  }
  p +
    facet_grid(rows = vars(.data$chrom), switch = "y") +
    scale_y_continuous(breaks = NULL, expand = expansion(mult = 0)) +
    scale_x_continuous(labels = scales::comma_format(suffix = " Mb")) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.spacing.y = unit(0.05, "lines")
    )
}


#' Q-Q plot of -log10(p-values).
#'
#' Plots the NB series (always) and, when `pval_perm` is present, the
#' permutation series too. Diagonal reference. Per the plan agent's
#' caveat: NB and permutation nulls differ by construction; do not
#' interpret tail divergence as model misfit per se.
plot_qq <- function(window_df,
                    title = "Hotspot p-value Q-Q plot",
                    subtitle = NULL) {
  build_qq_series <- function(p, name) {
    p <- p[!is.na(p) & p > 0]
    if (length(p) == 0L) return(NULL)
    tibble::tibble(
      observed = -log10(sort(p)),
      expected = -log10((seq_along(p) - 0.5) / length(p)),
      series   = name
    )
  }
  series_list <- list(build_qq_series(window_df$pval_nb, "NB"))
  if ("pval_perm" %in% colnames(window_df)) {
    series_list <- c(series_list,
                     list(build_qq_series(window_df$pval_perm, "Permutation")))
  }
  series_list <- Filter(Negate(is.null), series_list)
  if (length(series_list) == 0L) {
    return(.empty_plot(title, "No p-values to plot"))
  }
  df <- dplyr::bind_rows(series_list)
  ggplot(df, aes(x = .data$expected, y = .data$observed,
                 colour = .data$series)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey50",
                linetype = "dashed") +
    geom_point(size = 0.7, alpha = 0.7) +
    ggsci::scale_colour_nejm(name = "Method") +
    labs(
      title = title,
      subtitle = subtitle %||% "NB and permutation nulls differ; tail divergence is not necessarily misfit",
      x = expression(Expected~~-log[10](p)),
      y = expression(Observed~~-log[10](p))
    ) +
    theme_minimal()
}


#' Two-panel summary: hotspots-per-Mb per chromosome + width histogram.
#' Composed via `patchwork`.
plot_summary_panel <- function(hotspots, seqlengths,
                               title = "Hotspot summary",
                               subtitle = NULL) {
  if (length(hotspots) == 0L) {
    return(.empty_plot(title, "No hotspots passed thresholds"))
  }
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
    dplyr::arrange(dplyr::desc(.data$hotspots_per_mb)) %>%
    dplyr::mutate(chrom = factor(.data$chrom, levels = .data$chrom))

  p_density <- ggplot(per_chrom_density,
                      aes(x = .data$chrom, y = .data$hotspots_per_mb)) +
    geom_col(fill = ggsci::pal_futurama()(1), alpha = 0.85) +
    labs(
      title = "Hotspots per Mb",
      x = "Chromosome",
      y = "Hotspots / Mb"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

  p_widths <- ggplot(hs_df, aes(x = .data$width_kb)) +
    geom_histogram(bins = 30L, fill = ggsci::pal_futurama()(2)[2], alpha = 0.85) +
    labs(
      title = "Hotspot width distribution",
      x = "Width (kb)",
      y = "Count"
    ) +
    theme_minimal()

  (p_density / p_widths) +
    patchwork::plot_annotation(title = title, subtitle = subtitle)
}


#' Per-label histogram of permutation null counts (validation pass only).
plot_perm_histogram <- function(perm_values, observed_value, label, species,
                                threshold = NULL) {
  if (length(perm_values) == 0L) {
    return(.empty_plot(paste("Permutation histogram —", label),
                       "No permutation values"))
  }
  df <- tibble::tibble(perm_overlaps = as.numeric(perm_values))
  ggplot(df, aes(x = .data$perm_overlaps)) +
    geom_histogram(bins = 30L, fill = "lightblue", colour = "black") +
    geom_vline(xintercept = observed_value,
               colour = "red", linetype = "dashed", linewidth = 1) +
    labs(
      title = paste("Permutation null histogram —", species),
      subtitle = paste0("Label: ", label,
                        " | Observed = ", observed_value),
      x = "Number of overlaps under null",
      y = "Frequency"
    ) +
    theme_minimal()
}


#' Per-label kernel density of permutation null counts (validation pass only).
plot_perm_density <- function(perm_values, observed_value, label, species) {
  if (length(perm_values) == 0L) {
    return(.empty_plot(paste("Permutation density —", label),
                       "No permutation values"))
  }
  df <- tibble::tibble(perm_overlaps = as.numeric(perm_values))
  ggplot(df, aes(x = .data$perm_overlaps)) +
    geom_density(fill = "lightblue", alpha = 0.5) +
    geom_vline(xintercept = observed_value,
               colour = "red", linetype = "dashed", linewidth = 1) +
    labs(
      title = paste("Permutation null density —", species),
      subtitle = paste("Label:", label,
                       "| Observed =", observed_value),
      x = "Number of overlaps under null",
      y = "Density"
    ) +
    theme_minimal()
}


# Internal helper for "I have nothing meaningful to plot" placeholder pages.
.empty_plot <- function(title, message_text) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = message_text, size = 5, colour = "grey40") +
    xlim(0, 1) + ylim(0, 1) +
    labs(title = title) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5))
}


#' Multi-page PDF saver. Each element of `plots` becomes one page.
save_plots_pdf_pages <- function(plots, path, width, height) {
  grDevices::pdf(path, width = width, height = height)
  on.exit(grDevices::dev.off())
  for (p in plots) {
    print(p)
  }
  invisible(path)
}


#' Emit a 1-page blank-but-valid PDF for the Snakemake static-output contract
#' when a plot family is conditionally suppressed (e.g. permutation pass off).
save_blank_pdf <- function(path, title, message_text,
                           width = 8, height = 6) {
  save_plots_pdf_pages(list(.empty_plot(title, message_text)),
                       path = path, width = width, height = height)
}


# Local %||% (rlang's is shadowed when scripts source via `source()` rather
# than `library(rlang)` — we already library(rlang) at the top, but keeping
# this local fallback makes the module self-contained for unit testing).
`%||%` <- function(x, y) if (is.null(x) || (length(x) == 1L && is.na(x))) y else x
