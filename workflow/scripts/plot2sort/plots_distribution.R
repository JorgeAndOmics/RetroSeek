# =============================================================================
# plot2sort/plots_distribution.R — density-family builders
# =============================================================================
# Probe-keyed distributions over a continuous variable. None of these scale
# with species count — the x-axis is bitscore (density / raincloud) or
# query_coverage in [0, 1], and the fill is probe (~3 levels). Fixed canvas.


density_bitscore_plot <- function(data, q1, median, q3, x_scale = "linear",
                                  subset_label = NULL) {
  if (nrow(data) == 0L) return(empty_plot())
  manual_colours <- futurama_unlimited_palette(3, length(unique(data$probe)))

  # Pin max_bitscore eagerly so ggplot2's later scale build doesn't force it
  # in a deferred-evaluation context where it can mis-resolve.
  max_bs <- max(data$max_bitscore, na.rm = TRUE)

  # Unweighted density over max_bitscore — see commit history for the prior
  # weighted version. If per-range identity needs visual emphasis, prefer a
  # separate plot rather than a weight aesthetic.
  p <- ggplot(data, aes(x = max_bitscore, fill = probe, colour = probe)) +
    geom_density(alpha = 0.4, adjust = 3) +
    scale_fill_manual(values = manual_colours) +
    scale_color_manual(values = manual_colours) +
    geom_vline(xintercept = c(q1, median, q3),
               linetype = "dashed", linewidth = 0.2) +
    annotate("text", x = q1 + 20,     y = 1e-4, label = "Q1") +
    annotate("text", x = median + 20, y = 1e-4, label = "Q2") +
    annotate("text", x = q3 + 20,     y = 1e-4, label = "Q3") +
    labs(x = "HSP Bitscore (max per range)", y = "Density") +
    theme_minimal() +
    theme(
      text         = element_text(face = "bold"),
      axis.title.x = element_text(size = 15, margin = margin(b = 15)),
      axis.title.y = element_text(size = 15, margin = margin(l = 10)),
      axis.text    = element_text(size = 12)
    )
  p <- if (identical(x_scale, "log10")) {
    p + scale_x_log10()
  } else {
    p + scale_x_continuous(breaks = seq(0, max_bs, by = 100))
  }
  add_titles(p,
             title    = "Bitscore density",
             subtitle = paste0("Distribution of the strongest HSP bitscore per merged range, by probe ",
                               "(dashed lines: Q1 / median / Q3)"),
             subset_label = subset_label)
}


raincloud_bitscore_plot <- function(data, x_scale = "linear",
                                    subset_label = NULL) {
  if (nrow(data) == 0L) return(empty_plot())
  manual_colours <- futurama_unlimited_palette(3, length(unique(data$probe)))

  p <- ggplot(data, aes(y = probe, x = max_bitscore,
                        fill = probe, color = probe)) +
    ggdist::stat_halfeye(adjust = 0.5, justification = 0,
                         alpha = 0.5, .width = c(0.5, 0.95)) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
    ggdist::stat_dots(side = "right", dotsize = 0.1,
                      alpha = 0.01, binwidth = 0.2) +
    scale_fill_manual(values = manual_colours) +
    scale_color_manual(values = manual_colours) +
    theme_minimal() +
    labs(x = "HSP Bitscore (max per range)", y = "Probe") +
    theme(
      text       = element_text(face = "bold"),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 11)
    )
  p <- if (identical(x_scale, "log10")) p + scale_x_log10() else p
  add_titles(p,
             title    = "Bitscore distribution per probe",
             subtitle = "Halfeye + boxplot + dot strip of the max HSP bitscore per merged range",
             subset_label = subset_label)
}


# Distribution of `query_coverage` per probe — reveals under-aligned probes
# whose hits cover only a fraction of the probe sequence (candidates for pHMM
# follow-up). Input: per-range tibble (one row = one merged range).
query_coverage_plot <- function(data, subset_label = NULL) {
  if (nrow(data) == 0L) return(empty_plot())
  if (!"query_coverage" %in% colnames(data)) {
    return(empty_plot("query_coverage column missing"))
  }
  if (all(is.na(data$query_coverage))) {
    return(empty_plot("query_coverage all NA"))
  }
  manual_colours <- futurama_unlimited_palette(3, length(unique(data$probe)))

  p <- ggplot(data, aes(x = query_coverage, fill = probe, colour = probe)) +
    geom_density(alpha = 0.4, adjust = 1.5) +
    scale_fill_manual(values = manual_colours) +
    scale_color_manual(values = manual_colours) +
    labs(x = "Query coverage (alignment length / probe length)",
         y = "Density") +
    theme_minimal() +
    theme(
      text         = element_text(face = "bold"),
      axis.title.x = element_text(size = 15, margin = margin(b = 15)),
      axis.title.y = element_text(size = 15, margin = margin(l = 10)),
      axis.text    = element_text(size = 12)
    )
  add_titles(p,
             title    = "Probe query coverage density",
             subtitle = paste0("Per-range alignment length / probe length, by probe ",
                               "(1.0 = full-length probe coverage)"),
             subset_label = subset_label)
}
