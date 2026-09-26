# =============================================================================
# plot2sort/plots_distribution.R - distributions over a continuous measure
# =============================================================================
# Probe-keyed distributions of bitscore and query coverage. The x axis is the
# measure, so none of these grows with the species count. `colours` is the
# stage's probe palette (probe_colours() over every probe), so a probe keeps its
# colour on the main and the accessory pages alike.


# Density of the strongest HSP bitscore per merged range, by probe, with the
# quartiles marked.
density_bitscore_plot <- function(data, q1, median, q3, x_scale = "linear",
                                  subset_label = NULL, colours = NULL) {
  if (nrow(data) == 0L) return(empty_plot())
  if (is.null(colours)) colours <- probe_colours(data$probe)

  # Pin max_bitscore eagerly so ggplot2's later scale build doesn't force it
  # in a deferred-evaluation context where it can mis-resolve.
  max_bs <- max(data$max_bitscore, na.rm = TRUE)
  quartiles <- tibble::tibble(x = c(q1, median, q3), label = c("Q1", "Median", "Q3"))

  # Unweighted density over max_bitscore - see commit history for the prior
  # weighted version. If per-range identity needs visual emphasis, prefer a
  # separate plot rather than a weight aesthetic.
  p <- ggplot(data, aes(x = max_bitscore, fill = probe, colour = probe)) +
    geom_density(alpha = 0.35, adjust = 3) +
    geom_vline(data = quartiles, aes(xintercept = .data$x),
               linetype = "dashed", linewidth = 0.3, colour = .INK_SOFT) +
    geom_text(data = quartiles, aes(x = .data$x, y = Inf, label = .data$label),
              inherit.aes = FALSE, vjust = 1.5, hjust = -0.1, size = 3.2,
              colour = .INK_SOFT, family = .FONT) +
    scale_fill_manual(values = colours) +
    scale_colour_manual(values = colours) +
    labs(x = "Strongest HSP bitscore per range", y = "Density",
         fill = "Probe", colour = "Probe")
  p <- if (identical(x_scale, "log10")) {
    p + scale_x_log10()
  } else {
    p + scale_x_continuous(breaks = seq(0, max_bs, by = 100))
  }
  add_titles(p,
             title    = "Bitscore density",
             subtitle = paste(
               "The strongest HSP bitscore of each merged range, by probe.",
               "Dashed lines: the quartiles."
             ),
             subset_label = subset_label)
}


# Bitscore per probe as a raincloud: a half density, a box and the ranges as dots.
raincloud_bitscore_plot <- function(data, x_scale = "linear",
                                    subset_label = NULL, colours = NULL) {
  if (nrow(data) == 0L) return(empty_plot())
  if (is.null(colours)) colours <- probe_colours(data$probe)

  p <- ggplot(data, aes(y = probe, x = max_bitscore, fill = probe, colour = probe)) +
    ggdist::stat_halfeye(adjust = 0.5, justification = 0,
                         alpha = 0.5, .width = c(0.5, 0.95)) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7, colour = .INK_SOFT) +
    ggdist::stat_dots(side = "right", dotsize = 0.1,
                      alpha = 0.01, binwidth = 0.2) +
    scale_fill_manual(values = colours, guide = "none") +
    scale_colour_manual(values = colours, guide = "none") +
    labs(x = "Strongest HSP bitscore per range", y = NULL) +
    theme(axis.text.y = element_text(face = "italic"),
          panel.grid.major.y = element_blank())
  p <- if (identical(x_scale, "log10")) p + scale_x_log10() else p
  add_titles(
    p,
    title    = "Bitscore per probe",
    subtitle =
      "The strongest HSP bitscore of each merged range: density, box and ranges.",
    subset_label = subset_label
  )
}


# Distribution of `query_coverage` per probe - reveals under-aligned probes
# whose hits cover only a fraction of the probe sequence (candidates for pHMM
# follow-up). Input: per-range tibble (one row = one merged range).
query_coverage_plot <- function(data, subset_label = NULL, colours = NULL) {
  if (nrow(data) == 0L) return(empty_plot())
  if (!"query_coverage" %in% colnames(data)) {
    return(empty_plot("query_coverage column missing"))
  }
  if (all(is.na(data$query_coverage))) {
    return(empty_plot("query_coverage all NA"))
  }
  if (is.null(colours)) colours <- probe_colours(data$probe)

  p <- ggplot(data, aes(x = query_coverage, fill = probe, colour = probe)) +
    geom_density(alpha = 0.35, adjust = 1.5) +
    scale_fill_manual(values = colours) +
    scale_colour_manual(values = colours) +
    scale_x_continuous(labels = scales::percent) +
    labs(x = "Share of the probe covered by the alignment", y = "Density",
         fill = "Probe", colour = "Probe")
  add_titles(p,
             title    = "How much of each probe the hits cover",
             subtitle = paste(
               "Alignment length over probe length, per range and probe.",
               "100% is a full-length match."
             ),
             subset_label = subset_label)
}
