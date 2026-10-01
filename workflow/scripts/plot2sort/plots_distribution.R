# =============================================================================
# plot2sort/plots_distribution.R - distributions over a continuous measure
# =============================================================================
# Probe-keyed distributions of bitscore and query coverage. The x axis is the
# measure, so none of these grows with the species count. `colours` is the
# stage's probe palette (probe_colours() over every probe), so a probe keeps its
# colour on the main and the accessory pages alike.


# Density of the strongest HSP bitscore per merged range, by probe, with the
# quartiles marked. Overlaid curves stay readable for a few probes; past
# MAX_OVERLAID_PROBES each probe gets its own panel (the accessory set has ten).
# The quartile values go in the subtitle: close quartiles printed as labels on
# the plot landed on top of each other.
MAX_OVERLAID_PROBES <- 5L

# The bitscore axis: log10 on request; otherwise a tick per 100 bitscore on one
# panel, about four on a small panel, which has no room for more.
.bitscore_x_scale <- function(x_scale, many, max_bs) {
  if (identical(x_scale, "log10")) return(scale_x_log10())
  scale_x_continuous(breaks = if (many) scales::breaks_extended(n = 4) else
                       seq(0, max_bs, by = 100))
}

density_bitscore_plot <- function(data, q1, median, q3, x_scale = "linear",
                                  subset_label = NULL, colours = NULL) {
  if (nrow(data) == 0L) return(empty_plot())
  if (is.null(colours)) colours <- probe_colours(data$probe)

  # Pin max_bitscore eagerly so ggplot2's later scale build doesn't force it
  # in a deferred-evaluation context where it can mis-resolve.
  max_bs <- max(data$max_bitscore, na.rm = TRUE)
  many <- length(unique(data$probe)) > MAX_OVERLAID_PROBES

  # Unweighted density over max_bitscore - see commit history for the prior
  # weighted version. If per-range identity needs visual emphasis, prefer a
  # separate plot rather than a weight aesthetic.
  p <- ggplot(data, aes(x = max_bitscore, fill = probe, colour = probe)) +
    geom_density(alpha = 0.35, adjust = 3) +
    geom_vline(xintercept = c(q1, median, q3), linetype = "dashed",
               linewidth = 0.3, colour = .INK_SOFT) +
    scale_fill_manual(values = colours, guide = if (many) "none" else "legend") +
    scale_colour_manual(values = colours, guide = if (many) "none" else "legend") +
    labs(x = "Strongest HSP bitscore per range", y = "Density",
         fill = "Probe", colour = "Probe")
  if (many) p <- p + facet_wrap(~probe, ncol = 2, scales = "free_y")
  p <- p + .bitscore_x_scale(x_scale, many, max_bs)
  add_titles(p,
             title    = "Bitscore density",
             subtitle = sprintf(paste(
               "The strongest HSP bitscore of each merged range, by probe.",
               "Dashed lines: the quartiles (%s, %s and %s)."
             ), format(round(q1)), format(round(median)), format(round(q3))),
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
