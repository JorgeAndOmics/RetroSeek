# =============================================================================
# stage_plot_generator/plots_overlap.R — pre-reduction overlap / redundancy
# =============================================================================
# How much the UNREDUCED original tier (gr_virus) overlaps itself, and how much
# that redundancy collapses under reduction. These plots make the reduction
# step's effect visible. Same builder contract as the other plot modules —
# `(data, ..., subset_label, warning_caption)` → ggplot. Reuses plot2sort
# helpers. All four describe NON-reduced inputs except the explicit before/after
# comparisons; the orchestrator stamps the reduced-state note.

.REDUCE_FILL <- c(unreduced = "#7570b3", reduced = "#1b9e77")


# Per-probe stacked bars of the self-overlap degree — how many OTHER unreduced
# loci each locus overlaps. A tall right tail means heavy positional redundancy
# that reduction will collapse.
overlap_degree_plot <- function(overlap_df, subset_label = NULL,
                                warning_caption = NULL) {
  if (nrow(overlap_df) == 0L) return(empty_plot())
  d <- overlap_df %>% dplyr::count(probe, overlap_degree, name = "loci")
  ordered_probe <- order_by_count(d, "probe", weight = "loci")
  d <- d %>% dplyr::mutate(probe = factor(probe, levels = ordered_probe))

  p <- ggplot(d, aes(x = overlap_degree, y = loci, fill = probe)) +
    geom_col(colour = NA) +
    scale_x_continuous(breaks = scales::breaks_pretty()) +
    scale_y_continuous(labels = scales::label_comma()) +
    scale_fill_manual(values = futurama_unlimited_palette(12, length(ordered_probe))) +
    theme_minimal() +
    labs(x = "Overlapping unreduced loci (degree)", y = "Loci", fill = "Probe") +
    theme(text = element_text(face = "bold"))
  add_titles(
    p,
    title    = "Pre-reduction self-overlap degree",
    subtitle = "Other unreduced loci each locus overlaps — the redundancy reduction collapses",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}


# Histogram of the maximum reciprocal-overlap fraction per locus (over its
# overlapping partners), restricted to loci that overlap something. 1.0 = a
# locus fully contained in / containing a neighbour.
reciprocal_fraction_plot <- function(overlap_df, subset_label = NULL,
                                     warning_caption = NULL) {
  if (nrow(overlap_df) == 0L) return(empty_plot())
  d <- overlap_df %>% dplyr::filter(overlap_degree > 0L)
  if (nrow(d) == 0L) return(empty_plot("no overlapping loci"))

  p <- ggplot(d, aes(x = max_reciprocal_fraction)) +
    geom_histogram(bins = 30, fill = "#386cb0", colour = NA, alpha = 0.85) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(labels = scales::label_comma()) +
    theme_minimal() +
    labs(x = "Max reciprocal overlap fraction (intersection / min width)",
         y = "Loci") +
    theme(text = element_text(face = "bold"))
  add_titles(
    p,
    title    = "Pre-reduction reciprocal-overlap fraction",
    subtitle = "Among overlapping unreduced loci; 1.0 = full containment",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}


# Per-probe dodged bars of locus count before (unreduced gr_virus) vs after
# (globally reduced gr_global) reduction — the collapse fold per probe.
reduction_fold_plot <- function(hits_df, reduced_df, subset_label = NULL,
                                warning_caption = NULL) {
  if (nrow(hits_df) == 0L && nrow(reduced_df) == 0L) return(empty_plot())
  unred <- if (nrow(hits_df) > 0L) {
    hits_df %>% dplyr::count(probe, name = "count") %>%
      dplyr::mutate(stage = "unreduced")
  } else NULL
  red <- if (nrow(reduced_df) > 0L) {
    reduced_df %>% dplyr::count(probe, name = "count") %>%
      dplyr::mutate(stage = "reduced")
  } else NULL
  d <- dplyr::bind_rows(unred, red) %>%
    dplyr::mutate(stage = factor(stage, levels = c("unreduced", "reduced")))
  ordered_probe <- order_by_count(d, "probe", weight = "count")
  d <- d %>% dplyr::mutate(probe = factor(probe, levels = ordered_probe))

  p <- ggplot(d, aes(x = probe, y = count, fill = stage)) +
    geom_col(position = position_dodge(width = 0.8), colour = "black",
             linewidth = 0.2) +
    scale_fill_manual(values = .REDUCE_FILL) +
    scale_y_continuous(labels = scales::label_comma()) +
    theme_minimal() +
    labs(x = "Probe", y = "Loci", fill = NULL) +
    theme(text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  out <- add_titles(
    p,
    title    = "Reduction collapse fold per probe",
    subtitle = "Loci before (per-virus, unreduced) vs after (per-probe, globally reduced)",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
  attr(out, "intended_dims") <- auto_dims(length(ordered_probe), axis = "x")
  out
}


# Two bars: total range length (bp) summed across loci, unreduced vs reduced.
# The drop is the overlap-redundant length removed by global reduction.
coverage_before_after_plot <- function(coverage_df, subset_label = NULL,
                                       warning_caption = NULL) {
  if (nrow(coverage_df) == 0L) return(empty_plot())
  labels <- c(total_bp_unreduced = "unreduced", total_bp_reduced = "reduced")
  d <- coverage_df %>%
    dplyr::mutate(stage = factor(unname(labels[metric]),
                                 levels = c("unreduced", "reduced")))

  p <- ggplot(d, aes(x = stage, y = value, fill = stage)) +
    geom_col(colour = "black", linewidth = 0.2, width = 0.6) +
    scale_fill_manual(values = .REDUCE_FILL, guide = "none") +
    scale_y_continuous(labels = scales::label_comma()) +
    theme_minimal() +
    labs(x = NULL, y = "Total range length (bp)") +
    theme(text = element_text(face = "bold"))
  add_titles(
    p,
    title    = "Total range length before vs after reduction",
    subtitle = "Summed locus widths; the drop is overlap-redundant length collapsed",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}
