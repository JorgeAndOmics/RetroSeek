# =============================================================================
# plot2sort/plots_sankey.R — alluvial / Sankey builders
# =============================================================================
# All three Sankey wrappers share `.sankey_two_axis`. Each wrapper tags the
# returned ggplot with an `intended_dims` attribute so the orchestrator can
# size the canvas based on the larger of the two axis cardinalities.
#
# `top_n = NULL` (the default) shows every stratum — the existing collapsing
# behaviour is opt-in via `cfg$plots$sankey_top_n`. See helpers.R::collapse_long_tail.


# Generic two-axis Sankey builder shared by the three exposed wrappers below.
# Applies long-tail collapse + count ordering on each axis, then re-aggregates
# in case the "Other" stratum folded multiple rows into the same (axis_a, axis_b)
# combination.
.sankey_two_axis <- function(data, axis_a, axis_b, fill_axis,
                             top_n = NULL, other_label = "Other",
                             title = NULL, subtitle = NULL,
                             subset_label = NULL) {
  if (nrow(data) == 0L) return(empty_plot())

  data <- data %>%
    collapse_long_tail(axis_a, top_n = top_n,
                       other_label = other_label, weight = "count") %>%
    collapse_long_tail(axis_b, top_n = top_n,
                       other_label = other_label, weight = "count") %>%
    dplyr::group_by(.data[[axis_a]], .data[[axis_b]]) %>%
    dplyr::summarise(count = sum(count), .groups = "drop")

  ordered_a <- order_by_count(data, axis_a, weight = "count")
  ordered_b <- order_by_count(data, axis_b, weight = "count")

  data <- data %>%
    dplyr::mutate(
      "{axis_a}" := factor(.data[[axis_a]], levels = ordered_a),
      "{axis_b}" := factor(.data[[axis_b]], levels = ordered_b)
    )

  n_strata       <- max(length(ordered_a), length(ordered_b))
  manual_colours <- futurama_unlimited_palette(12, n_strata)

  p <- ggplot(data, aes(axis1 = .data[[axis_a]],
                        axis2 = .data[[axis_b]],
                        y = count)) +
    geom_alluvium(aes(fill = .data[[fill_axis]], alpha = count)) +
    geom_stratum(alpha = 0, color = "black", linewidth = 0.2) +
    geom_text(
      aes(
        size  = after_stat(count),
        label = sprintf("%s [%d]", after_stat(stratum), after_stat(count))
      ),
      stat     = "stratum",
      nudge_x  = -0.15,
      fontface = "bold",
      hjust    = 0
    ) +
    scale_size_continuous(range = c(2, 5)) +
    scale_alpha_continuous(range = c(0.4, 1)) +
    scale_fill_manual(values = manual_colours) +
    theme_void() +
    theme(legend.position = "none")
  out <- if (!is.null(title)) {
    add_titles(p, title = title, subtitle = subtitle,
               subset_label = subset_label)
  } else {
    p
  }
  # Sankey strata stack vertically; the canvas grows along height.
  attr(out, "intended_dims") <- auto_dims(n_strata, axis = "y")
  out
}


sankey_species_probe_plot <- function(data, top_n = NULL,
                                      other_label = "Other",
                                      subset_label = NULL) {
  .sankey_two_axis(data, "species", "probe", fill_axis = "species",
                   top_n = top_n, other_label = other_label,
                   title    = "Flow: species → probe",
                   subtitle = "Range-count alluvium between species and probe",
                   subset_label = subset_label)
}


sankey_label_probe_plot <- function(data, top_n = NULL,
                                    other_label = "Other",
                                    subset_label = NULL) {
  .sankey_two_axis(data, "label", "probe", fill_axis = "label",
                   top_n = top_n, other_label = other_label,
                   title    = "Flow: label → probe",
                   subtitle = "Range-count alluvium between category label and probe",
                   subset_label = subset_label)
}


sankey_species_label_plot <- function(data, top_n = NULL,
                                      other_label = "Other",
                                      subset_label = NULL) {
  .sankey_two_axis(data, "species", "label", fill_axis = "species",
                   top_n = top_n, other_label = other_label,
                   title    = "Flow: species → label",
                   subtitle = "Range-count alluvium between species and category label",
                   subset_label = subset_label)
}
