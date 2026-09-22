# =============================================================================
# plot2sort/plots_sankey.R - alluvial pages
# =============================================================================
# All three alluvial wrappers share `.sankey_two_axis`. Flows are coloured by
# the axis that carries a meaning (the probe, or the lineage), never by host:
# hosts have no colour of their own. A host axis follows the canonical order
# (ctx, see tree_axis.R), top first; other axes are ordered by count.
#
# `top_n = NULL` (the default) shows every stratum - the existing collapsing
# behaviour is opt-in via `cfg$plots$sankey_top_n`. See helpers.R::collapse_long_tail.


# Generic two-axis alluvium shared by the three wrappers below. Applies the
# long-tail collapse on each axis, then re-aggregates in case the "Other"
# stratum folded several rows into the same (axis_a, axis_b) combination.
#   colours      named colours for the fill axis's values ("Other" is added grey)
#   axis_titles  the words printed under the two axes
.sankey_two_axis <- function(data, axis_a, axis_b, fill_axis, colours,
                             axis_titles, top_n = NULL, other_label = "Other",
                             title = NULL, subtitle = NULL,
                             subset_label = NULL, ctx = NULL) {
  if (nrow(data) == 0L) return(empty_plot())

  data <- data %>%
    collapse_long_tail(axis_a, top_n = top_n,
                       other_label = other_label, weight = "count") %>%
    collapse_long_tail(axis_b, top_n = top_n,
                       other_label = other_label, weight = "count") %>%
    dplyr::group_by(.data[[axis_a]], .data[[axis_b]]) %>%
    dplyr::summarise(count = sum(count), .groups = "drop")

  # Strata stack in factor order, first level on top.
  axis_levels <- function(axis) {
    values <- as.character(data[[axis]])
    if (!identical(axis, "species")) return(order_by_count(data, axis, weight = "count"))
    top_first <- rev(species_order(values, ctx$species_tree, ctx$species_order))
    kept <- top_first[top_first %in% values]
    c(kept[!grepl("^Other", kept)], kept[grepl("^Other", kept)])
  }
  data <- data %>%
    dplyr::mutate(
      "{axis_a}" := factor(.data[[axis_a]], levels = axis_levels(axis_a)),
      "{axis_b}" := factor(.data[[axis_b]], levels = axis_levels(axis_b))
    )
  others <- grep("^Other", unique(as.character(data[[fill_axis]])), value = TRUE)
  colours <- c(colours, stats::setNames(rep(.GREY_OTHER, length(others)), others))

  p <- ggplot(data, aes(axis1 = .data[[axis_a]], axis2 = .data[[axis_b]], y = count)) +
    geom_alluvium(aes(fill = .data[[fill_axis]]), alpha = 0.7, width = 0.2) +
    geom_stratum(fill = .PAPER, colour = .GREY_MID, width = 0.2) +
    # Strata under 2% of the ranges are too thin to hold a label.
    geom_text(
      aes(label = after_stat(ifelse(prop > 0.02,
                                    sprintf("%s (%s)", stratum, scales::comma(count)),
                                    ""))),
      stat = "stratum", size = 3, family = .FONT, colour = .INK
    ) +
    scale_x_discrete(limits = axis_titles, expand = c(0.12, 0.12)) +
    scale_fill_manual(values = colours) +
    theme_retroseek_blank() +
    theme(legend.position = "none",
          axis.text.x = element_text(face = "bold", colour = .INK))
  if (is.null(title)) return(p)
  add_titles(p, title = title, subtitle = subtitle, subset_label = subset_label)
}


sankey_species_probe_plot <- function(data, top_n = NULL, other_label = "Other",
                                      subset_label = NULL, ctx = NULL, colours = NULL) {
  if (is.null(colours)) colours <- probe_colours(data$probe)
  .sankey_two_axis(data, "species", "probe", fill_axis = "probe",
                   colours = colours, axis_titles = c("Host", "Probe"),
                   top_n = top_n, other_label = other_label,
                   title    = "From hosts to probes",
                   subtitle = "Merged ranges from each host to the probe that found them, coloured by probe.",
                   subset_label = subset_label, ctx = ctx)
}


sankey_label_probe_plot <- function(data, top_n = NULL, other_label = "Other",
                                    subset_label = NULL, ctx = NULL, colours = NULL) {
  .sankey_two_axis(data, "label", "probe", fill_axis = "label",
                   colours = taxon_colours(data$label), axis_titles = c("Lineage", "Probe"),
                   top_n = top_n, other_label = other_label,
                   title    = "From lineages to probes",
                   subtitle = "Merged ranges from each lineage to the probe that found them, coloured by lineage.",
                   subset_label = subset_label, ctx = ctx)
}


sankey_species_label_plot <- function(data, top_n = NULL, other_label = "Other",
                                      subset_label = NULL, ctx = NULL, colours = NULL) {
  .sankey_two_axis(data, "species", "label", fill_axis = "label",
                   colours = taxon_colours(data$label), axis_titles = c("Host", "Lineage"),
                   top_n = top_n, other_label = other_label,
                   title    = "From hosts to lineages",
                   subtitle = "Merged ranges from each host to the lineage of the probe virus, coloured by lineage.",
                   subset_label = subset_label, ctx = ctx)
}
