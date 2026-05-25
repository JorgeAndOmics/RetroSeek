# =============================================================================
# plot2sort/plots_categorical.R — categorical / proportional plot builders
# =============================================================================
# Plots whose readability scales with input cardinality. Each builder reports
# an `intended_dims` attribute on its return value so the orchestrator can
# pass the right canvas to ggsave (see save_plot in io.R). Returning the
# dims via attr keeps the builder signature pinned to (data, ..., subset_label)
# and avoids a tuple return that would complicate testing.


# Threshold past which per-tile heatmap labels become unreadable. Below this,
# `geom_text(count)` is layered; above, only the fill gradient + the axis
# ticks carry information.
.HEATMAP_LABEL_THRESHOLD <- 40L

# Threshold past which a 45° x-axis label collides with its neighbours. Above
# this, rotate the labels to 90°.
.BAR_VERTICAL_LABEL_THRESHOLD <- 50L


bar_plot <- function(data, subset_label = NULL) {
  if (nrow(data) == 0L) return(empty_plot())
  data <- data %>%
    group_by(species, label) %>%
    summarise(count = sum(count), .groups = "drop")

  ordered_species <- order_by_count(data, "species", weight = "count")
  data <- data %>%
    mutate(species = factor(species, levels = ordered_species))

  n_species <- length(ordered_species)
  label_angle <- if (n_species > .BAR_VERTICAL_LABEL_THRESHOLD) 90 else 45

  p <- ggplot(data) +
    aes(x = species, y = count, fill = label, alpha = count) +
    geom_col(color = "black", linewidth = 0.2) +
    scale_fill_futurama() +
    scale_alpha_continuous(range = c(0.5, 1)) +
    theme_void() +
    labs(fill = "Group", y = "Count") +
    theme(
      axis.title  = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 11, angle = label_angle, hjust = 1),
      axis.text.y = element_text(size = 11)
    ) +
    guides(alpha = "none")
  out <- add_titles(p,
                    title    = "Range counts per species",
                    subtitle = "Stacked by category label; species ordered by total count",
                    subset_label = subset_label)
  attr(out, "intended_dims") <- auto_dims(n_species, axis = "x")
  out
}


balloon_virus_species_plot <- function(data, subset_label = NULL) {
  if (nrow(data) == 0L) return(empty_plot())
  manual_colours <- futurama_unlimited_palette(5, length(unique(data$probe)))

  ordered_species <- order_by_count(data, "species", weight = "count")
  n_species <- length(ordered_species)
  data <- data %>%
    dplyr::mutate(species = factor(species, levels = rev(ordered_species)))

  p <- ggplot(data, aes(x = abbreviation, y = species)) +
    geom_point(aes(color = probe, fill = probe,
                   alpha = count, size = count^3), shape = 15) +
    facet_wrap(~ probe) +
    scale_color_manual(values = manual_colours) +
    scale_fill_manual(values = manual_colours) +
    scale_alpha_continuous(range = c(0.4, 1)) +
    theme_minimal() +
    labs(x = "Virus", y = "Species") +
    theme(
      text        = element_text(face = "bold"),
      axis.title  = element_text(size = 12),
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 9)
    ) +
    guides(size = "none")
  out <- add_titles(p,
                    title    = "Hits across virus × species",
                    subtitle = "Bubble area scales with range count; faceted by probe",
                    subset_label = subset_label)
  attr(out, "intended_dims") <- auto_dims(n_species, axis = "y")
  out
}


# Cross-species heatmap — `geom_tile`, fill = hit count, both axes ordered by
# marginal totals. Missing (probe, species) cells expand to zero so users see
# absences rather than missing columns. Per-cell text labels are dropped past
# `.HEATMAP_LABEL_THRESHOLD` species — the gradient still conveys magnitude
# and the labels would be unreadable.
heatmap_probe_species_plot <- function(data, subset_label = NULL) {
  if (nrow(data) == 0L) return(empty_plot())

  counts <- data %>%
    dplyr::count(species, probe, name = "count") %>%
    tidyr::complete(species, probe, fill = list(count = 0L))

  ordered_species <- order_by_count(counts, "species", weight = "count")
  ordered_probe   <- order_by_count(counts, "probe",   weight = "count")

  counts <- counts %>%
    dplyr::mutate(
      species = factor(species, levels = ordered_species),
      probe   = factor(probe,   levels = ordered_probe)
    )

  n_species <- length(ordered_species)
  show_cell_labels <- n_species <= .HEATMAP_LABEL_THRESHOLD

  p <- ggplot(counts, aes(x = species, y = probe, fill = count)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient(low = "#fff7ec", high = "#7f0000",
                        trans = "sqrt") +
    labs(x = "Species", y = "Probe", fill = "Hits") +
    theme_minimal() +
    theme(
      text        = element_text(face = "bold"),
      axis.title  = element_text(size = 12),
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10)
    )
  if (show_cell_labels) {
    p <- p + geom_text(aes(label = count), size = 3, fontface = "bold")
  }
  subtitle_suffix <- if (show_cell_labels) "" else
    sprintf(" (cell labels hidden — %d species exceeds threshold of %d)",
            n_species, .HEATMAP_LABEL_THRESHOLD)
  out <- add_titles(p,
                    title    = "Hits across probe × species",
                    subtitle = paste0("Tile fill = hit count (sqrt scale); axes ordered by marginal totals; ",
                                      "zero-hit cells visible", subtitle_suffix),
                    subset_label = subset_label)
  attr(out, "intended_dims") <- auto_dims(n_species, axis = "x")
  out
}


# Waffle chart — virus-genus proportions. Each square represents `unit_hits`
# range counts. `unit_hits = NULL` (or any value that would yield more than
# `target_squares` total tiles) auto-derives a unit so the waffle stays
# legible. waffle::waffle() degrades badly past a few thousand squares —
# without this guard, production-scale inputs render as an unreadable block.
# Canvas is fixed; the waffle does its own auto-scaling internally so a
# wider PNG would just add whitespace.
waffle_virus_plot <- function(data, unit_hits = NULL,
                              target_squares = 400L,
                              subset_label = NULL) {
  if (nrow(data) == 0L) return(empty_plot())
  if (!requireNamespace("waffle", quietly = TRUE)) {
    return(empty_plot("waffle package unavailable"))
  }

  counts <- data %>%
    dplyr::count(virus, name = "count") %>%
    dplyr::arrange(dplyr::desc(count), virus)
  total_hits <- sum(counts$count)

  user_unit <- if (is.null(unit_hits)) NA_integer_ else as.integer(unit_hits)
  auto_unit <- max(1L, as.integer(ceiling(total_hits / target_squares)))
  unit_hits <- if (is.na(user_unit) || user_unit < auto_unit) auto_unit else user_unit
  auto_applied <- is.na(user_unit) || user_unit < auto_unit

  counts <- counts %>%
    dplyr::mutate(squares = pmax(1L, as.integer(count %/% unit_hits)))

  square_vec <- counts$squares
  names(square_vec) <- counts$virus
  manual_colours <- futurama_unlimited_palette(12, length(square_vec))

  unit_caption <- if (unit_hits > 1L) {
    sprintf("1 square = %d hits", unit_hits)
  } else {
    "1 square = 1 hit"
  }
  caption <- if (auto_applied && unit_hits > 1L) {
    sprintf("%s (auto-scaled from %d total hits)", unit_caption, total_hits)
  } else {
    unit_caption
  }

  p <- waffle::waffle(square_vec,
                      rows       = max(1L, floor(sqrt(sum(square_vec)))),
                      colors     = manual_colours,
                      legend_pos = "right") +
    labs(caption = caption) +
    theme(plot.caption = element_text(face = "italic"))
  add_titles(p,
             title    = "Hits per virus",
             subtitle = "Waffle of merged-range counts by virus; tile area is proportional to hit count",
             subset_label = subset_label)
}
