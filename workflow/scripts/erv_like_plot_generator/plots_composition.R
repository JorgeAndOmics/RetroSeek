# =============================================================================
# erv_like_plot_generator/plots_composition.R — composition & completeness
# =============================================================================
# Builders describing WHICH main probes each ERV-like candidate carries and how
# complete it is. Same builder contract as the other plot modules —
# `(data, ..., subset_label, warning_caption)` → ggplot, with an `intended_dims`
# attr for auto-scaled (categorical-axis) plots. Reuses plot2sort/helpers.R.

.FULL_FILL <- c(partial = "#d95f02", full = "#1b9e77")


# `is_full` → an ordered factor "partial" / "full" for stable fill + legend.
.full_factor <- function(is_full) {
  factor(ifelse(is_full, "full", "partial"), levels = c("partial", "full"))
}


# Bar of how often each *set* of main probes is co-assembled (an UpSet-style
# summary without the set-membership matrix — no extra dependency). Each bar is
# one distinct `probes_present` combination (e.g. "GAG; POL", "GAG; POL; ENV").
probe_combination_plot <- function(loci_df, top_n = NULL, other_label = "Other",
                                   subset_label = NULL, warning_caption = NULL) {
  if (nrow(loci_df) == 0L) return(empty_plot())
  d <- loci_df %>% dplyr::count(probes_present, name = "count")
  d <- collapse_long_tail(d, "probes_present", top_n, other_label,
                          weight = "count") %>%
    dplyr::group_by(probes_present) %>%
    dplyr::summarise(count = sum(count), .groups = "drop")
  ord <- order_by_count(d, "probes_present", weight = "count")
  d <- d %>% dplyr::mutate(probes_present = factor(probes_present, levels = ord))

  p <- ggplot(d, aes(x = probes_present, y = count)) +
    geom_col(fill = "#386cb0", colour = "black", linewidth = 0.2) +
    scale_y_continuous(labels = scales::label_comma()) +
    theme_minimal() +
    labs(x = "Main-probe combination", y = "ERV-like candidates") +
    theme(text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  out <- add_titles(
    p,
    title    = "ERV-like probe composition",
    subtitle = "Frequency of each set of main probes co-assembled into a candidate",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
  attr(out, "intended_dims") <- auto_dims(nrow(d), axis = "x")
  out
}


# Bar of candidates by the number of distinct main probes present, split into
# full vs partial. Directly reads off the completeness of the cohort.
completeness_plot <- function(loci_df, subset_label = NULL,
                              warning_caption = NULL) {
  if (nrow(loci_df) == 0L) return(empty_plot())
  d <- loci_df %>%
    dplyr::mutate(full = .full_factor(is_full)) %>%
    dplyr::count(n_main_present, full, name = "count")

  p <- ggplot(d, aes(x = factor(n_main_present), y = count, fill = full)) +
    geom_col(colour = "black", linewidth = 0.2) +
    scale_fill_manual(values = .FULL_FILL) +
    scale_y_continuous(labels = scales::label_comma()) +
    theme_minimal() +
    labs(x = "Distinct main probes present", y = "ERV-like candidates",
         fill = NULL) +
    theme(text = element_text(face = "bold"))
  add_titles(
    p,
    title    = "ERV-like completeness",
    subtitle = "Distinct main probes per candidate; full = meets the completeness threshold",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}


# Stacked bar of full vs partial candidates per virus. `virus` may be a
# multi-value cell under list/concatenate aggregation, so split first.
full_partial_by_virus_plot <- function(loci_df, sep = "; ", subset_label = NULL,
                                       warning_caption = NULL) {
  if (nrow(loci_df) == 0L) return(empty_plot())
  d <- loci_df %>%
    tidyr::separate_rows(virus, sep = sep) %>%
    dplyr::mutate(full = .full_factor(is_full)) %>%
    dplyr::count(virus, full, name = "count")
  ordered_virus <- order_by_count(d, "virus", weight = "count")
  d <- d %>% dplyr::mutate(virus = factor(virus, levels = ordered_virus))

  p <- ggplot(d, aes(x = virus, y = count, fill = full)) +
    geom_col(colour = "black", linewidth = 0.2) +
    scale_fill_manual(values = .FULL_FILL) +
    scale_y_continuous(labels = scales::label_comma()) +
    theme_minimal() +
    labs(x = "Virus", y = "ERV-like candidates", fill = NULL) +
    theme(text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  out <- add_titles(
    p,
    title    = "ERV-like completeness by virus",
    subtitle = "Full vs partial candidates per virus group",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
  attr(out, "intended_dims") <- auto_dims(length(ordered_virus), axis = "x")
  out
}


# Tile heatmap of main-probe presence across viruses. Each candidate's
# `probes_present` set is exploded to one row per (probe, virus); the tile fill
# is the candidate count carrying that probe in that virus.
composition_heatmap_plot <- function(loci_df, sep = "; ", subset_label = NULL,
                                     warning_caption = NULL) {
  if (nrow(loci_df) == 0L) return(empty_plot())
  d <- loci_df %>%
    tidyr::separate_rows(probes_present, sep = sep) %>%
    tidyr::separate_rows(virus, sep = sep) %>%
    dplyr::count(probes_present, virus, name = "count")
  ordered_virus <- order_by_count(d, "virus", weight = "count")
  ordered_probe <- order_by_count(d, "probes_present", weight = "count")
  d <- d %>% dplyr::mutate(
    virus          = factor(virus, levels = ordered_virus),
    probes_present = factor(probes_present, levels = ordered_probe)
  )

  p <- ggplot(d, aes(x = probes_present, y = virus, fill = count)) +
    geom_tile(colour = "grey90") +
    geom_text(aes(label = count), size = 3, fontface = "bold") +
    scale_fill_gradient(low = "#deebf7", high = "#08519c") +
    theme_minimal() +
    labs(x = "Main probe", y = "Virus", fill = "Candidates") +
    theme(text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  out <- add_titles(
    p,
    title    = "ERV-like composition: probe × virus",
    subtitle = "Candidates carrying each main probe, per virus group",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
  attr(out, "intended_dims") <- auto_dims(length(ordered_virus), axis = "y")
  out
}
