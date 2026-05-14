# =============================================================================
# stage_plot_generator/plots_funnel.R — pipeline refinement funnel
# =============================================================================
# Builders visualising how range counts collapse across the pipeline stages
# (homology → first-reduced → candidate → valid). Input is the tidy
# (genome, stage, count) tibble from load_manifest_counts() — no pipeline
# re-computation; the counts already live in the per-genome manifest YAMLs.
#
# Counts span orders of magnitude (≈10^5 homology hits → ≈10^2 valid ranges),
# so the y-axis is log10. Genomes whose count hits 0 at a stage simply do not
# render a point there — acceptable for a funnel overview.


# Per-genome funnel: one line per genome across the ordered stages.
refinement_funnel_plot <- function(counts_df, subset_label = NULL,
                                   warning_caption = NULL) {
  if (nrow(counts_df) == 0L) return(empty_plot())
  n_genomes <- dplyr::n_distinct(counts_df$genome)
  p <- ggplot(counts_df, aes(x = stage, y = count,
                             group = genome, colour = genome)) +
    geom_line(linewidth = 0.8, alpha = 0.8) +
    geom_point(size = 2) +
    scale_y_log10(labels = scales::label_comma()) +
    scale_colour_manual(values = futurama_unlimited_palette(12, n_genomes)) +
    theme_minimal() +
    labs(x = "Pipeline stage", y = "Range count (log10)", colour = "Genome") +
    theme(text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 25, hjust = 1))
  add_titles(
    p,
    title    = "Refinement funnel per genome",
    subtitle = "Range counts collapsing homology → first-reduced → candidate → valid",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}


# Aggregate cohort funnel: summed counts per stage (bars) with the individual
# genomes overlaid as jittered points, so both the cohort total and its spread
# are visible in one panel.
aggregate_funnel_plot <- function(counts_df, subset_label = NULL,
                                  warning_caption = NULL) {
  if (nrow(counts_df) == 0L) return(empty_plot())
  totals <- counts_df %>%
    dplyr::group_by(stage) %>%
    dplyr::summarise(total = sum(count, na.rm = TRUE), .groups = "drop")
  p <- ggplot(totals, aes(x = stage, y = total)) +
    geom_col(fill = "#386cb0", colour = "black", linewidth = 0.2,
             alpha = 0.85) +
    geom_point(data = counts_df, aes(x = stage, y = count),
               position = position_jitter(width = 0.12, height = 0),
               size = 1.4, alpha = 0.6, colour = "#252525") +
    scale_y_log10(labels = scales::label_comma()) +
    theme_minimal() +
    labs(x = "Pipeline stage", y = "Range count (log10)") +
    theme(text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 25, hjust = 1))
  add_titles(
    p,
    title    = "Refinement funnel across the cohort",
    subtitle = "Bars = summed counts across genomes; points = individual genomes",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}
