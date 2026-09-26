# =============================================================================
# stage_plot_generator/plots_funnel.R - pipeline refinement funnel
# =============================================================================
# Builders visualising how range counts collapse across the pipeline stages.
# Input is the tidy (genome, stage, count) tibble from load_counts_table() - no
# pipeline re-computation; the counts already live in the per-genome counts
# tables. Axes are linear: each step removes a share of the ranges, which a log
# axis would flatten.


# Per-genome funnel: one small multiple per genome, in the canonical species
# order (ctx, see tree_axis.R), with the stages as rows in pipeline order.
refinement_funnel_plot <- function(counts_df, subset_label = NULL,
                                   warning_caption = NULL, ctx = NULL) {
  if (nrow(counts_df) == 0L) return(empty_plot())
  top_first <- rev(species_order(counts_df$genome, ctx$species_tree, ctx$species_order))
  d <- counts_df %>%
    dplyr::mutate(genome = factor(genome, levels = top_first),
                  stage = forcats::fct_rev(stage))
  p <- ggplot(d, aes(x = count, y = stage)) +
    geom_col(fill = .TIER_COLOUR[["ltr-flanked"]], width = 0.65) +
    geom_text(aes(label = scales::comma(count)), hjust = -0.1, size = 3,
              family = .FONT) +
    facet_wrap(~ genome, scales = "free_x") +
    scale_x_continuous(labels = scales::label_comma(), n.breaks = 3,
                       expand = expansion(mult = c(0, 0.35))) +
    labs(x = "Ranges", y = NULL) +
    theme(panel.grid.major.y = element_blank(),
          strip.text = element_text(face = "bold.italic"))
  add_titles(
    p,
    title    = "The refinement funnel per genome",
    subtitle =
      "Ranges left after quality filtering, the first reduction, and the LTR overlap.",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}


# Aggregate cohort funnel: summed counts per stage (bars) with the individual
# genomes overlaid as points, so both the cohort total and its spread show.
aggregate_funnel_plot <- function(counts_df, subset_label = NULL,
                                  warning_caption = NULL) {
  if (nrow(counts_df) == 0L) return(empty_plot())
  totals <- counts_df %>%
    dplyr::group_by(stage) %>%
    dplyr::summarise(total = sum(count, na.rm = TRUE), .groups = "drop")
  p <- ggplot(totals, aes(x = forcats::fct_rev(stage), y = total)) +
    geom_col(fill = .GREY_OTHER, width = 0.65) +
    geom_point(data = counts_df, aes(x = forcats::fct_rev(stage), y = count),
               position = position_jitter(width = 0.12, height = 0, seed = 1),
               size = 1.8, colour = .TIER_COLOUR[["ltr-flanked"]]) +
    coord_flip() +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(x = NULL, y = "Ranges") +
    theme(panel.grid.major.y = element_blank())
  add_titles(
    p,
    title    = "The refinement funnel across the study",
    subtitle = "Bars: ranges summed over every genome. Points: each genome.",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}
