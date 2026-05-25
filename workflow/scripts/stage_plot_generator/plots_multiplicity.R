# =============================================================================
# stage_plot_generator/plots_multiplicity.R — locus multiplicity metrics
# =============================================================================
# How many primitive sequences collapse into each reduced locus:
#   M1  raw threshold-passing tBLASTn hits per first-reduced (gr_virus) locus
#       — the `n_hits` column, plottable stratified by tier because
#       candidate / valid are nested subsets of the original tier.
#   M2  per-virus (gr_virus) loci per globally-reduced (gr_global) locus
#       — the `n_loci` column on the reduced dataframe.
# Both metrics span a wide range, so the count axis is log10.

.TIER_FILL <- c(original = "#7570b3", candidate = "#d95f02", valid = "#1b9e77")


# Overlaid frequency bars of M1 (n_hits) for the original / candidate / valid
# tiers. n_hits is a small-integer count (typically single digits), so a
# discrete frequency bar chart is the honest geom — a binned histogram on a
# log axis would scatter ~4 values into mostly-empty bins. Tiers are nested
# subsets, so the bars are drawn overlaid (position identity); if the valid-
# tier bars sit right of the original-tier ones, evidence depth predicts
# survival through refinement.
multiplicity_m1_plot <- function(hits_df, subset_label = NULL,
                                 warning_caption = NULL) {
  if (nrow(hits_df) == 0L) return(empty_plot())
  tiers <- dplyr::bind_rows(
    hits_df %>% dplyr::transmute(n_hits, tier = "original"),
    hits_df %>% dplyr::filter(is_candidate) %>%
      dplyr::transmute(n_hits, tier = "candidate"),
    hits_df %>% dplyr::filter(is_valid) %>%
      dplyr::transmute(n_hits, tier = "valid")
  ) %>%
    dplyr::mutate(tier = factor(tier,
                                levels = c("original", "candidate", "valid"))) %>%
    dplyr::count(tier, n_hits, name = "loci")

  p <- ggplot(tiers, aes(x = n_hits, y = loci, fill = tier)) +
    geom_col(position = "identity", alpha = 0.55, colour = NA) +
    scale_x_continuous(breaks = scales::breaks_pretty()) +
    scale_y_continuous(labels = scales::label_comma()) +
    scale_fill_manual(values = .TIER_FILL) +
    theme_minimal() +
    labs(x = "Raw tBLASTn hits per locus (M1)", y = "Loci", fill = "Tier") +
    theme(text = element_text(face = "bold"))
  add_titles(
    p,
    title    = "Pre-reduction multiplicity (M1) by tier",
    subtitle = "Raw hits collapsed per first-reduced locus; tiers are nested subsets",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}


# Histogram of M2 (n_loci) — per-virus loci collapsed into each per-probe
# global locus by reduce_global.
multiplicity_m2_plot <- function(reduced_df, subset_label = NULL,
                                 warning_caption = NULL) {
  if (nrow(reduced_df) == 0L) return(empty_plot())
  p <- ggplot(reduced_df, aes(x = n_loci)) +
    geom_histogram(bins = 40, fill = "#386cb0", colour = NA, alpha = 0.85) +
    scale_x_log10(labels = scales::label_comma()) +
    theme_minimal() +
    labs(x = "First-reduced loci per global locus (M2, log10)",
         y = "Global loci") +
    theme(text = element_text(face = "bold"))
  add_titles(
    p,
    title    = "Global-reduction multiplicity (M2)",
    subtitle = "Per-virus loci collapsed into each per-probe global locus",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}
