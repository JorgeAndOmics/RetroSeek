# =============================================================================
# stage_plot_generator/plots_multiplicity.R - locus multiplicity metrics
# =============================================================================
# How many primitive sequences collapse into each reduced locus:
#   M1  raw threshold-passing tBLASTn hits per first-reduced (gr_virus) locus
#       - the `n_hits` column, plottable stratified by tier because
#       candidate is a nested subset of the original tier.
#   M2  per-virus (gr_virus) loci per globally-reduced (gr_global) locus
#       - the `n_loci` column on the reduced dataframe.


# M1 (n_hits) for the original tier and the LTR-overlapping candidates, side by
# side. n_hits is a small integer, so a bar per value is the honest geom. Shown
# as shares within each tier so the tiers' different sizes do not swamp the
# comparison: if candidates sit further right, evidence depth predicts survival.
multiplicity_m1_plot <- function(hits_df, subset_label = NULL,
                                 warning_caption = NULL) {
  if (nrow(hits_df) == 0L) return(empty_plot())
  tiers <- dplyr::bind_rows(
    hits_df %>% dplyr::transmute(n_hits, tier = "original"),
    hits_df %>%
      dplyr::filter(is_ltr_flanked) %>%
      dplyr::transmute(n_hits, tier = "candidate")
  ) %>%
    dplyr::mutate(tier = factor(tier, levels = c("original", "candidate"))) %>%
    dplyr::count(tier, n_hits, name = "loci") %>%
    dplyr::group_by(tier) %>%
    dplyr::mutate(share = loci / sum(loci)) %>%
    dplyr::ungroup()

  p <- ggplot(tiers, aes(x = n_hits, y = share, fill = tier)) +
    geom_col(position = position_dodge(width = 0.8, preserve = "single"),
             width = 0.75) +
    # A tick per value while there are few, so the first bar (1 hit) is labelled.
    scale_x_continuous(breaks = function(lim) {
      if (diff(lim) <= 12) seq(ceiling(lim[1]), floor(lim[2])) else
        scales::breaks_pretty()(lim)
    }) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c(original = .GREY_MID,
                                 candidate = .TIER_COLOUR[["ltr-flanked"]]),
                      labels = c(original = "Every first-reduced locus",
                                 candidate = "Overlapping an LTR element")) +
    labs(x = "Raw tBLASTn hits per locus (M1)", y = "Share of the tier's loci",
         fill = NULL) +
    theme(panel.grid.major.x = element_blank())
  add_titles(
    p,
    title    = "Hits merged into each locus (M1)",
    subtitle = "Raw hits collapsed into each first-reduced locus, for all loci and for the LTR-overlapping ones.",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}


# M2 (n_loci): per-virus loci collapsed into each per-probe global locus.
multiplicity_m2_plot <- function(reduced_df, subset_label = NULL,
                                 warning_caption = NULL) {
  if (nrow(reduced_df) == 0L) return(empty_plot())
  p <- ggplot(reduced_df, aes(x = n_loci)) +
    geom_histogram(bins = 40, fill = .DATA_COLOUR, colour = NA) +
    scale_x_log10(labels = scales::label_comma()) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(x = "First-reduced loci per global locus (M2, log scale)", y = "Global loci")
  add_titles(
    p,
    title    = "Loci merged by the global reduction (M2)",
    subtitle = "Per-virus loci collapsed into each per-probe global locus.",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}
