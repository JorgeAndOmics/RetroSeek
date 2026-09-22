# =============================================================================
# stage_plot_generator/plots_concordance.R - homology <-> LTR integration
# =============================================================================
# Builders for the homology-vs-LTRdigest middle stage: where tBLASTn loci sit
# relative to LTR retrotransposons, and how each probe's loci thin out through
# the candidate step. Builder contract: `(data, ..., subset_label,
# warning_caption)` -> ggplot. Probes are rows (gene names in italics), ordered
# by count, largest on top.


# Probes as horizontal rows, the largest on top, with the house treatment of a
# categorical row axis. Shared by every probe-axis bar in this stage.
.probe_rows <- function(p) {
  p + coord_flip() +
    theme(axis.text.y = element_text(face = "italic"),
          panel.grid.major.y = element_blank())
}


# Per probe, the inside / flanking / disjoint breakdown of gr_virus loci
# relative to LTRdigest retrotransposons: the spatial basis of the candidate-hit
# selection step. Closer means darker on the ramp.
concordance_plot <- function(hits_df, subset_label = NULL,
                             warning_caption = NULL) {
  if (nrow(hits_df) == 0L) return(empty_plot())
  conc_levels <- c("inside", "flanking", "disjoint")
  d <- hits_df %>%
    dplyr::mutate(concordance = factor(concordance, levels = conc_levels)) %>%
    dplyr::count(probe, concordance, name = "count")
  ordered_probe <- rev(order_by_count(d, "probe", weight = "count"))
  d <- d %>% dplyr::mutate(probe = factor(probe, levels = ordered_probe))

  p <- ggplot(d, aes(x = probe, y = count, fill = concordance)) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
    scale_fill_manual(values = c(inside = seq_colours(4)[4],
                                 flanking = seq_colours(4)[2],
                                 disjoint = .GREY_OTHER),
                      labels = c(inside = "Inside an element",
                                 flanking = "Within 1 kb of one",
                                 disjoint = "Farther away")) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(x = NULL, y = "First-reduced loci", fill = NULL)
  add_titles(
    .probe_rows(p),
    title    = "Where homology loci sit relative to LTR elements",
    subtitle = "Each probe's first-reduced loci, by their position relative to LTRdigest retrotransposons.",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}


# Per probe, the loci found by homology and those that also overlap an LTR
# element (the candidates). Candidates nest inside homology loci, so the second
# bar never exceeds the first.
probe_yield_plot <- function(hits_df, subset_label = NULL,
                             warning_caption = NULL) {
  if (nrow(hits_df) == 0L) return(empty_plot())
  d <- hits_df %>%
    dplyr::group_by(probe) %>%
    dplyr::summarise(
      homology  = dplyr::n(),
      candidate = sum(is_ltr_flanked),
      .groups   = "drop"
    ) %>%
    tidyr::pivot_longer(c(homology, candidate),
                        names_to = "stage", values_to = "count") %>%
    dplyr::mutate(stage = factor(stage, levels = c("candidate", "homology")))
  ordered_probe <- d %>%
    dplyr::filter(stage == "homology") %>%
    dplyr::arrange(count, dplyr::desc(probe)) %>%
    dplyr::pull(probe) %>%
    as.character()
  d <- d %>% dplyr::mutate(probe = factor(probe, levels = ordered_probe))

  p <- ggplot(d, aes(x = probe, y = count, fill = stage)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    scale_fill_manual(values = c(homology = .GREY_MID,
                                 candidate = .TIER_COLOUR[["ltr-flanked"]]),
                      labels = c(homology = "Found by homology",
                                 candidate = "Also overlapping an LTR element"),
                      breaks = c("homology", "candidate")) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(x = NULL, y = "Loci", fill = NULL)
  add_titles(
    .probe_rows(p),
    title    = "What each probe yields",
    subtitle = "Loci found by homology, and those that also overlap an LTR element: the candidates.",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}
