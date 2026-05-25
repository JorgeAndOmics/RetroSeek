# =============================================================================
# stage_plot_generator/plots_concordance.R — homology ↔ LTR integration
# =============================================================================
# Builders for the homology-vs-LTRdigest middle stage: where tBLASTn loci sit
# relative to LTR retrotransposons, and how each probe's loci thin out through
# the candidate / valid refinement steps. Same builder contract as
# plot2sort/plots_*.R — `(data, ..., subset_label, warning_caption)` →
# ggplot, with an `intended_dims` attr for auto-scaled (probe-axis) plots.

.CONCORDANCE_FILL <- c(inside = "#1b9e77", flanking = "#d95f02",
                       disjoint = "#7570b3")
.STAGE_FILL       <- c(homology = "#7570b3", candidate = "#d95f02",
                       valid = "#1b9e77")


# Stacked bar: per probe, the inside / flanking / disjoint breakdown of
# gr_virus loci relative to LTRdigest retrotransposons. Directly visualises
# the spatial basis of the candidate-hit selection step.
concordance_plot <- function(hits_df, subset_label = NULL,
                             warning_caption = NULL) {
  if (nrow(hits_df) == 0L) return(empty_plot())
  conc_levels <- c("inside", "flanking", "disjoint")
  d <- hits_df %>%
    dplyr::mutate(concordance = factor(concordance, levels = conc_levels)) %>%
    dplyr::count(probe, concordance, name = "count")
  ordered_probe <- order_by_count(d, "probe", weight = "count")
  d <- d %>% dplyr::mutate(probe = factor(probe, levels = ordered_probe))

  p <- ggplot(d, aes(x = probe, y = count, fill = concordance)) +
    geom_col(colour = "black", linewidth = 0.2) +
    scale_fill_manual(values = .CONCORDANCE_FILL) +
    theme_minimal() +
    labs(x = "Probe", y = "gr_virus loci", fill = "Concordance") +
    theme(text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  out <- add_titles(
    p,
    title    = "Homology-locus concordance with LTR elements",
    subtitle = "inside = overlaps a retrotransposon; flanking = within 1 kb; disjoint = farther",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
  attr(out, "intended_dims") <- auto_dims(length(ordered_probe), axis = "x")
  out
}


# Grouped bar: per probe, locus count surviving each refinement stage
# (homology → candidate → valid). `is_candidate` / `is_valid` are nested
# flags, so the bars are monotonically non-increasing within a probe.
probe_yield_plot <- function(hits_df, subset_label = NULL,
                             warning_caption = NULL) {
  if (nrow(hits_df) == 0L) return(empty_plot())
  stage_levels <- c("homology", "candidate", "valid")
  d <- hits_df %>%
    dplyr::group_by(probe) %>%
    dplyr::summarise(
      homology  = dplyr::n(),
      candidate = sum(is_candidate),
      valid     = sum(is_valid),
      .groups   = "drop"
    ) %>%
    tidyr::pivot_longer(c(homology, candidate, valid),
                        names_to = "stage", values_to = "count") %>%
    dplyr::mutate(stage = factor(stage, levels = stage_levels))
  ordered_probe <- d %>%
    dplyr::filter(stage == "homology") %>%
    dplyr::arrange(dplyr::desc(count), probe) %>%
    dplyr::pull(probe) %>%
    as.character()
  d <- d %>% dplyr::mutate(probe = factor(probe, levels = ordered_probe))

  p <- ggplot(d, aes(x = probe, y = count, fill = stage)) +
    geom_col(position = position_dodge(width = 0.8), colour = "black",
             linewidth = 0.2) +
    scale_fill_manual(values = .STAGE_FILL) +
    theme_minimal() +
    labs(x = "Probe", y = "Loci", fill = "Stage") +
    theme(text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  out <- add_titles(
    p,
    title    = "Per-probe yield through the refinement funnel",
    subtitle = "Loci surviving homology → candidate (LTR-overlapping) → valid (domain-matched)",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
  attr(out, "intended_dims") <- auto_dims(length(ordered_probe), axis = "x")
  out
}
