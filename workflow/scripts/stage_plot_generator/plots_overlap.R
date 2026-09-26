# =============================================================================
# stage_plot_generator/plots_overlap.R - pre-reduction overlap / redundancy
# =============================================================================
# How much the UNREDUCED original tier (gr_virus) overlaps itself, and how much
# that redundancy collapses under reduction. These plots make the reduction
# step's effect visible. Same builder contract as the other plot modules.
# Before/after pairs share one convention: before in grey, after in the data
# colour.

.REDUCE_FILL <- function() c(unreduced = .GREY_MID, reduced = .DATA_COLOUR)
.REDUCE_LABELS <- c(unreduced = "Before reduction", reduced = "After reduction")


# Self-overlap degree: how many OTHER unreduced loci each locus overlaps, stacked
# by probe. A long right tail is the positional redundancy reduction collapses.
overlap_degree_plot <- function(overlap_df, subset_label = NULL,
                                warning_caption = NULL, colours = NULL) {
  if (nrow(overlap_df) == 0L) return(empty_plot())
  if (is.null(colours)) colours <- probe_colours(overlap_df$probe)
  d <- overlap_df %>% dplyr::count(probe, overlap_degree, name = "loci")

  p <- ggplot(d, aes(x = overlap_degree, y = loci, fill = probe)) +
    geom_col(colour = NA) +
    scale_x_continuous(breaks = scales::breaks_pretty()) +
    scale_y_continuous(labels = scales::label_comma()) +
    scale_fill_manual(values = colours) +
    labs(x = "Other unreduced loci each locus overlaps", y = "Loci", fill = "Probe")
  add_titles(
    p,
    title    = "How much the loci overlap before reduction",
    subtitle = "The number of other unreduced loci each locus overlaps: the redundancy that reduction collapses.",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}


# The largest reciprocal-overlap fraction per locus, among loci that overlap
# something. 1.0 = a locus fully contained in, or containing, a neighbour.
reciprocal_fraction_plot <- function(overlap_df, subset_label = NULL,
                                     warning_caption = NULL) {
  if (nrow(overlap_df) == 0L) return(empty_plot())
  d <- overlap_df %>% dplyr::filter(overlap_degree > 0L)
  if (nrow(d) == 0L) return(empty_plot("no overlapping loci"))

  p <- ggplot(d, aes(x = max_reciprocal_fraction)) +
    geom_histogram(bins = 30, boundary = 0, fill = .DATA_COLOUR, colour = .PAPER,
                   linewidth = 0.2) +
    # Zoom rather than limit the scale: scale limits would drop the edge bins.
    coord_cartesian(xlim = c(0, 1)) +
    scale_x_continuous(labels = scales::percent) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(x = "Largest reciprocal overlap (intersection over the smaller locus)",
         y = "Loci")
  add_titles(
    p,
    title    = "How deeply overlapping loci overlap",
    subtitle = "Among unreduced loci that overlap another. 100% is full containment.",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}


# Per probe, locus count before (unreduced gr_virus) and after (globally reduced
# gr_global) reduction: the collapse per probe.
reduction_fold_plot <- function(hits_df, reduced_df, subset_label = NULL,
                                warning_caption = NULL) {
  if (nrow(hits_df) == 0L && nrow(reduced_df) == 0L) return(empty_plot())
  unred <- if (nrow(hits_df) > 0L) {
    hits_df %>%
      dplyr::count(probe, name = "count") %>%
      dplyr::mutate(stage = "unreduced")
  } else NULL
  red <- if (nrow(reduced_df) > 0L) {
    reduced_df %>%
      dplyr::count(probe, name = "count") %>%
      dplyr::mutate(stage = "reduced")
  } else NULL
  d <- dplyr::bind_rows(unred, red) %>%
    dplyr::mutate(stage = factor(stage, levels = c("reduced", "unreduced")))
  ordered_probe <- rev(order_by_count(d, "probe", weight = "count"))
  d <- d %>% dplyr::mutate(probe = factor(probe, levels = ordered_probe))

  p <- ggplot(d, aes(x = probe, y = count, fill = stage)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    scale_fill_manual(values = .REDUCE_FILL(), labels = .REDUCE_LABELS,
                      breaks = c("unreduced", "reduced")) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(x = NULL, y = "Loci", fill = NULL)
  add_titles(
    .probe_rows(p),
    title    = "How far reduction collapses each probe",
    subtitle =
      "Loci before (per virus) and after (per probe, globally reduced) the reduction.",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}


# Total range length (bp) summed across loci, before and after reduction. The
# drop is the overlap-redundant length removed by global reduction.
coverage_before_after_plot <- function(coverage_df, subset_label = NULL,
                                       warning_caption = NULL) {
  if (nrow(coverage_df) == 0L) return(empty_plot())
  keys <- c(total_bp_unreduced = "unreduced", total_bp_reduced = "reduced")
  d <- coverage_df %>%
    dplyr::mutate(stage = factor(unname(keys[metric]),
                                 levels = c("reduced", "unreduced")))

  p <- ggplot(d, aes(x = stage, y = value, fill = stage)) +
    geom_col(width = 0.6, show.legend = FALSE) +
    geom_text(aes(label = scales::comma(value)), hjust = -0.1, size = 3.2,
              family = .FONT) +
    coord_flip() +
    scale_x_discrete(labels = .REDUCE_LABELS) +
    scale_fill_manual(values = .REDUCE_FILL()) +
    scale_y_continuous(labels = scales::label_comma(),
                       expand = expansion(mult = c(0, 0.2))) +
    labs(x = NULL, y = "Total range length (bp)") +
    theme(panel.grid.major.y = element_blank())
  add_titles(
    p,
    title    = "Total range length before and after reduction",
    subtitle =
      "Locus widths summed; the drop is the overlapping length the reduction removed.",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}
