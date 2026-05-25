# =============================================================================
# erv_like_plot_generator/plots_order_structure.R — gene order & structure
# =============================================================================
# Builders describing the spatial structure of ERV-like candidates: the
# observed probe order, whether it is canonical, candidate length, member count,
# and the inter-probe gap distribution (a tuning aid for max_join_distance).
# Same builder contract as the other plot modules. Reuses plot2sort/helpers.R.

.CANON_FILL <- c("canonical (kept)"     = "#1b9e77",
                 "rearranged (kept)"    = "#d95f02",
                 "rearranged (dropped)" = "#7570b3")


# Bar of observed 5'->3' main-probe orders (e.g. "GAG; POL; ENV"), most frequent
# first. Long tail optionally folded into "Other (k)" via collapse_long_tail.
gene_order_plot <- function(loci_df, top_n = NULL, other_label = "Other",
                            subset_label = NULL, warning_caption = NULL) {
  if (nrow(loci_df) == 0L) return(empty_plot())
  d <- loci_df %>% dplyr::count(observed_order, name = "count")
  d <- collapse_long_tail(d, "observed_order", top_n, other_label,
                          weight = "count") %>%
    dplyr::group_by(observed_order) %>%
    dplyr::summarise(count = sum(count), .groups = "drop")
  ord <- order_by_count(d, "observed_order", weight = "count")
  d <- d %>% dplyr::mutate(observed_order = factor(observed_order, levels = ord))

  p <- ggplot(d, aes(x = observed_order, y = count)) +
    geom_col(fill = "#386cb0", colour = "black", linewidth = 0.2) +
    scale_y_continuous(labels = scales::label_comma()) +
    theme_minimal() +
    labs(x = "Observed main-probe order (5'→3')", y = "ERV-like candidates") +
    theme(text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  out <- add_titles(
    p,
    title    = "ERV-like observed gene order",
    subtitle = "Frequency of each observed 5'→3' main-probe order",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
  attr(out, "intended_dims") <- auto_dims(nrow(d), axis = "x")
  out
}


# Canonical vs rearranged. Robust to `require_canonical_order`: retained
# candidates split into canonical / rearranged by `is_canonical`, plus the
# rearranged candidates that were *dropped* by the filter (from the counts
# table). When the filter is off, dropped = 0 and both retained bars show;
# when on, the rearranged-kept bar is empty and the dropped bar carries them.
canonical_plot <- function(loci_df, dropped_noncanonical = 0L,
                           subset_label = NULL, warning_caption = NULL) {
  dropped_noncanonical <- as.integer(dropped_noncanonical %||% 0L)
  if (nrow(loci_df) == 0L && dropped_noncanonical == 0L) return(empty_plot())
  kept_canon <- if (nrow(loci_df) > 0L) sum(loci_df$is_canonical, na.rm = TRUE) else 0L
  kept_rear  <- if (nrow(loci_df) > 0L) sum(!loci_df$is_canonical, na.rm = TRUE) else 0L
  d <- tibble::tibble(
    category = factor(names(.CANON_FILL), levels = names(.CANON_FILL)),
    count    = c(kept_canon, kept_rear, dropped_noncanonical)
  )

  p <- ggplot(d, aes(x = category, y = count, fill = category)) +
    geom_col(colour = "black", linewidth = 0.2) +
    scale_fill_manual(values = .CANON_FILL, guide = "none") +
    scale_y_continuous(labels = scales::label_comma()) +
    theme_minimal() +
    labs(x = NULL, y = "ERV-like candidates") +
    theme(text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 20, hjust = 1))
  add_titles(
    p,
    title    = "ERV-like canonical vs rearranged gene order",
    subtitle = "Canonical = main probes in main_probes order; dropped = removed by require_canonical_order",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}


# Histogram of candidate span (bp), split by full vs partial. Full ERVs should
# cluster near the expected full-length ERV size.
candidate_length_plot <- function(loci_df, subset_label = NULL,
                                  warning_caption = NULL) {
  if (nrow(loci_df) == 0L) return(empty_plot())
  d <- loci_df %>%
    dplyr::mutate(full = factor(ifelse(is_full, "full", "partial"),
                                levels = c("partial", "full")))
  p <- ggplot(d, aes(x = width, fill = full)) +
    geom_histogram(bins = 40, colour = NA, alpha = 0.8, position = "stack") +
    scale_x_continuous(labels = scales::label_comma()) +
    scale_fill_manual(values = c(partial = "#d95f02", full = "#1b9e77")) +
    theme_minimal() +
    labs(x = "ERV-like candidate span (bp)", y = "Candidates", fill = NULL) +
    theme(text = element_text(face = "bold"))
  add_titles(
    p,
    title    = "ERV-like candidate length distribution",
    subtitle = "Span from first to last constituent main-probe locus",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}


# Frequency bar of the number of constituent loci per candidate (small integer).
n_loci_plot <- function(loci_df, subset_label = NULL, warning_caption = NULL) {
  if (nrow(loci_df) == 0L) return(empty_plot())
  d <- loci_df %>% dplyr::count(n_loci, name = "count")
  p <- ggplot(d, aes(x = factor(n_loci), y = count)) +
    geom_col(fill = "#386cb0", colour = "black", linewidth = 0.2) +
    scale_y_continuous(labels = scales::label_comma()) +
    theme_minimal() +
    labs(x = "Constituent loci per candidate", y = "ERV-like candidates") +
    theme(text = element_text(face = "bold"))
  add_titles(
    p,
    title    = "ERV-like constituent-locus count",
    subtitle = "Number of main-probe loci chained into each candidate",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}


# Histogram of the gap (bp) between adjacent main-probe loci within candidates,
# with the configured `max_join_distance` marked. A tuning aid: candidates
# piling up against the cutoff suggest raising it would capture more; a
# distribution well below it means the cutoff is not binding.
interprobe_gap_plot <- function(members_df, max_join_distance = NULL,
                                subset_label = NULL, warning_caption = NULL) {
  if (nrow(members_df) == 0L) return(empty_plot())
  d <- members_df %>% dplyr::filter(!is.na(gap_to_prev))
  if (nrow(d) == 0L) return(empty_plot("no adjacent member pairs"))

  p <- ggplot(d, aes(x = gap_to_prev)) +
    geom_histogram(bins = 40, fill = "#386cb0", colour = NA, alpha = 0.85) +
    scale_x_continuous(labels = scales::label_comma()) +
    theme_minimal() +
    labs(x = "Gap between adjacent main-probe loci (bp)", y = "Member pairs") +
    theme(text = element_text(face = "bold"))
  if (!is.null(max_join_distance)) {
    p <- p +
      geom_vline(xintercept = max_join_distance, linetype = "dashed",
                 colour = "#b22222", linewidth = 0.5) +
      annotate("text", x = max_join_distance, y = Inf,
               label = sprintf("max_join_distance = %s", max_join_distance),
               hjust = 1.05, vjust = 1.5, size = 3.5, fontface = "bold",
               colour = "#b22222")
  }
  add_titles(
    p,
    title    = "ERV-like inter-probe gap distribution",
    subtitle = "Gaps between chained adjacent main-probe loci, vs the join cutoff",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}
