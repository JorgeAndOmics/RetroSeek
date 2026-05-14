# =============================================================================
# stage_plot_generator/plots_structure.R — LTR structural completeness
# =============================================================================
# Builders characterising the LTRdigest retrotransposon calls themselves: how
# structurally complete they are (flanking LTRs / internal Pfam domains / TSDs)
# and which probe-domain combinations they carry.


# Bar: retrotransposon counts by a 0–3 completeness score. One point each for
# both flanking LTRs present, ≥1 probe-assigned Pfam domain, and a target-site
# duplication. The per-feature rates are surfaced in the subtitle.
ltr_completeness_plot <- function(ltr_df, subset_label = NULL,
                                  warning_caption = NULL) {
  if (nrow(ltr_df) == 0L) return(empty_plot())
  d <- ltr_df %>%
    dplyr::mutate(
      completeness = as.integer(has_both_ltrs) +
                     as.integer(n_domains > 0L) +
                     as.integer(has_tsd),
      completeness = factor(completeness, levels = 0:3)
    ) %>%
    dplyr::count(completeness, name = "count", .drop = FALSE)

  both_rate <- mean(ltr_df$has_both_ltrs)
  dom_rate  <- mean(ltr_df$n_domains > 0L)
  tsd_rate  <- mean(ltr_df$has_tsd)

  p <- ggplot(d, aes(x = completeness, y = count, fill = completeness)) +
    geom_col(colour = "black", linewidth = 0.2) +
    scale_fill_brewer(palette = "YlGnBu") +
    theme_minimal() +
    labs(x = "Completeness score (0–3)", y = "Retrotransposons") +
    theme(text = element_text(face = "bold"), legend.position = "none")
  add_titles(
    p,
    title    = "LTR structural completeness",
    subtitle = sprintf(
      "Score = both flanking LTRs (%.0f%%) + ≥1 probe-domain (%.0f%%) + TSD (%.0f%%)",
      100 * both_rate, 100 * dom_rate, 100 * tsd_rate),
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}


# Bar: which probe-domain combinations occur across retrotransposons. The
# long tail of rare combinations is folded into a single "Other (k)" stratum
# via collapse_long_tail (shared with plot2sort).
domain_composition_plot <- function(ltr_df, subset_label = NULL,
                                    top_n = 20L, warning_caption = NULL) {
  if (nrow(ltr_df) == 0L) return(empty_plot())
  d <- ltr_df %>%
    dplyr::mutate(domain_probes = dplyr::if_else(is.na(domain_probes),
                                                 "(no domains)", domain_probes)) %>%
    dplyr::count(domain_probes, name = "count") %>%
    collapse_long_tail("domain_probes", top_n = top_n, weight = "count") %>%
    dplyr::group_by(domain_probes) %>%
    dplyr::summarise(count = sum(count), .groups = "drop")
  ordered <- order_by_count(d, "domain_probes", weight = "count")
  d <- d %>% dplyr::mutate(domain_probes = factor(domain_probes, levels = ordered))

  p <- ggplot(d, aes(x = domain_probes, y = count, fill = count)) +
    geom_col(colour = "black", linewidth = 0.2) +
    scale_fill_gradient(low = "#fff7ec", high = "#7f0000", trans = "sqrt") +
    theme_minimal() +
    labs(x = "Probe-domain combination", y = "Retrotransposons", fill = "Count") +
    theme(text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  out <- add_titles(
    p,
    title    = "Retrotransposon probe-domain composition",
    subtitle = "Distinct sets of probe-assigned Pfam domains per retrotransposon",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
  attr(out, "intended_dims") <- auto_dims(length(ordered), axis = "x")
  out
}
