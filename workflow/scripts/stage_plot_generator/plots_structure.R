# =============================================================================
# stage_plot_generator/plots_structure.R - LTR structural components
# =============================================================================
# Builders characterising the LTRdigest retrotransposon calls themselves: how
# structurally complete they are (flanking LTRs / Pfam domains / TSD / PPT) and
# which probe-domain combinations they carry.


# Faceted bar panel: retrotransposon counts per structural component. Replaces
# a former 0-3 composite "completeness score" - under default LTRharvest config
# `has_both_ltrs` and `has_tsd` are ~constant (LTRharvest only emits 2-LTR,
# TSD-flanked elements), so summing them into a score added no discrimination.
# Each component is shown on its own facet instead; the genuinely-varying
# signals (`n_selected_domains`, `n_domains_total`, `has_ppt`) carry the
# information, and the near-constant ones are reported honestly as rates in the
# subtitle rather than dressed up as a score.
ltr_structure_components_plot <- function(ltr_df, subset_label = NULL,
                                          warning_caption = NULL) {
  if (nrow(ltr_df) == 0L) return(empty_plot())

  # Bucket a count vector to "0","1","2","3","4+" for compact, readable facets.
  bucket <- function(x) {
    b <- pmin(as.integer(x), 4L)
    ifelse(b == 4L, "4+", as.character(b))
  }
  yesno <- function(x) ifelse(x, "yes", "no")

  components <- c(
    "Flanking LTRs (n)", "Probe-assigned domains (n)",
    "Pfam domains, total (n)", "Polypurine tract (PPT)",
    "Target-site duplication (TSD)"
  )
  d <- dplyr::bind_rows(
    tibble::tibble(component = components[1], category = bucket(ltr_df$n_flanking_ltrs)),
    tibble::tibble(component = components[2], category = bucket(ltr_df$n_selected_domains)),
    tibble::tibble(component = components[3], category = bucket(ltr_df$n_domains_total)),
    tibble::tibble(component = components[4], category = yesno(ltr_df$has_ppt)),
    tibble::tibble(component = components[5], category = yesno(ltr_df$has_tsd))
  ) %>%
    dplyr::count(component, category, name = "count") %>%
    dplyr::mutate(
      category  = factor(category,
                          levels = c("0", "1", "2", "3", "4+", "no", "yes")),
      component = factor(component, levels = components)
    )

  both_rate <- mean(ltr_df$has_both_ltrs)
  tsd_rate  <- mean(ltr_df$has_tsd)
  ppt_rate  <- mean(ltr_df$has_ppt)

  p <- ggplot(d, aes(x = category, y = count, fill = component)) +
    geom_col(colour = "black", linewidth = 0.2) +
    facet_wrap(~ component, scales = "free_x", ncol = 3) +
    scale_fill_brewer(palette = "YlGnBu") +
    theme_minimal() +
    labs(x = "Per-retrotransposon value", y = "Retrotransposons") +
    theme(text = element_text(face = "bold"), legend.position = "none")
  add_titles(
    p,
    title    = "LTR structural components",
    subtitle = sprintf(
      paste0("Both flanking LTRs %.0f%% - TSD %.0f%% - PPT %.0f%% - ",
             "the first two are ~constant under default LTRharvest config"),
      100 * both_rate, 100 * tsd_rate, 100 * ppt_rate),
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}


# Bar: which curated domain-class combinations occur across retrotransposons.
# An element carrying only L1 machinery now reads as `non_ltr` rather than
# being counted as POL by the retired name regex (ADR-015). The
# long tail of rare combinations is folded into a single "Other (k)" stratum
# via collapse_long_tail (shared with plot2sort).
domain_composition_plot <- function(ltr_df, subset_label = NULL,
                                    top_n = 20L, warning_caption = NULL) {
  if (nrow(ltr_df) == 0L) return(empty_plot())
  d <- ltr_df %>%
    dplyr::mutate(domain_classes = dplyr::if_else(is.na(domain_classes),
                                                 "(no domains)", domain_classes)) %>%
    dplyr::count(domain_classes, name = "count") %>%
    collapse_long_tail("domain_classes", top_n = top_n, weight = "count") %>%
    dplyr::group_by(domain_classes) %>%
    dplyr::summarise(count = sum(count), .groups = "drop")
  ordered <- order_by_count(d, "domain_classes", weight = "count")
  d <- d %>% dplyr::mutate(domain_classes = factor(domain_classes, levels = ordered))

  p <- ggplot(d, aes(x = domain_classes, y = count, fill = count)) +
    geom_col(colour = "black", linewidth = 0.2) +
    scale_fill_gradient(low = "#fff7ec", high = "#7f0000", trans = "sqrt") +
    theme_minimal() +
    labs(x = "Domain-class combination", y = "Retrotransposons", fill = "Count") +
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
