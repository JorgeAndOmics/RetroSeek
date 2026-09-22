# =============================================================================
# stage_plot_generator/plots_structure.R - LTR structural components
# =============================================================================
# Builders characterising the LTRdigest retrotransposon calls themselves: how
# structurally complete they are (flanking LTRs / Pfam domains / TSD / PPT) and
# which domain combinations they carry.


# One small multiple per structural component. Replaces a former 0-3 composite
# "completeness score": under default LTRharvest config `has_both_ltrs` and
# `has_tsd` are ~constant (LTRharvest only emits 2-LTR, TSD-flanked elements), so
# summing them into a score added no discrimination. The genuinely varying
# signals (`n_domains_total`, `has_ppt`) carry the information, and the
# near-constant ones are reported as rates in the subtitle.
ltr_structure_components_plot <- function(ltr_df, subset_label = NULL,
                                          warning_caption = NULL) {
  if (nrow(ltr_df) == 0L) return(empty_plot())

  # Bucket a count vector to "0","1","2","3","4+" for compact, readable facets.
  bucket <- function(x) {
    b <- pmin(as.integer(x), 4L)
    ifelse(b == 4L, "4+", as.character(b))
  }
  yesno <- function(x) ifelse(x, "Yes", "No")

  components <- c(
    "Flanking LTRs", "Pfam domains",
    "Polypurine tract", "Target-site duplication"
  )
  d <- dplyr::bind_rows(
    tibble::tibble(component = components[1], category = bucket(ltr_df$n_flanking_ltrs)),
    tibble::tibble(component = components[2], category = bucket(ltr_df$n_domains_total)),
    tibble::tibble(component = components[3], category = yesno(ltr_df$has_ppt)),
    tibble::tibble(component = components[4], category = yesno(ltr_df$has_tsd))
  ) %>%
    dplyr::count(component, category, name = "count") %>%
    dplyr::mutate(
      category  = factor(category, levels = c("0", "1", "2", "3", "4+", "No", "Yes")),
      component = factor(component, levels = components)
    )

  p <- ggplot(d, aes(x = category, y = count)) +
    geom_col(fill = .DATA_COLOUR, width = 0.7) +
    facet_wrap(~ component, scales = "free_x", nrow = 1) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(x = NULL, y = "Retrotransposons") +
    theme(panel.grid.major.x = element_blank())
  add_titles(
    p,
    title    = "What each LTR retrotransposon carries",
    subtitle = sprintf(
      paste("Both flanking LTRs in %.0f%%, a target-site duplication in %.0f%%, a",
            "polypurine tract in %.0f%%. The first two are near constant under the",
            "default LTRharvest settings."),
      100 * mean(ltr_df$has_both_ltrs), 100 * mean(ltr_df$has_tsd),
      100 * mean(ltr_df$has_ppt)),
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}


# Which Pfam domain combinations occur across retrotransposons, by raw name.
# Unclassified on purpose: what a domain means is decided once, in the scan, at
# locus grain (ADR-016). The long tail of rare combinations is folded into one
# "Other (k)" bar via collapse_long_tail (shared with plot2sort).
domain_composition_plot <- function(ltr_df, subset_label = NULL,
                                    top_n = 20L, warning_caption = NULL) {
  if (nrow(ltr_df) == 0L) return(empty_plot())
  d <- ltr_df %>%
    dplyr::mutate(element_domains = dplyr::if_else(is.na(element_domains),
                                                 "No domains", element_domains)) %>%
    dplyr::count(element_domains, name = "count") %>%
    collapse_long_tail("element_domains", top_n = top_n, weight = "count") %>%
    dplyr::group_by(element_domains) %>%
    dplyr::summarise(count = sum(count), .groups = "drop")
  ordered <- order_by_count(d, "element_domains", weight = "count")
  # Largest on top, "Other" at the bottom whatever its size.
  ordered <- rev(c(ordered[!grepl("^Other", ordered)], ordered[grepl("^Other", ordered)]))
  d <- d %>% dplyr::mutate(element_domains = factor(element_domains, levels = ordered),
                           other = grepl("^Other", element_domains))

  p <- ggplot(d, aes(x = element_domains, y = count, fill = other)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    coord_flip() +
    # Combinations can list a dozen domains: wrap them so the bars keep their room.
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, 60)) +
    scale_fill_manual(values = c(`FALSE` = .DATA_COLOUR, `TRUE` = .GREY_OTHER)) +
    scale_y_continuous(labels = scales::label_comma(), n.breaks = 4) +
    labs(x = NULL, y = "Retrotransposons") +
    theme(panel.grid.major.y = element_blank())
  add_titles(
    p,
    title    = "Pfam domain combinations",
    subtitle = sprintf("The %d most common sets of Pfam domains per retrotransposon, the rest pooled.",
                       top_n),
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}
