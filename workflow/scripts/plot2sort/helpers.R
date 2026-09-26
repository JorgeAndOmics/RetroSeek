# =============================================================================
# plot2sort/helpers.R - pure utilities shared by every plot builder
# =============================================================================
# Ordering + long-tail collapse + empty placeholder + aggregation + titles +
# panel rendering. Colour, type and output format live in style.R. None of
# these functions touch disk or external services; they are deterministic
# transforms of in-memory tables (or, for the ggplot helpers, of ggplot objects).


# Return the levels of `col` ordered by total count (or summed `weight`)
# descending, ties broken alphabetically. Used to factor a column so plots
# render largest-first.
order_by_count <- function(df, col, weight = NULL) {
  if (nrow(df) == 0L) return(character(0))
  if (is.null(weight)) {
    df %>%
      dplyr::count(.data[[col]], name = ".n") %>%
      dplyr::arrange(dplyr::desc(.n), .data[[col]]) %>%
      dplyr::pull(.data[[col]]) %>%
      as.character()
  } else {
    df %>%
      dplyr::group_by(.data[[col]]) %>%
      dplyr::summarise(.n = sum(.data[[weight]], na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(.n), .data[[col]]) %>%
      dplyr::pull(.data[[col]]) %>%
      as.character()
  }
}


# TRUE when `top_n` asks for a cut: a positive number, not NULL or NA.
.keeps_top_n <- function(top_n) {
  !is.null(top_n) && !is.na(top_n) && top_n > 0L
}

# Keep only the top-N strata of `col` by count (or summed `weight`); fold the
# rest into a single labelled "Other (k)" stratum that records how many strata
# were collapsed. `top_n = NULL | NA | <=0` short-circuits and returns `df`
# unchanged - the default "show every stratum" behaviour.
#
# Note: returns a possibly-non-aggregated frame; callers that grouped on `col`
# should re-aggregate after calling this so duplicate "Other" rows fold.
collapse_long_tail <- function(df, col, top_n, other_label = "Other",
                               weight = NULL) {
  if (!.keeps_top_n(top_n) || nrow(df) == 0L) return(df)
  ranking <- order_by_count(df, col, weight = weight)
  if (length(ranking) <= top_n) return(df)
  keep <- ranking[seq_len(top_n)]
  k_collapsed <- length(ranking) - top_n
  label_with_count <- sprintf("%s (%d)", other_label, k_collapsed)
  df %>%
    dplyr::mutate(
      "{col}" := dplyr::if_else(
        as.character(.data[[col]]) %in% keep,
        as.character(.data[[col]]),
        label_with_count
      )
    )
}


# Build every page a panel registry declares, in registry order. Each entry's
# `data` picks its tier scope: "loci" is LTR-flanked only, "combined" is both
# tiers. Losing that distinction would quietly mix orphans into the composition
# and mosaic pages, which deliberately exclude them.
render_panel <- function(registry, loci, combined, ctx) {
  lapply(registry, function(e) {
    e$build(if (identical(e$data, "loci")) loci else combined, ctx)
  })
}


# Placeholder for zero-row inputs, so a stage still writes its page and the DAG
# keeps flowing on genomes with no hits. Styled like every other page.
empty_plot <- function(label = "No data") {
  ggplot() +
    theme_void(base_family = .FONT) +
    labs(title = label) +
    theme(
      plot.title.position = "plot",
      plot.title      = element_text(hjust = 0, face = "bold", size = 15,
                                     colour = .INK_SOFT),
      plot.background = element_rect(fill = .PAPER, colour = NA),
      plot.margin     = margin(14, 18, 12, 14)
    )
}


# TRUE for NULL or an empty string: nothing to print.
.is_blank <- function(x) is.null(x) || !nzchar(x)

# The subtitle with `lead` (a genome, a probe set) in front, "Lead. Subtitle";
# either alone when the other is blank.
.lead_subtitle <- function(lead, subtitle) {
  if (.is_blank(lead)) return(subtitle)
  if (.is_blank(subtitle)) return(lead)
  paste0(lead, ". ", subtitle)
}

# Title and subtitle for a plot, in the house style (style.R): left-aligned,
# sentence case. `subset_label` names what the page is about (a genome, a probe
# set, a segment); it leads the subtitle rather than being glued onto the title
# with a dash, so titles stay short and identical across genomes. Genome stems
# must already be readable names here (display_species), never file names.
# `warning_caption`, when supplied, stamps a caveat that travels with the page.
add_titles <- function(p, title, subtitle, subset_label = NULL,
                       warning_caption = NULL) {
  p <- p + labs(title = title, subtitle = .lead_subtitle(subset_label, subtitle)) +
    # A void-theme page (sankey, placeholder) is otherwise transparent, which some
    # viewers compose on black and which hides the title.
    theme(plot.background = element_rect(fill = .PAPER, colour = NA))
  stamp_warning_caption(p, warning_caption)
}


# Append a neutral tier / reduced-state note to a finished plot's subtitle so
# each PNG is self-documenting about which range tier it shows and whether those
# ranges are reduced (overlaps merged) or non-reduced. Stamped at the
# orchestrator's emit() choke point - no per-builder change needed. Reads the
# builder's existing subtitle from the ggplot object and appends to it; the
# subtitle theme set by add_titles() then styles the whole line uniformly.
stamp_tier_note <- function(p, tier) {
  if (.is_blank(tier)) return(p)
  existing <- p$labels$subtitle
  combined <- if (.is_blank(existing)) tier else paste0(existing, "\n", tier)
  p + labs(subtitle = combined)
}


# Stamp a bold red warning caption onto a finished plot, or return it
# unchanged when `caption` is NULL/empty. Single styling source shared by
# add_titles() (stage_plot_generator builders) and plot2sort.R's emit()
# wrapper, so the entry-explosion caveat looks identical everywhere.
stamp_warning_caption <- function(p, caption) {
  if (.is_blank(caption)) return(p)
  p +
    labs(caption = caption) +
    theme(plot.caption = element_text(hjust = 0, face = "bold",
                                      colour = .WARNING_INK, size = 9.5,
                                      margin = margin(t = 8)))
}


# "name=strategy" when a column aggregates to several values per locus, else NULL.
.multi_value <- function(name, strategy) {
  if (is.null(strategy) || !strategy %in% c("list", "concatenate")) return(NULL)
  sprintf("%s=%s", name, strategy)
}

# Build the entry-explosion warning caption from a parsed config, or return
# NULL when `virus`/`label` use a singular aggregation strategy. `list` and
# `concatenate` produce multi-value cells: `concatenate` explodes one locus
# into N plot rows (count inflation); `list` leaves a compound "A; B; C"
# category label. Either way the aggregate plots are not statistically
# meaningful. Used by plot2sort.R and stage_plot_generator.R.
aggregation_warning <- function(cfg) {
  agg <- cfg$parameters$aggregation
  if (is.null(agg)) return(NULL)
  offenders <- c(.multi_value("virus", agg$virus), .multi_value("label", agg$label))
  if (length(offenders) == 0L) return(NULL)
  sprintf(paste0("Caution: multi-value aggregation is active (%s), so plot counts ",
                 "may be inflated by entry explosion."),
          paste(offenders, collapse = ", "))
}


# Quartile summary over per-range max_bitscore (the strongest single alignment
# inside each merged range). Replaces an earlier mean-of-bitscores which
# averaged on a log-scaled quantity.
q_stats <- function(df) {
  list(
    mean   = mean(df$max_bitscore),
    q1     = quantile(df$max_bitscore, 0.25),
    median = quantile(df$max_bitscore, 0.50),
    q3     = quantile(df$max_bitscore, 0.75)
  )
}


# Aggregate range-level data to per-(species, virus, probe, label, abbreviation)
# counts. Used as the input shape for bar / balloon / sankey plots.
group_count <- function(df) {
  df %>%
    group_by(species, virus, probe, label, abbreviation) %>%
    summarise(count = n(), .groups = "drop")
}
