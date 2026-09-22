# =============================================================================
# plot2sort/helpers.R - pure utilities shared by every plot builder
# =============================================================================
# Palette + ordering + long-tail collapse + empty placeholder + aggregation +
# title styling + dimension auto-scaling. None of these functions touch disk
# or external services; they are deterministic transforms of in-memory tables
# (or, for the ggplot helpers, of ggplot objects).


# Generate `output_colour_number` interpolated colours from the planet-express
# Futurama palette. Used wherever a categorical fill needs more shades than
# the base palette provides (~12).
futurama_unlimited_palette <- function(input_colour_number = 12, output_colour_number) {
  planet_express <- pal_futurama("planetexpress")(input_colour_number)
  output_colour  <- colorRampPalette(planet_express)(output_colour_number)
  return(output_colour)
}


# IGV categorical palette that never runs out. ggsci's IGV palette caps at 51
# discrete colours, so `scale_fill_igv()` errors when a fill has more levels
# (e.g. 102 host species). Return the exact IGV colours for n <= 51 (identical
# appearance to scale_fill_igv) and interpolate beyond, so high-cardinality
# per-species panels render instead of aborting the whole plot stage.
igv_unlimited_palette <- function(n) {
  n <- max(as.integer(n), 1L)
  igv <- ggsci::pal_igv("default")(min(n, 51L))
  if (n <= 51L) igv else grDevices::colorRampPalette(igv)(n)
}


# Map genome FASTA stems to their display names from the config `species:` map
# (stem -> "Display name"). Unmapped stems pass through unchanged, so the plot
# still renders if a genome is missing from the map. Vectorised; the caller
# assigns the result back to a `species` / `genome` column before plotting.
relabel_species <- function(values, species_map) {
  values <- as.character(values)
  if (is.null(species_map) || length(species_map) == 0L) return(values)
  vapply(values, function(v) {
    nm <- species_map[[v]]
    if (is.null(nm) || !nzchar(as.character(nm))) v else as.character(nm)
  }, character(1), USE.NAMES = FALSE)
}


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


# Keep only the top-N strata of `col` by count (or summed `weight`); fold the
# rest into a single labelled "Other (k)" stratum that records how many strata
# were collapsed. `top_n = NULL | NA | <=0` short-circuits and returns `df`
# unchanged - the default "show every stratum" behaviour.
#
# Note: returns a possibly-non-aggregated frame; callers that grouped on `col`
# should re-aggregate after calling this so duplicate "Other" rows fold.
collapse_long_tail <- function(df, col, top_n, other_label = "Other",
                               weight = NULL) {
  if (is.null(top_n) || is.na(top_n) || top_n <= 0L) return(df)
  if (nrow(df) == 0L) return(df)
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


# Placeholder for zero-row inputs, so a stage still writes its page and the DAG
# keeps flowing on genomes with no hits. Styled like every other page.
empty_plot <- function(label = "No data") {
  ggplot() +
    theme_void(base_family = .FONT) +
    labs(title = label) +
    theme(
      plot.title.position = "plot",
      plot.title      = element_text(hjust = 0, face = "bold", size = 15, colour = .INK_SOFT),
      plot.background = element_rect(fill = .PAPER, colour = NA),
      plot.margin     = margin(14, 18, 12, 14)
    )
}


# Title and subtitle for a plot, in the house style (style.R): left-aligned,
# sentence case. `subset_label` names what the page is about (a genome, a probe
# set, a segment); it leads the subtitle rather than being glued onto the title
# with a dash, so titles stay short and identical across genomes. Genome stems
# must already be readable names here (display_species), never file names.
# `warning_caption`, when supplied, stamps a caveat that travels with the page.
add_titles <- function(p, title, subtitle, subset_label = NULL,
                       warning_caption = NULL) {
  lead <- if (!is.null(subset_label) && nzchar(subset_label)) subset_label else NULL
  full_subtitle <- if (!is.null(lead) && !is.null(subtitle) && nzchar(subtitle)) {
    paste0(lead, ". ", subtitle)
  } else if (!is.null(lead)) {
    lead
  } else {
    subtitle
  }
  p <- p + labs(title = title, subtitle = full_subtitle) +
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
  if (is.null(tier) || !nzchar(tier)) return(p)
  existing <- p$labels$subtitle
  combined <- if (!is.null(existing) && nzchar(existing)) {
    paste0(existing, "\n", tier)
  } else {
    tier
  }
  p + labs(subtitle = combined)
}


# Stamp a bold red warning caption onto a finished plot, or return it
# unchanged when `caption` is NULL/empty. Single styling source shared by
# add_titles() (stage_plot_generator builders) and plot2sort.R's emit()
# wrapper, so the entry-explosion caveat looks identical everywhere.
stamp_warning_caption <- function(p, caption) {
  if (is.null(caption) || !nzchar(caption)) return(p)
  p +
    labs(caption = caption) +
    theme(plot.caption = element_text(hjust = 0, face = "bold",
                                      colour = .WARNING_INK, size = 9.5,
                                      margin = margin(t = 8)))
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
  multi <- c("list", "concatenate")
  offenders <- c(
    if (!is.null(agg$virus) && agg$virus %in% multi) sprintf("virus=%s", agg$virus),
    if (!is.null(agg$label) && agg$label %in% multi) sprintf("label=%s", agg$label)
  )
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


# Auto-scale the (width, height) of a ggsave canvas based on the cardinality
# of a categorical axis. The grow-with-N axis is parameterised so the same
# helper covers x-axis (bar / heatmap) and y-axis (balloon) plots.
#
#   n            number of strata that will appear on the scaled axis.
#   axis         "x" -> width grows with n, height stays at base_h.
#                "y" -> height grows with n, width stays at base_w.
#   base_w/h     fallback canvas (inches) for small inputs.
#   per_stratum  inches added per stratum past the `base_strata` floor.
#   base_strata  number of strata that fit in the base canvas; below this,
#                the canvas stays at base_w / base_h.
#   cap          hard upper bound (inches) so PNGs stay renderable. At
#                300 dpi the default 60in x 18,000 px is the practical limit.
#
# Returns list(w, h) of doubles in inches.
auto_dims <- function(n, axis = c("x", "y"),
                      base_w = 15, base_h = 12,
                      per_stratum = 0.18,
                      base_strata = 12L,
                      cap = 60) {
  axis  <- match.arg(axis)
  base  <- if (axis == "x") base_w else base_h
  extra <- max(0L, n - base_strata) * per_stratum
  scaled <- min(cap, base + extra)
  if (axis == "x") list(w = scaled, h = base_h)
  else             list(w = base_w, h = scaled)
}


# Point size for a categorical axis carrying `n` tick labels. auto_dims() grows
# the CANVAS but not the TEXT, so at high cardinality (102 host genomes) labels
# still collide on a wider page. Shrink linearly from `base_size` once past
# `base_strata`, with a legibility floor - below ~5pt a label is unreadable
# anyway, and the canvas growth has to carry the rest.
categorical_text_size <- function(n, base_size = 11, base_strata = 12L,
                                  floor_size = 5) {
  if (n <= base_strata) return(base_size)
  max(floor_size, base_size - (n - base_strata) * 0.06)
}


# Scale a finished plot to the cardinality of its categorical axis: attach the
# `intended_dims` attribute save_plot() reads, AND apply the matching theme so
# the text scales with the canvas. Rotates x tick labels once they are too dense
# to sit side by side (y labels read horizontally at any n, so they are left
# alone). Wraps auto_dims() rather than replacing it - the 12 existing call
# sites keep their signature.
#
#   n     number of strata on the scaled axis (species, probes, taxa, ...).
#   axis  "x" (width grows) or "y" (height grows).
#   ...   passed through to auto_dims (base_w/base_h/per_stratum/cap).
#
# Returns the plot with `intended_dims` attached, so emit()/save_plot pick the
# canvas up automatically.
scale_categorical_axis <- function(p, n, axis = c("x", "y"),
                                   rotate_at = 20L, ...) {
  axis <- match.arg(axis)
  size <- categorical_text_size(n)
  p <- if (axis == "x" && n >= rotate_at) {
    p + theme(axis.text.x = element_text(size = size, angle = 90,
                                         hjust = 1, vjust = 0.5))
  } else if (axis == "x") {
    p + theme(axis.text.x = element_text(size = size))
  } else {
    p + theme(axis.text.y = element_text(size = size))
  }
  attr(p, "intended_dims") <- auto_dims(n, axis = axis, ...)
  p
}
