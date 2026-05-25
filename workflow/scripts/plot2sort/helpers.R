# =============================================================================
# plot2sort/helpers.R — pure utilities shared by every plot builder
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
# unchanged — the default "show every stratum" behaviour.
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


# Placeholder ggplot for zero-row inputs. Keeps the Snakemake DAG flowing on
# genomes / probe_types with no hits — the rule still produces an output PNG.
# Forces white plot.background so the diagnostic title remains legible in
# viewers that compose transparent PNGs on a dark canvas (same reason as
# add_titles).
empty_plot <- function(label = "no data") {
  ggplot() +
    theme_void() +
    labs(title = label) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}


# Attach a centred title + subtitle to a plot. `subset_label` (e.g. "Main",
# "Accessory") is prepended to the title so the same builder can produce
# differently-named PNGs without duplicating logic. `warning_caption`, when
# supplied, stamps a bold red caption on the plot — used to flag multi-value
# aggregation (entry explosion) so the caveat travels with the PNG artifact.
add_titles <- function(p, title, subtitle, subset_label = NULL,
                       warning_caption = NULL) {
  full_title <- if (!is.null(subset_label) && nzchar(subset_label)) {
    sprintf("%s — %s", subset_label, title)
  } else {
    title
  }
  p <- p +
    labs(title = full_title, subtitle = subtitle) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5, size = 16,
                                     margin = margin(b = 4)),
      plot.subtitle   = element_text(hjust = 0.5, size = 11,
                                     margin = margin(b = 10)),
      # theme_void()-based builders (sankey, bar) leave plot.background as
      # element_blank(), which serialises as transparent in the PNG. Some
      # viewers compose transparent against a dark canvas, making the black
      # title text invisible. Forcing white here pins the background across
      # every builder so titles are always legible.
      plot.background = element_rect(fill = "white", colour = NA)
    )
  stamp_warning_caption(p, warning_caption)
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
                                      colour = "#b22222", size = 10,
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
  sprintf(paste0("⚠ multi-value aggregation active (%s) — plot counts ",
                 "may be inflated by entry explosion; interpret with caution"),
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
#   axis         "x" → width grows with n, height stays at base_h.
#                "y" → height grows with n, width stays at base_w.
#   base_w/h     fallback canvas (inches) for small inputs.
#   per_stratum  inches added per stratum past the `base_strata` floor.
#   base_strata  number of strata that fit in the base canvas; below this,
#                the canvas stays at base_w / base_h.
#   cap          hard upper bound (inches) so PNGs stay renderable. At
#                300 dpi the default 60in × 18,000 px is the practical limit.
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
