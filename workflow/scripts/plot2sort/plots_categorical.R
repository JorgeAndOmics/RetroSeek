# =============================================================================
# plot2sort/plots_categorical.R - counts per host, lineage, virus and probe
# =============================================================================
# The categorical homology pages. Hosts are always rows, in the canonical order
# (host tree, else config) carried by `ctx` (panel_ctx() in tree_axis.R), beside
# the host tree when one is configured. Lineages wear their fixed genus colours
# and viruses shades of their genus colour (style.R). Other axes (probes,
# viruses) are ordered by count.


# Threshold past which per-tile heatmap labels become unreadable. Below this,
# `geom_text(count)` is layered; above, only the fill gradient + the axis
# ticks carry information.
.HEATMAP_LABEL_THRESHOLD <- 40L


# Viruses in legend order: grouped by genus (most abundant genus first), then by
# count within the genus, "Other" last. `totals` has virus, label and n.
.virus_levels <- function(totals) {
  genus_rank <- match(totals$label, taxon_levels(totals$label, totals$n))
  totals$virus[order(grepl("^Other", totals$virus), genus_rank, -totals$n)]
}


# Ranges per host, stacked by the lineage label of the probe virus.
bar_plot <- function(data, subset_label = NULL, ctx = NULL) {
  if (nrow(data) == 0L) return(empty_plot())
  data <- data %>%
    group_by(species, label) %>%
    summarise(count = sum(count), .groups = "drop")
  data$label <- taxon_factor(data$label, data$count)

  p <- ggplot(data, aes(x = species, y = count, fill = label)) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
    scale_fill_taxon(data$label, data$count) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(x = NULL, y = "Merged ranges", fill = NULL)
  p <- add_titles(p,
                  title    = "Ranges per host, by lineage",
                  subtitle = "Merged ranges per host, by the lineage of the probe virus that found them.",
                  subset_label = subset_label)
  on_rows(p, data$species, ctx)
}


# Ranges per host coloured by the DETECTED VIRUS (the probeset `virus`, finer
# than the lineage label bar_plot uses). The long virus tail is folded into
# "Other (k)" so the stack stays legible.
# `colours` is the stage's virus palette (virus_colours() over every virus), so a
# virus keeps its shade on every page; computed from this page when NULL.
bar_virus_plot <- function(data, top_n = 15L, subset_label = NULL, ctx = NULL,
                           colours = NULL) {
  if (nrow(data) == 0L) return(empty_plot())
  if (!"label" %in% names(data)) data$label <- NA_character_
  if (is.null(colours)) colours <- virus_colours(data$virus, data$label, data$count)
  data <- collapse_long_tail(data, "virus", top_n = top_n, weight = "count")
  totals <- data %>%
    group_by(virus) %>%
    summarise(n = sum(count), label = dplyr::first(label), .groups = "drop")
  others <- grep("^Other", totals$virus, value = TRUE)
  colours <- c(colours, stats::setNames(rep(.GREY_OTHER, length(others)), others))
  levels_ <- .virus_levels(totals)
  data <- data %>%
    group_by(species, virus) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    mutate(virus = factor(virus, levels = levels_))

  p <- ggplot(data, aes(x = species, y = count, fill = virus)) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
    scale_fill_manual(values = colours, breaks = levels_) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(x = NULL, y = "Merged ranges", fill = NULL) +
    guides(fill = guide_legend(ncol = 3))
  p <- add_titles(p,
                  title    = "Ranges per host, by virus",
                  subtitle = sprintf(paste("The %d probe viruses finding the most ranges, the rest",
                                           "pooled. Each virus is a shade of its lineage's colour."),
                                     top_n),
                  subset_label = subset_label)
  on_rows(p, data$species, ctx)
}


# Ranges per virus and host: bubble area is the range count, colour the virus's
# lineage. Probes are pooled; the heatmap page splits hosts by probe.
balloon_virus_species_plot <- function(data, subset_label = NULL, ctx = NULL) {
  if (nrow(data) == 0L) return(empty_plot())
  d <- data %>%
    group_by(species, abbreviation, label) %>%
    summarise(count = sum(count), .groups = "drop")
  totals <- d %>% group_by(virus = abbreviation) %>%
    summarise(n = sum(count), label = dplyr::first(label), .groups = "drop")
  d$abbreviation <- factor(d$abbreviation, levels = .virus_levels(totals))

  p <- ggplot(d, aes(x = abbreviation, y = species, size = count, colour = label)) +
    geom_point(alpha = 0.85) +
    scale_size_area(max_size = 11, labels = scales::label_comma(), name = "Ranges") +
    scale_colour_manual(values = taxon_colours(d$label), labels = taxon_labels,
                        breaks = taxon_levels(d$label, d$count), name = NULL) +
    labs(x = NULL, y = NULL) +
    # Virus abbreviations are a dense axis of short codes: the one place a
    # tilted label is the lesser evil.
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major = element_line(colour = .GRID, linewidth = 0.3))
  p <- add_titles(p,
                  title    = "Ranges per virus and host",
                  subtitle = paste("Bubble area is the number of merged ranges; colour the virus's",
                                   "lineage. All probes of the set pooled."),
                  subset_label = subset_label)
  on_rows(p, d$species, ctx, axis = "y")
}


# Probe by host heatmap: fill is the range count, zero-hit cells stay visible.
# Probes ordered by total count; hosts are rows in canonical order. Cell labels
# are dropped past `.HEATMAP_LABEL_THRESHOLD` hosts, where they would overprint.
heatmap_probe_species_plot <- function(data, subset_label = NULL, ctx = NULL) {
  if (nrow(data) == 0L) return(empty_plot())

  counts <- data %>%
    dplyr::count(species, probe, name = "count") %>%
    tidyr::complete(species, probe, fill = list(count = 0L))
  ordered_probe <- order_by_count(counts, "probe", weight = "count")
  counts <- counts %>% dplyr::mutate(probe = factor(probe, levels = ordered_probe))
  counts$ink <- ink_on_ramp(counts$count, trans = sqrt)

  n_species <- dplyr::n_distinct(counts$species)
  show_cell_labels <- n_species <= .HEATMAP_LABEL_THRESHOLD

  p <- ggplot(counts, aes(x = probe, y = species, fill = count)) +
    geom_tile(colour = .PAPER, linewidth = 0.6) +
    scale_fill_ramp(trans = "sqrt", labels = scales::label_comma(), name = "Ranges") +
    scale_colour_identity() +
    labs(x = NULL, y = NULL) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(face = "italic"))
  if (show_cell_labels) {
    p <- p + geom_text(aes(label = scales::comma(count), colour = ink), size = 3,
                       family = .FONT)
  }
  subtitle <- paste("Merged ranges per probe and host, square-root colour scale;",
                    "hosts with no hit for a probe stay visible.")
  if (!show_cell_labels) {
    subtitle <- paste(subtitle, sprintf("Cell numbers hidden: %d hosts.", n_species))
  }
  p <- add_titles(p, title = "Ranges per probe and host", subtitle = subtitle,
                  subset_label = subset_label)
  on_rows(p, counts$species, ctx, axis = "y")
}


# Waffle chart - virus proportions. Each square represents `unit_hits` ranges.
# `unit_hits = NULL` (or any value that would yield more than `target_squares`
# total tiles) auto-derives a unit so the waffle stays legible. waffle::waffle()
# degrades badly past a few thousand squares - without this guard,
# production-scale inputs render as an unreadable block.
waffle_virus_plot <- function(data, unit_hits = NULL,
                              target_squares = 400L,
                              subset_label = NULL, colours = NULL) {
  if (nrow(data) == 0L) return(empty_plot())
  if (!requireNamespace("waffle", quietly = TRUE)) {
    return(empty_plot("waffle package unavailable"))
  }
  if (!"label" %in% names(data)) data$label <- NA_character_

  counts <- data %>%
    dplyr::group_by(virus) %>%
    dplyr::summarise(count = dplyr::n(), label = dplyr::first(label), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(count), virus)
  total_hits <- sum(counts$count)

  user_unit <- if (is.null(unit_hits)) NA_integer_ else as.integer(unit_hits)
  auto_unit <- max(1L, as.integer(ceiling(total_hits / target_squares)))
  unit_hits <- if (is.na(user_unit) || user_unit < auto_unit) auto_unit else user_unit
  auto_applied <- is.na(user_unit) || user_unit < auto_unit

  counts <- counts %>%
    dplyr::mutate(squares = pmax(1L, as.integer(count %/% unit_hits)))
  square_vec <- stats::setNames(counts$squares, counts$virus)
  if (is.null(colours)) colours <- virus_colours(counts$virus, counts$label, counts$count)

  unit_caption <- if (unit_hits > 1L) {
    sprintf("1 square = %d hits", unit_hits)
  } else {
    "1 square = 1 hit"
  }
  caption <- if (auto_applied && unit_hits > 1L) {
    sprintf("%s (auto-scaled from %d total hits)", unit_caption, total_hits)
  } else {
    unit_caption
  }

  p <- waffle::waffle(square_vec,
                      rows       = max(1L, floor(sqrt(sum(square_vec)))),
                      colors     = unname(colours[counts$virus]),
                      size       = 0.4,
                      legend_pos = "right") +
    labs(caption = caption) +
    theme_retroseek_blank() +
    theme(legend.position = "right")
  add_titles(p,
             title    = "Ranges per virus",
             subtitle = "Each square is a fixed number of merged ranges; each virus a shade of its lineage's colour.",
             subset_label = subset_label)
}
