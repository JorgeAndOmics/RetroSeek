# ----------------------------
# PLOT GENERATION SCRIPT
# ----------------------------
#
# Builds RetroSeek's global PNG plots from the per-genome plot dataframes
# that ranges_analysis.R writes to results/tables/plot_dataframes/.
#
# Structure:
#   1. Pure helpers   — ordering, long-tail collapse, empty-input placeholder.
#   2. Plot builders  — one function per plot, all take a tibble + config knobs.
#   3. main()         — CLI entry; reads parquets, joins species map, calls builders.
#   4. Bottom guard   — main() runs only under Rscript (sys.nframe() == 0L);
#                       testthat sources this file and only sees the function defs.

options(warn = -1)  # Suppress warnings

suppressMessages({
  library(argparse)     # Command-line argument parser
  library(arrow)        # Parquet I/O
  library(tidyverse)    # Data manipulation and visualisation
  library(yaml)         # YAML config
  library(ggsci)        # Scientific colour palettes (planet express)
  library(ggalluvial)   # Sankey / alluvial geoms
  library(ggdist)       # Raincloud / halfeye / dots
  library(scales)       # Axis labellers
})


# ============================================================================
# 1. PURE HELPERS
# ============================================================================

# Generate `output_colour_number` interpolated colours from the planet-express
# Futurama palette. Used wherever a categorical fill needs more shades than the
# base palette provides (~12).
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
empty_plot <- function(label = "no data") {
  ggplot() +
    theme_void() +
    labs(title = label) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
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


# ============================================================================
# 2. PLOT BUILDERS
# ============================================================================

density_bitscore_plot <- function(data, q1, median, q3, x_scale = "linear") {
  if (nrow(data) == 0L) return(empty_plot())
  manual_colours <- futurama_unlimited_palette(3, length(unique(data$probe)))

  # Unweighted density over max_bitscore — see commit history for the prior
  # weighted version. If per-range identity needs visual emphasis, prefer a
  # separate plot rather than a weight aesthetic.
  p <- ggplot(data, aes(x = max_bitscore, fill = probe, colour = probe)) +
    geom_density(alpha = 0.4, adjust = 3) +
    scale_fill_manual(values = manual_colours) +
    scale_color_manual(values = manual_colours) +
    geom_vline(xintercept = c(q1, median, q3),
               linetype = "dashed", linewidth = 0.2) +
    annotate("text", x = q1 + 20,     y = 1e-4, label = "Q1") +
    annotate("text", x = median + 20, y = 1e-4, label = "Q2") +
    annotate("text", x = q3 + 20,     y = 1e-4, label = "Q3") +
    labs(x = "HSP Bitscore (max per range)", y = "Density") +
    theme_minimal() +
    theme(
      text         = element_text(face = "bold"),
      axis.title.x = element_text(size = 15, margin = margin(b = 15)),
      axis.title.y = element_text(size = 15, margin = margin(l = 10)),
      axis.text    = element_text(size = 12)
    )
  if (identical(x_scale, "log10")) {
    p + scale_x_log10()
  } else {
    p + scale_x_continuous(
      breaks = seq(0, max(data$max_bitscore, na.rm = TRUE), by = 100)
    )
  }
}


raincloud_bitscore_plot <- function(data, x_scale = "linear") {
  if (nrow(data) == 0L) return(empty_plot())
  manual_colours <- futurama_unlimited_palette(3, length(unique(data$probe)))

  p <- ggplot(data, aes(y = probe, x = max_bitscore,
                        fill = probe, color = probe)) +
    ggdist::stat_halfeye(adjust = 0.5, justification = 0,
                         alpha = 0.5, .width = c(0.5, 0.95)) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
    ggdist::stat_dots(side = "right", dotsize = 0.1,
                      alpha = 0.01, binwidth = 0.2) +
    scale_fill_manual(values = manual_colours) +
    scale_color_manual(values = manual_colours) +
    theme_minimal() +
    labs(x = "HSP Bitscore (max per range)", y = "Probe") +
    theme(
      text       = element_text(face = "bold"),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 11)
    )
  if (identical(x_scale, "log10")) p + scale_x_log10() else p
}


bar_plot <- function(data) {
  if (nrow(data) == 0L) return(empty_plot())
  data <- data %>%
    group_by(species, label) %>%
    summarise(count = sum(count), .groups = "drop")

  ordered_species <- order_by_count(data, "species", weight = "count")
  data <- data %>%
    mutate(species = factor(species, levels = ordered_species))

  ggplot(data) +
    aes(x = species, y = count, fill = label, alpha = count) +
    geom_col(color = "black", linewidth = 0.2) +
    scale_fill_futurama() +
    scale_alpha_continuous(range = c(0.5, 1)) +
    theme_void() +
    labs(fill = "Group", y = "Count") +
    theme(
      axis.title  = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 11)
    ) +
    guides(alpha = "none")
}


balloon_virus_species_plot <- function(data) {
  if (nrow(data) == 0L) return(empty_plot())
  manual_colours <- futurama_unlimited_palette(5, length(unique(data$probe)))

  ordered_species <- order_by_count(data, "species", weight = "count")
  data <- data %>%
    dplyr::mutate(species = factor(species, levels = rev(ordered_species)))

  ggplot(data, aes(x = abbreviation, y = species)) +
    geom_point(aes(color = probe, fill = probe,
                   alpha = count, size = count^3), shape = 15) +
    facet_wrap(~ probe) +
    scale_color_manual(values = manual_colours) +
    scale_fill_manual(values = manual_colours) +
    scale_alpha_continuous(range = c(0.4, 1)) +
    theme_minimal() +
    labs(x = "Virus", y = "Species") +
    theme(
      text        = element_text(face = "bold"),
      axis.title  = element_text(size = 12),
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 9)
    ) +
    guides(size = "none")
}


# Generic two-axis Sankey builder shared by the three exposed wrappers below.
# Applies long-tail collapse + count ordering on each axis, then re-aggregates
# in case "Other" stratum folded multiple rows into the same (axis_a, axis_b)
# combination.
.sankey_two_axis <- function(data, axis_a, axis_b, fill_axis,
                             top_n = NULL, other_label = "Other") {
  if (nrow(data) == 0L) return(empty_plot())

  data <- data %>%
    collapse_long_tail(axis_a, top_n = top_n,
                       other_label = other_label, weight = "count") %>%
    collapse_long_tail(axis_b, top_n = top_n,
                       other_label = other_label, weight = "count") %>%
    dplyr::group_by(.data[[axis_a]], .data[[axis_b]]) %>%
    dplyr::summarise(count = sum(count), .groups = "drop")

  ordered_a <- order_by_count(data, axis_a, weight = "count")
  ordered_b <- order_by_count(data, axis_b, weight = "count")

  data <- data %>%
    dplyr::mutate(
      "{axis_a}" := factor(.data[[axis_a]], levels = ordered_a),
      "{axis_b}" := factor(.data[[axis_b]], levels = ordered_b)
    )

  num_colours    <- max(length(ordered_a), length(ordered_b))
  manual_colours <- futurama_unlimited_palette(12, num_colours)

  ggplot(data, aes(axis1 = .data[[axis_a]],
                   axis2 = .data[[axis_b]],
                   y = count)) +
    geom_alluvium(aes(fill = .data[[fill_axis]], alpha = count)) +
    geom_stratum(alpha = 0, color = "black", linewidth = 0.2) +
    geom_text(
      aes(
        size  = after_stat(count),
        label = sprintf("%s [%d]", after_stat(stratum), after_stat(count))
      ),
      stat      = "stratum",
      nudge_x   = -0.15,
      fontface  = "bold",
      hjust     = 0,
      direction = "y"
    ) +
    scale_size_continuous(range = c(2, 5)) +
    scale_alpha_continuous(range = c(0.4, 1)) +
    scale_fill_manual(values = manual_colours) +
    theme_void() +
    theme(legend.position = "none")
}


sankey_species_probe_plot <- function(data, top_n = NULL,
                                      other_label = "Other") {
  .sankey_two_axis(data, "species", "probe", fill_axis = "species",
                   top_n = top_n, other_label = other_label)
}


sankey_label_probe_plot <- function(data, top_n = NULL,
                                    other_label = "Other") {
  .sankey_two_axis(data, "label", "probe", fill_axis = "label",
                   top_n = top_n, other_label = other_label)
}


sankey_species_label_plot <- function(data, top_n = NULL,
                                      other_label = "Other") {
  .sankey_two_axis(data, "species", "label", fill_axis = "species",
                   top_n = top_n, other_label = other_label)
}


# Distribution of `query_coverage` per probe — reveals under-aligned probes
# whose hits cover only a fraction of the probe sequence (candidates for pHMM
# follow-up). Input: per-range tibble (one row = one merged range).
query_coverage_plot <- function(data) {
  if (nrow(data) == 0L) return(empty_plot())
  if (!"query_coverage" %in% colnames(data)) {
    return(empty_plot("query_coverage column missing"))
  }
  if (all(is.na(data$query_coverage))) {
    return(empty_plot("query_coverage all NA"))
  }
  manual_colours <- futurama_unlimited_palette(3, length(unique(data$probe)))

  ggplot(data, aes(x = query_coverage, fill = probe, colour = probe)) +
    geom_density(alpha = 0.4, adjust = 1.5) +
    scale_fill_manual(values = manual_colours) +
    scale_color_manual(values = manual_colours) +
    labs(x = "Query coverage (alignment length / probe length)",
         y = "Density") +
    theme_minimal() +
    theme(
      text         = element_text(face = "bold"),
      axis.title.x = element_text(size = 15, margin = margin(b = 15)),
      axis.title.y = element_text(size = 15, margin = margin(l = 10)),
      axis.text    = element_text(size = 12)
    )
}


# Cross-species heatmap — `geom_tile`, fill = hit count, both axes ordered by
# marginal totals. Missing (probe, species) cells expand to zero so users see
# absences rather than missing columns.
heatmap_probe_species_plot <- function(data) {
  if (nrow(data) == 0L) return(empty_plot())

  counts <- data %>%
    dplyr::count(species, probe, name = "count") %>%
    tidyr::complete(species, probe, fill = list(count = 0L))

  ordered_species <- counts %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(.t = sum(count), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(.t), species) %>%
    dplyr::pull(species) %>%
    as.character()
  ordered_probe <- counts %>%
    dplyr::group_by(probe) %>%
    dplyr::summarise(.t = sum(count), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(.t), probe) %>%
    dplyr::pull(probe) %>%
    as.character()

  counts <- counts %>%
    dplyr::mutate(
      species = factor(species, levels = ordered_species),
      probe   = factor(probe,   levels = ordered_probe)
    )

  ggplot(counts, aes(x = species, y = probe, fill = count)) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_text(aes(label = count), size = 3, fontface = "bold") +
    scale_fill_gradient(low = "#fff7ec", high = "#7f0000",
                        trans = "sqrt") +
    labs(x = "Species", y = "Probe", fill = "Hits") +
    theme_minimal() +
    theme(
      text        = element_text(face = "bold"),
      axis.title  = element_text(size = 12),
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10)
    )
}


# Waffle chart — virus-genus proportions. Each square represents `unit_hits`
# range counts; on huge inputs increase the knob to keep the waffle legible.
waffle_virus_plot <- function(data, unit_hits = 1L) {
  if (nrow(data) == 0L) return(empty_plot())
  if (!requireNamespace("waffle", quietly = TRUE)) {
    return(empty_plot("waffle package unavailable"))
  }
  unit_hits <- max(1L, as.integer(unit_hits))

  counts <- data %>%
    dplyr::count(virus, name = "count") %>%
    dplyr::mutate(squares = pmax(1L, as.integer(count %/% unit_hits))) %>%
    dplyr::arrange(dplyr::desc(count), virus)

  square_vec <- counts$squares
  names(square_vec) <- counts$virus
  manual_colours <- futurama_unlimited_palette(12, length(square_vec))

  caption <- if (unit_hits > 1L) {
    sprintf("1 square = %d hits", unit_hits)
  } else {
    "1 square = 1 hit"
  }

  waffle::waffle(square_vec, rows = max(1L, floor(sqrt(sum(square_vec)))),
                 colors = manual_colours, legend_pos = "right") +
    labs(title = "Hits per virus", caption = caption) +
    theme(plot.caption = element_text(face = "italic"))
}


# ============================================================================
# 3. MAIN ENTRY (CLI)
# ============================================================================

main <- function() {
  parser <- ArgumentParser(
    description = "Generate global plots from RetroSeek per-genome Parquet results"
  )

  parser$add_argument("--input",  required = TRUE,
                      help = "Input directory with one {genome}.parquet per species; each carries a probe_type column (main | accessory).")
  parser$add_argument("--output", required = TRUE,
                      help = "Directory to save output plots")
  parser$add_argument("--config", required = TRUE,
                      help = "YAML config file with plot parameters")

  args <- parser$parse_args()

  cfg <- yaml::read_yaml(args$config)
  plot_dpi      <- cfg$plots$dpi    %||% 300
  plot_height   <- cfg$plots$height %||% 12
  plot_width    <- cfg$plots$width  %||% 15
  x_scale       <- cfg$plots$bitscore_x_scale %||% "linear"
  top_n         <- cfg$plots$sankey_top_n            # NULL = show all
  other_label   <- cfg$plots$sankey_other_label %||% "Other"
  unit_hits     <- cfg$plots$waffle_unit_hits %||% 1L

  cat("Generating RetroSeek plots...\n")

  parquet_files <- list.files(args$input, pattern = "\\.parquet$",
                              full.names = TRUE)
  if (length(parquet_files) == 0) stop("No per-genome plot parquet files found.")

  all.full <- purrr::map_dfr(parquet_files, arrow::read_parquet)
  if (!"probe_type" %in% colnames(all.full)) {
    stop("plot2sort.R: input parquet(s) missing probe_type column. ",
         "Rebuild ranges_analysis outputs — the plot dataframe contract changed.")
  }
  all.main      <- all.full %>% filter(probe_type == "main")
  all.accessory <- all.full %>% filter(probe_type == "accessory")

  # Map accession → species name from the YAML config map.
  species_map <- cfg$species
  species_df <- tibble(
    species      = names(species_map),
    species_name = unname(unlist(species_map))
  )
  attach_species_name <- function(df) {
    df %>%
      left_join(species_df, by = "species") %>%
      select(-species) %>%
      rename(species = species_name)
  }
  all.main      <- attach_species_name(all.main)
  all.accessory <- attach_species_name(all.accessory)
  all.full      <- attach_species_name(all.full)

  bit.main      <- if (nrow(all.main)      > 0L) q_stats(all.main)
                   else list(q1 = NA, median = NA, q3 = NA, mean = NA)
  bit.accessory <- if (nrow(all.accessory) > 0L) q_stats(all.accessory)
                   else list(q1 = NA, median = NA, q3 = NA, mean = NA)

  all.counted_probe       <- group_count(all.full)
  main.counted_probe      <- group_count(all.main)
  accessory.counted_probe <- group_count(all.accessory)

  # Sankey input frames (species×probe, label×probe, species×label).
  data.main.sankey.species.probe <- main.counted_probe %>%
    group_by(species, probe) %>% summarise(count = sum(count), .groups = "drop")
  data.accessory.sankey.species.probe <- accessory.counted_probe %>%
    group_by(species, probe) %>% summarise(count = sum(count), .groups = "drop")
  data.main.sankey.label.probe <- main.counted_probe %>%
    group_by(label, probe) %>% summarise(count = sum(count), .groups = "drop")
  data.accessory.sankey.label.probe <- accessory.counted_probe %>%
    group_by(label, probe) %>% summarise(count = sum(count), .groups = "drop")
  data.main.sankey.species.label <- main.counted_probe %>%
    group_by(species, label) %>% summarise(count = sum(count), .groups = "drop")
  data.accessory.sankey.species.label <- accessory.counted_probe %>%
    group_by(species, label) %>% summarise(count = sum(count), .groups = "drop")

  save_plot <- function(name, plot,
                        w = plot_width, h = plot_height, r = plot_dpi) {
    ggsave(filename = file.path(args$output, name), plot = plot,
           width = w, height = h, dpi = r)
  }

  # Density / raincloud / bar / balloon
  save_plot("main_density.png",
            density_bitscore_plot(all.main,
                                  bit.main$q1, bit.main$median, bit.main$q3,
                                  x_scale = x_scale))
  save_plot("accessory_density.png",
            density_bitscore_plot(all.accessory,
                                  bit.accessory$q1, bit.accessory$median, bit.accessory$q3,
                                  x_scale = x_scale))
  save_plot("main_raincloud.png",
            raincloud_bitscore_plot(all.main, x_scale = x_scale))
  save_plot("accessory_raincloud.png",
            raincloud_bitscore_plot(all.accessory, x_scale = x_scale))
  save_plot("full_bar.png",      bar_plot(all.counted_probe))
  save_plot("main_bar.png",      bar_plot(main.counted_probe))
  save_plot("accessory_bar.png", bar_plot(accessory.counted_probe))
  save_plot("main_balloon.png",      balloon_virus_species_plot(main.counted_probe))
  save_plot("accessory_balloon.png", balloon_virus_species_plot(accessory.counted_probe))

  # Sankey (a / b / c × main / accessory)
  save_plot("main_sankey_a.png",
            sankey_species_probe_plot(data.main.sankey.species.probe,
                                      top_n = top_n, other_label = other_label))
  save_plot("main_sankey_b.png",
            sankey_label_probe_plot(data.main.sankey.label.probe,
                                    top_n = top_n, other_label = other_label))
  save_plot("main_sankey_c.png",
            sankey_species_label_plot(data.main.sankey.species.label,
                                      top_n = top_n, other_label = other_label))
  save_plot("accessory_sankey_a.png",
            sankey_species_probe_plot(data.accessory.sankey.species.probe,
                                      top_n = top_n, other_label = other_label))
  save_plot("accessory_sankey_b.png",
            sankey_label_probe_plot(data.accessory.sankey.label.probe,
                                    top_n = top_n, other_label = other_label))
  save_plot("accessory_sankey_c.png",
            sankey_species_label_plot(data.accessory.sankey.species.label,
                                      top_n = top_n, other_label = other_label))

  # New plots: query coverage, probe×species heatmap, virus-genus waffle
  save_plot("main_query_coverage.png",       query_coverage_plot(all.main))
  save_plot("accessory_query_coverage.png",  query_coverage_plot(all.accessory))
  save_plot("main_heatmap_probe_species.png",      heatmap_probe_species_plot(all.main))
  save_plot("accessory_heatmap_probe_species.png", heatmap_probe_species_plot(all.accessory))
  save_plot("main_waffle_virus.png",
            waffle_virus_plot(all.main, unit_hits = unit_hits))
  save_plot("accessory_waffle_virus.png",
            waffle_virus_plot(all.accessory, unit_hits = unit_hits))
}


# ============================================================================
# 4. ENTRY-POINT GUARD
# ============================================================================
# Only invoke main() under `Rscript plot2sort.R …`. testthat sources this file
# inside test functions where sys.nframe() > 0, so unit tests get the helpers
# and builders without firing the CLI.
if (sys.nframe() == 0L) main()
