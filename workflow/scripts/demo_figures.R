# =============================================================================
# demo_figures.R — anonymised README showcase figures
# =============================================================================
# Regenerates the figures embedded in README.md from REAL RetroSeek output, with
# host species / provirus / lineage names replaced by neutral placeholders
# (Species A.., Provirus A.., Lineage A..) so the public repo never ships real
# (possibly unpublished) research identifiers. Gene/probe names (POL/GAG/ENV)
# and all numeric values are kept, so the figures stay biologically legible.
#
# Reuses the production plot builders in plot2sort/*.R and
# erv_like_plot_generator/*.R — the demo figures are therefore visually
# identical to real pipeline output, only relabelled and rendered on a lighter
# canvas. Source parquets live outside the repo (gitignored results), so a clean
# clone cannot regenerate these; the committed PNGs are the artifact and this
# script is the refresh tool.
#
# Usage:
#   Rscript demo_figures.R \
#     --input  results/tables/ranges_analysis \  # *.final_loci + *.erv_like_loci parquets
#     --output data/images \
#     --config data/config/config.yaml
#
# The relabel scheme is deterministic (sorted unique value -> letter) and lives
# in make_label_map() below — adjust the prefixes there if desired.

suppressMessages({
  library(argparse)
  library(arrow)
  library(tidyverse)
  library(yaml)
  library(ggsci)
  library(ggalluvial)
  library(ggdist)
  library(scales)
})


# ----------------------------------------------------------------------------
# Locate sibling scripts + source the production builders (no duplication)
# ----------------------------------------------------------------------------
.resolve_script_dir <- function() {
  for (i in rev(seq_len(sys.nframe()))) {
    fr <- tryCatch(sys.frame(i), error = function(e) NULL)
    if (is.null(fr)) next
    ofile <- tryCatch(fr$ofile, error = function(e) NULL)
    if (!is.null(ofile)) return(dirname(normalizePath(ofile, mustWork = FALSE)))
  }
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0L) return(dirname(sub("^--file=", "", file_arg[1])))
  "scripts"
}
.script_dir <- .resolve_script_dir()
source(file.path(.script_dir, "plot2sort", "helpers.R"))
source(file.path(.script_dir, "plot2sort", "io.R"))
source(file.path(.script_dir, "plot2sort", "plots_distribution.R"))
source(file.path(.script_dir, "plot2sort", "plots_categorical.R"))
source(file.path(.script_dir, "plot2sort", "plots_sankey.R"))
source(file.path(.script_dir, "erv_like_plot_generator", "io.R"))
source(file.path(.script_dir, "erv_like_plot_generator", "plots_composition.R"))


# ----------------------------------------------------------------------------
# Instrumentation — io.R helpers call log_section(), so define it (as the
# production orchestrators do) before any loader runs.
# ----------------------------------------------------------------------------
.t0 <- Sys.time()
log_section <- function(name) {
  elapsed <- as.numeric(difftime(Sys.time(), .t0, units = "secs"))
  message(sprintf("[%6.2fs] > %s", elapsed, name))
}


# ----------------------------------------------------------------------------
# Anonymisation
# ----------------------------------------------------------------------------
# Deterministic value -> placeholder map for one column: sorted unique values
# map to "<prefix> A", "<prefix> B", ... (NA dropped). Used for the identifying
# columns only; probe/gene names and numbers are never mapped.
make_label_map <- function(values, prefix) {
  u <- sort(unique(values[!is.na(values)]))
  if (length(u) > length(LETTERS)) {
    stop(sprintf(
      paste0(
        "Cannot anonymise %d unique '%s' values with single-letter labels (max %d). ",
        "Demo figures target a small, representative run — subset --input to a ",
        "lighter cohort (a 100-category showcase plot is unreadable anyway)."
      ),
      length(u), prefix, length(LETTERS)
    ), call. = FALSE)
  }
  setNames(paste(prefix, LETTERS[seq_along(u)]), u)
}

# Apply a named map to a vector; values absent from the map pass through.
apply_map <- function(x, map) {
  out <- unname(map[as.character(x)])
  ifelse(is.na(out), as.character(x), out)
}

# Relabel the identifying columns of a data frame in place.
anonymise <- function(df, maps) {
  if ("species" %in% names(df)) df$species <- apply_map(df$species, maps$species)
  if ("virus" %in% names(df)) df$virus <- apply_map(df$virus, maps$virus)
  if ("abbreviation" %in% names(df)) df$abbreviation <- apply_map(df$abbreviation, maps$abbreviation)
  if ("label" %in% names(df)) df$label <- apply_map(df$label, maps$label)
  df
}


# ----------------------------------------------------------------------------
# main()
# ----------------------------------------------------------------------------
main <- function() {
  parser <- ArgumentParser(description = "Regenerate anonymised README demo figures")
  parser$add_argument("--input", default = "results/tables/ranges_analysis",
                      help = "Dir with *.final_loci.parquet and *.erv_like_loci.parquet")
  parser$add_argument("--output", default = "data/images",
                      help = "Output directory for the demo PNGs")
  parser$add_argument("--config", default = "data/config/config.yaml",
                      help = "YAML config (read for plot parameters)")
  args <- parser$parse_args()

  cfg <- yaml::read_yaml(args$config)
  # Demo render: lighter canvas + dpi than production so README PNGs stay small.
  plot_dpi <- 150
  plot_width <- 11
  plot_height <- 7
  x_scale <- cfg$plots$bitscore_x_scale %||% "linear"
  top_n <- cfg$plots$sankey_top_n
  other_label <- cfg$plots$sankey_other_label %||% "Other"
  unit_hits <- cfg$plots$waffle_unit_hits
  sep <- cfg$parameters$aggregation$concat_separator %||% "; "

  dir.create(args$output, showWarnings = FALSE, recursive = TRUE)

  emit <- function(name, plot) {
    dims <- attr(plot, "intended_dims")
    plot <- stamp_tier_note(plot, "demo data — anonymised")
    save_plot(name, plot, args$output,
              dims = dims, base_w = plot_width, base_h = plot_height, dpi = plot_dpi)
    message("wrote ", file.path(args$output, name))
  }

  # ---- Load + anonymise -----------------------------------------------------
  all.full <- load_plot_dataframes(args$input)
  verify_required_columns(
    all.full,
    c("probe_type", "species", "virus", "probe", "label",
      "abbreviation", "max_bitscore", "query_coverage"),
    source_label = "input parquets"
  )

  # Maps built from the FULL value sets so every panel uses the same scheme.
  maps <- list(
    species      = make_label_map(all.full$species,      "Species"),
    virus        = make_label_map(all.full$virus,        "Provirus"),
    abbreviation = make_label_map(all.full$abbreviation, "Pv"),
    label        = make_label_map(all.full$label,        "Lineage")
  )
  all.full <- anonymise(all.full, maps)
  all.main <- all.full %>% dplyr::filter(probe_type == "main")

  bit.main <- if (nrow(all.main) > 0L) q_stats(all.main)
              else list(q1 = NA, median = NA, q3 = NA, mean = NA)

  all.counted_probe <- group_count(all.full)
  main.counted_probe <- group_count(all.main)

  pair_counts <- function(df, ax_a, ax_b) {
    df %>%
      dplyr::group_by(.data[[ax_a]], .data[[ax_b]]) %>%
      dplyr::summarise(count = sum(count), .groups = "drop")
  }

  # ---- Six README figures ---------------------------------------------------
  # 1. Sankey A: species -> probe
  emit("sankey_a.png",
       sankey_species_probe_plot(pair_counts(main.counted_probe, "species", "probe"),
                                 top_n = top_n, other_label = other_label,
                                 subset_label = "Main"))
  # 2. Waffle: provirus proportions
  emit("waffle.png",
       waffle_virus_plot(all.main, unit_hits = unit_hits, subset_label = "Main"))
  # 3. Bubble / balloon: provirus x species
  emit("balloon.png",
       balloon_virus_species_plot(main.counted_probe, subset_label = "Main"))
  # 4. Raincloud: bitscore distribution by probe
  emit("raincloud.png",
       raincloud_bitscore_plot(all.main, x_scale = x_scale, subset_label = "Main"))
  # 5. Range counts per species (bar, stacked by lineage)
  emit("bar.png",
       bar_plot(all.counted_probe, subset_label = "All probes"))

  # 6. ERV-like composition heatmap: probe x provirus
  loci <- load_erv_like_loci(args$input)
  loci <- anonymise(loci, maps)
  emit("erv_like_heatmap.png",
       composition_heatmap_plot(loci, sep = sep))

  message("Done — wrote 6 anonymised demo figures to ", args$output)
}

if (sys.nframe() == 0L) main()
