# =============================================================================
# demo_figures.R - anonymised README showcase figures
# =============================================================================
# Regenerates the figures embedded in README.md from REAL RetroSeek output, with
# host species / provirus / lineage names replaced by neutral placeholders
# (Species A.., Provirus A.., Lineage A..) so the public repo never ships real
# (possibly unpublished) research identifiers. Gene/probe names (POL/GAG/ENV)
# and all numeric values are kept, so the figures stay biologically legible.
#
# Reuses the production plot builders in plot2sort/*.R and
# erv_like_plot_generator.R - the demo figures are therefore visually identical
# to real pipeline output, only relabelled and rendered on a lighter canvas.
# Source parquets live outside the repo (gitignored results), so a clean clone
# cannot regenerate these; the committed PNGs are the artifact and this script is
# the refresh tool.
#
# Usage:
#   Rscript demo_figures.R \
#     --input  results/tables/ranges_analysis \  # *.final_loci parquets
#     --output data/images \
#     --config data/config/config.yaml
# The erv-like heatmap is read from the sibling results/tables/taxonomy_classification/.
#
# The relabel scheme is deterministic (sorted unique value -> letter) and lives
# in make_label_map() below - adjust the prefixes there if desired.
#
# The figures stay PNG (GitHub renders them inline in the README) but follow the
# house style through the shared builders. Lineages are coloured by their genus
# (style.R), so an anonymised "Lineage A" has no genus colour and draws grey:
# colouring it by its real genus would give the name away, since the genus
# colours are documented.

suppressMessages({
  library(argparse)
  library(arrow)
  library(tidyverse)
  library(yaml)
  library(ggalluvial)
  library(ggdist)
  library(scales)
})


# ----------------------------------------------------------------------------
# Locate sibling scripts + source the production builders (no duplication)
# ----------------------------------------------------------------------------
# Where this script lives, so it can source its siblings. The file name travels
# in `ofile` when the script is source()d (testthat does, several frames deep)
# and in `--file=` when Rscript runs it. The scripts that tests source carry
# this copy; a script cannot source a shared helper before it knows where it
# lives.
.resolve_script_dir <- function() {
  for (frame in rev(sys.frames())) {
    if (!is.null(frame$ofile))
      return(dirname(normalizePath(frame$ofile, mustWork = FALSE)))
  }
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) return(dirname(sub("^--file=", "", file_arg[1])))
  "scripts"
}
.script_dir <- .resolve_script_dir()
# log.R: the line contract and run_main (ADR-021).
# style.R: palette, theme, labels and stage PDFs.
source(file.path(.script_dir, "utils", "log.R"))
source(file.path(.script_dir, "plot2sort", "style.R"))
source(file.path(.script_dir, "plot2sort", "helpers.R"))
source(file.path(.script_dir, "plot2sort", "io.R"))
# tree_axis.R: species rows (the builders use it).
source(file.path(.script_dir, "plot2sort", "tree_axis.R"))
source(file.path(.script_dir, "plot2sort", "plots_distribution.R"))
source(file.path(.script_dir, "plot2sort", "plots_categorical.R"))
source(file.path(.script_dir, "plot2sort", "plots_sankey.R"))
# The erv-like structural panel now reads the genus-founded loci table; reuse its
# loader + composition-heatmap builder (the file's main() stays dormant when sourced).
source(file.path(.script_dir, "taxonomy", "erv_like_plot_generator.R"))


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
        "Demo figures target a small, representative run: subset --input to a ",
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
  if ("abbreviation" %in% names(df))
    df$abbreviation <- apply_map(df$abbreviation, maps$abbreviation)
  if ("label" %in% names(df)) df$label <- apply_map(df$label, maps$label)
  df
}


# ----------------------------------------------------------------------------
# Entry point
# ----------------------------------------------------------------------------
main <- function() {
  parser <- ArgumentParser(description = "Regenerate anonymised README demo figures")
  parser$add_argument(
    "--input", default = "results/tables/ranges_analysis",
    help =
      "Dir with *.final_loci.parquet (taxon loci come from taxonomy_classification/)"
  )
  parser$add_argument("--output", default = "data/images",
                      help = "Output directory for the demo PNGs")
  parser$add_argument("--config", default = "data/config/config.yaml",
                      help = "YAML config (read for plot parameters)")
  parser$add_argument("--log", default = NULL,
                      help = "job log file; the Snakemake log: path")
  args <- parser$parse_args()
  log_job(args$log, "demo_figures")

  use_retroseek_style()
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
    plot <- stamp_tier_note(plot, "Demo data, anonymised.")
    save_plot(name, plot, args$output,
              base_w = plot_width, base_h = plot_height, dpi = plot_dpi)
  }

  # ---- Load + anonymise -----------------------------------------------------
  all_full <- load_plot_dataframes(args$input)
  verify_required_columns(
    all_full,
    c("probe_type", "species", "virus", "probe", "label",
      "abbreviation", "max_bitscore", "query_coverage"),
    source_label = "input parquets"
  )

  # Maps built from the FULL value sets so every panel uses the same scheme.
  maps <- list(
    species      = make_label_map(all_full$species,      "Species"),
    virus        = make_label_map(all_full$virus,        "Provirus"),
    abbreviation = make_label_map(all_full$abbreviation, "Pv"),
    label        = make_label_map(all_full$label,        "Lineage")
  )
  all_full <- anonymise(all_full, maps)
  all_main <- all_full %>% dplyr::filter(probe_type == "main")


  all_counted_probe <- group_count(all_full)
  main_counted_probe <- group_count(all_main)

  pair_counts <- function(df, ax_a, ax_b) {
    df %>%
      dplyr::group_by(.data[[ax_a]], .data[[ax_b]]) %>%
      dplyr::summarise(count = sum(count), .groups = "drop")
  }

  # ---- Six README figures ---------------------------------------------------
  # 1. Sankey A: species -> probe
  emit("sankey_a.png",
       sankey_species_probe_plot(pair_counts(main_counted_probe, "species", "probe"),
                                 top_n = top_n, other_label = other_label,
                                 subset_label = "Main"))
  # 2. Waffle: provirus proportions
  emit("waffle.png",
       waffle_virus_plot(all_main, unit_hits = unit_hits, subset_label = "Main"))
  # 3. Bubble / balloon: provirus x species
  emit("balloon.png",
       balloon_virus_species_plot(main_counted_probe, subset_label = "Main"))
  # 4. Raincloud: bitscore distribution by probe
  emit("raincloud.png",
       raincloud_bitscore_plot(all_main, x_scale = x_scale, subset_label = "Main"))
  # 5. Range counts per species (bar, stacked by lineage)
  emit("bar.png",
       bar_plot(all_counted_probe, subset_label = "All probes"))

  # 6. ERV-like composition heatmap: taxon x gene, from the taxon-founded loci
  #    table (taxonomy_classify output). Taxon + gene names are public taxonomy,
  #    so nothing here needs anonymising. The table lives in a sibling dir.
  taxonomy_dir <- file.path(dirname(args$input), "taxonomy_classification")
  loci <- load_taxon_loci(taxonomy_dir)
  emit("erv_like_heatmap.png", composition_heatmap_plot(loci))

  log_ok("wrote 6 anonymised demo figures to %s", args$output)
}

if (sys.nframe() == 0L) run_main(main)
