# =============================================================================
# plot2sort.R — orchestrator
# =============================================================================
# Builds RetroSeek's global PNG plots from the per-genome `final_loci` tables
# that `ranges_analysis.R` writes to `data/tables/ranges_analysis/`. Heavy
# lifting is in `workflow/scripts/plot2sort/*.R`; this file argument-parses,
# sources the modules, prepares the joined per-probe-type frames, emits each PNG.
#
# Phases:
#   1. Load + validate the per-genome parquets (verifies probe_type +
#      query_coverage columns are present).
#   2. Attach human-readable species names from the YAML species map.
#   3. Aggregate per-probe-type frames + bitscore quartile summaries.
#   4. Render distribution plots (density, raincloud, query coverage).
#   5. Render categorical plots (bar, balloon, heatmap, waffle).
#   6. Render alluvial / Sankey plots (species×probe, label×probe, species×label).
#
# `testthat` sources this file. The bottom-of-file `if (sys.nframe() == 0L) main()`
# guard prevents the CLI block from firing during sourcing — testthat then sees
# every helper + builder via the sourced sub-modules.

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


# ----------------------------------------------------------------------------
# Locate sibling scripts + source modules
# ----------------------------------------------------------------------------
.resolve_script_dir <- function() {
  # Walk the call stack: the most recent source() frame carries `ofile`.
  # Required for testthat (which source()s us several frames deep) and also
  # works under Rscript (no frames have ofile, so we fall through).
  for (i in rev(seq_len(sys.nframe()))) {
    fr <- tryCatch(sys.frame(i), error = function(e) NULL)
    if (is.null(fr)) next
    ofile <- tryCatch(fr$ofile, error = function(e) NULL)
    if (!is.null(ofile)) {
      return(dirname(normalizePath(ofile, mustWork = FALSE)))
    }
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


# ----------------------------------------------------------------------------
# Pipeline instrumentation — same idiom as ranges_analysis.R / hotspot_detector.R
# ----------------------------------------------------------------------------
.t0 <- Sys.time()
log_section <- function(name) {
  elapsed <- as.numeric(difftime(Sys.time(), .t0, units = "secs"))
  message(sprintf("[%6.2fs] > %s", elapsed, name))
}


# ----------------------------------------------------------------------------
# main()
# ----------------------------------------------------------------------------
main <- function() {
  parser <- ArgumentParser(
    description = "Generate global plots from RetroSeek per-genome Parquet results"
  )
  parser$add_argument("--input",  required = TRUE,
                      help = "Directory with the per-genome ranges-analysis tables; reads {genome}.final_loci.parquet (carries a probe_type column: main | accessory).")
  parser$add_argument("--output", required = TRUE,
                      help = "Directory to save output plots")
  parser$add_argument("--config", required = TRUE,
                      help = "YAML config file with plot parameters")
  args <- parser$parse_args()

  cfg <- yaml::read_yaml(args$config)
  plot_dpi    <- cfg$plots$dpi    %||% 300
  plot_height <- cfg$plots$height %||% 12
  plot_width  <- cfg$plots$width  %||% 15
  x_scale     <- cfg$plots$bitscore_x_scale %||% "linear"
  top_n       <- cfg$plots$sankey_top_n            # NULL = show all
  other_label <- cfg$plots$sankey_other_label %||% "Other"
  unit_hits   <- cfg$plots$waffle_unit_hits        # NULL → auto-derive

  dir.create(args$output, showWarnings = FALSE, recursive = TRUE)
  log_section(sprintf("RetroSeek plot generation (output: %s)", args$output))

  # Multi-value aggregation warning: when virus/label use list/concatenate,
  # plot counts are inflated by entry explosion. Log it once and stamp every
  # emitted plot so the caveat travels with the PNG artifact.
  warn_caption <- aggregation_warning(cfg)
  if (!is.null(warn_caption)) {
    warning(warn_caption, call. = FALSE)
    log_section(paste0("  ", warn_caption))
  }

  # Local save wrapper that:
  #   1) captures output_dir + base canvas + dpi from config,
  #   2) extracts the per-builder `intended_dims` attribute (set by
  #      auto-scaling builders) — read BEFORE stamping, since `+` drops attrs,
  #   3) stamps the multi-value aggregation warning caption when active,
  #   4) passes the captured dims to io.R::save_plot.
  emit <- function(name, plot) {
    dims <- attr(plot, "intended_dims")
    plot <- stamp_warning_caption(plot, warn_caption)
    save_plot(name, plot, args$output,
              dims = dims, base_w = plot_width, base_h = plot_height,
              dpi  = plot_dpi)
  }

  # ---------- Phase 1: load + validate ---------------------------------------
  all.full <- load_plot_dataframes(args$input)
  log_section("Verifying required columns")
  verify_required_columns(
    all.full,
    c("probe_type", "species", "virus", "probe", "label",
      "abbreviation", "max_bitscore", "query_coverage"),
    source_label = "input parquets"
  )

  all.main      <- all.full %>% dplyr::filter(probe_type == "main")
  all.accessory <- all.full %>% dplyr::filter(probe_type == "accessory")
  log_section(sprintf(
    "Split by probe_type: main = %d rows, accessory = %d rows",
    nrow(all.main), nrow(all.accessory)
  ))

  # ---------- Phase 2: attach human-readable species names -------------------
  log_section("Attaching species names from config map")
  species_map   <- cfg$species
  all.main      <- attach_species_name(all.main,      species_map)
  all.accessory <- attach_species_name(all.accessory, species_map)
  all.full      <- attach_species_name(all.full,      species_map)

  # ---------- Phase 3: aggregations + bitscore quartiles ---------------------
  log_section("Computing bitscore quartiles + per-probe counts")
  bit.main      <- if (nrow(all.main)      > 0L) q_stats(all.main)
                   else list(q1 = NA, median = NA, q3 = NA, mean = NA)
  bit.accessory <- if (nrow(all.accessory) > 0L) q_stats(all.accessory)
                   else list(q1 = NA, median = NA, q3 = NA, mean = NA)

  all.counted_probe       <- group_count(all.full)
  main.counted_probe      <- group_count(all.main)
  accessory.counted_probe <- group_count(all.accessory)

  # Sankey input frames: collapse counted_probe to a 2-axis count for each pair.
  pair_counts <- function(df, ax_a, ax_b) {
    df %>%
      dplyr::group_by(.data[[ax_a]], .data[[ax_b]]) %>%
      dplyr::summarise(count = sum(count), .groups = "drop")
  }
  data.main.sankey.species.probe      <- pair_counts(main.counted_probe,      "species", "probe")
  data.accessory.sankey.species.probe <- pair_counts(accessory.counted_probe, "species", "probe")
  data.main.sankey.label.probe        <- pair_counts(main.counted_probe,      "label",   "probe")
  data.accessory.sankey.label.probe   <- pair_counts(accessory.counted_probe, "label",   "probe")
  data.main.sankey.species.label      <- pair_counts(main.counted_probe,      "species", "label")
  data.accessory.sankey.species.label <- pair_counts(accessory.counted_probe, "species", "label")

  # ---------- Phase 4: distribution plots ------------------------------------
  log_section("Rendering distribution plots (density, raincloud, query coverage)")
  emit("main_density.png",
       density_bitscore_plot(all.main,
                             bit.main$q1, bit.main$median, bit.main$q3,
                             x_scale = x_scale, subset_label = "Main"))
  emit("accessory_density.png",
       density_bitscore_plot(all.accessory,
                             bit.accessory$q1, bit.accessory$median, bit.accessory$q3,
                             x_scale = x_scale, subset_label = "Accessory"))
  emit("main_raincloud.png",
       raincloud_bitscore_plot(all.main, x_scale = x_scale,
                               subset_label = "Main"))
  emit("accessory_raincloud.png",
       raincloud_bitscore_plot(all.accessory, x_scale = x_scale,
                               subset_label = "Accessory"))
  emit("main_query_coverage.png",
       query_coverage_plot(all.main,      subset_label = "Main"))
  emit("accessory_query_coverage.png",
       query_coverage_plot(all.accessory, subset_label = "Accessory"))

  # ---------- Phase 5: categorical plots -------------------------------------
  log_section("Rendering categorical plots (bar, balloon, heatmap, waffle)")
  emit("full_bar.png",      bar_plot(all.counted_probe,       subset_label = "All probes"))
  emit("main_bar.png",      bar_plot(main.counted_probe,      subset_label = "Main"))
  emit("accessory_bar.png", bar_plot(accessory.counted_probe, subset_label = "Accessory"))
  emit("main_balloon.png",
       balloon_virus_species_plot(main.counted_probe,      subset_label = "Main"))
  emit("accessory_balloon.png",
       balloon_virus_species_plot(accessory.counted_probe, subset_label = "Accessory"))
  emit("main_heatmap_probe_species.png",
       heatmap_probe_species_plot(all.main,      subset_label = "Main"))
  emit("accessory_heatmap_probe_species.png",
       heatmap_probe_species_plot(all.accessory, subset_label = "Accessory"))
  emit("main_waffle_virus.png",
       waffle_virus_plot(all.main, unit_hits = unit_hits,
                         subset_label = "Main"))
  emit("accessory_waffle_virus.png",
       waffle_virus_plot(all.accessory, unit_hits = unit_hits,
                         subset_label = "Accessory"))

  # ---------- Phase 6: alluvial / Sankey -------------------------------------
  log_section("Rendering alluvial plots (sankey × {a, b, c} × {main, accessory})")
  emit("main_sankey_a.png",
       sankey_species_probe_plot(data.main.sankey.species.probe,
                                 top_n = top_n, other_label = other_label,
                                 subset_label = "Main"))
  emit("main_sankey_b.png",
       sankey_label_probe_plot(data.main.sankey.label.probe,
                               top_n = top_n, other_label = other_label,
                               subset_label = "Main"))
  emit("main_sankey_c.png",
       sankey_species_label_plot(data.main.sankey.species.label,
                                 top_n = top_n, other_label = other_label,
                                 subset_label = "Main"))
  emit("accessory_sankey_a.png",
       sankey_species_probe_plot(data.accessory.sankey.species.probe,
                                 top_n = top_n, other_label = other_label,
                                 subset_label = "Accessory"))
  emit("accessory_sankey_b.png",
       sankey_label_probe_plot(data.accessory.sankey.label.probe,
                               top_n = top_n, other_label = other_label,
                               subset_label = "Accessory"))
  emit("accessory_sankey_c.png",
       sankey_species_label_plot(data.accessory.sankey.species.label,
                                 top_n = top_n, other_label = other_label,
                                 subset_label = "Accessory"))

  log_section(sprintf("Done — wrote 21 PNGs to %s", args$output))
}


# ----------------------------------------------------------------------------
# Entry-point guard
# ----------------------------------------------------------------------------
# Only invoke main() under `Rscript plot2sort.R …`. testthat sources this file
# inside test functions where sys.nframe() > 0, so unit tests get the helpers
# and builders without firing the CLI.
if (sys.nframe() == 0L) main()
