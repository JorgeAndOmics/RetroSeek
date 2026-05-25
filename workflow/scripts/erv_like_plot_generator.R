# =============================================================================
# erv_like_plot_generator.R — orchestrator
# =============================================================================
# Builds RetroSeek's ERV-like panel: plots describing the composite ERV-like
# candidates assembled by range_analysis/erv_assembly.R (composition,
# completeness, observed gene order, canonical-vs-rearranged, length, member
# count, inter-probe gap). Input is the per-genome `erv_like_loci` +
# `erv_like_members` parquet tables (plus the `counts` table for the dropped-
# non-canonical tally) that ranges_analysis.R writes to
# `data/tables/ranges_analysis/`.
#
# Shared infrastructure (theme, add_titles, stamp_tier_note, auto_dims,
# empty_plot, save_plot, the aggregation-warning helper) is reused from
# plot2sort/*.R — not duplicated — so the panel matches the existing plots.
#
# Phases:
#   1. Load the erv_like_loci + erv_like_members tables + dropped-noncanon count.
#   2. Composition & completeness (probe combinations, completeness, by-virus,
#      probe × virus heatmap).
#   3. Order & structure (observed gene order, canonical vs rearranged, length,
#      member count, inter-probe gap vs join cutoff).
#
# `testthat` sources this file; the `if (sys.nframe() == 0L) main()` guard keeps
# the CLI block from firing during sourcing.

suppressMessages({
  library(argparse)     # Command-line argument parser
  library(arrow)        # Parquet I/O
  library(tidyverse)    # Data manipulation and visualisation
  library(yaml)         # YAML config
  library(ggsci)        # Scientific colour palettes (planet express)
  library(scales)       # Axis labellers
})


# ----------------------------------------------------------------------------
# Locate sibling scripts + source modules
# ----------------------------------------------------------------------------
.resolve_script_dir <- function() {
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
# Reused from plot2sort: helpers (theme, add_titles, stamp_tier_note, auto_dims,
# empty_plot, order_by_count, collapse_long_tail, aggregation_warning) + io
# (save_plot). Then the erv-like-specific modules.
source(file.path(.script_dir, "plot2sort", "helpers.R"))
source(file.path(.script_dir, "plot2sort", "io.R"))
source(file.path(.script_dir, "erv_like_plot_generator", "io.R"))
source(file.path(.script_dir, "erv_like_plot_generator", "plots_composition.R"))
source(file.path(.script_dir, "erv_like_plot_generator", "plots_order_structure.R"))


# ----------------------------------------------------------------------------
# Pipeline instrumentation — same idiom as the other plot generators.
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
    description = "Generate RetroSeek ERV-like panel (composite candidate plots)"
  )
  parser$add_argument("--input", required = TRUE,
                      help = paste("Directory with per-genome ranges-analysis parquet",
                                   "tables (erv_like_loci / erv_like_members / counts)."))
  parser$add_argument("--output", required = TRUE,
                      help = "Directory to save output plots.")
  parser$add_argument("--config", required = TRUE,
                      help = "YAML config file with plot parameters.")
  args <- parser$parse_args()

  cfg <- yaml::read_yaml(args$config)
  plot_dpi    <- cfg$plots$dpi    %||% 300
  plot_height <- cfg$plots$height %||% 12
  plot_width  <- cfg$plots$width  %||% 15
  top_n       <- cfg$plots$sankey_top_n              # NULL = show all
  other_label <- cfg$plots$sankey_other_label %||% "Other"
  sep         <- cfg$parameters$aggregation$concat_separator %||% "; "
  max_join    <- cfg$parameters$erv_like$max_join_distance   %||% 1500

  dir.create(args$output, showWarnings = FALSE, recursive = TRUE)
  log_section(sprintf("RetroSeek erv-like plot generation (output: %s)", args$output))

  warn_caption <- aggregation_warning(cfg)
  if (!is.null(warn_caption)) {
    warning(warn_caption, call. = FALSE)
    log_section(paste0("  ", warn_caption))
  }

  # Every erv-like plot describes the erv_like candidate tier; stamp it on each.
  # Read intended_dims BEFORE stamping, since `+` drops attributes.
  emit <- function(name, plot) {
    dims <- attr(plot, "intended_dims")
    plot <- stamp_tier_note(plot, "erv-like candidates")
    plot <- stamp_warning_caption(plot, warn_caption)
    save_plot(name, plot, args$output,
              dims = dims, base_w = plot_width, base_h = plot_height,
              dpi  = plot_dpi)
  }

  # ---------- Phase 1: load --------------------------------------------------
  loci    <- load_erv_like_loci(args$input)
  members <- load_erv_like_members(args$input)
  dropped <- load_dropped_noncanonical(args$input)
  log_section(sprintf("Loaded %d candidates, %d members, %d dropped (non-canonical)",
                      nrow(loci), nrow(members), dropped))

  # ---------- Phase 2: composition & completeness ----------------------------
  log_section("Rendering composition & completeness plots")
  emit("erv_like_probe_combinations.png",
       probe_combination_plot(loci, top_n = top_n, other_label = other_label))
  emit("erv_like_completeness.png", completeness_plot(loci))
  emit("erv_like_full_partial_by_virus.png",
       full_partial_by_virus_plot(loci, sep = sep))
  emit("erv_like_composition_heatmap.png",
       composition_heatmap_plot(loci, sep = sep))

  # ---------- Phase 3: order & structure -------------------------------------
  log_section("Rendering order & structure plots")
  emit("erv_like_gene_order.png",
       gene_order_plot(loci, top_n = top_n, other_label = other_label))
  emit("erv_like_canonical.png",
       canonical_plot(loci, dropped_noncanonical = dropped))
  emit("erv_like_length_distribution.png", candidate_length_plot(loci))
  emit("erv_like_n_loci.png", n_loci_plot(loci))
  emit("erv_like_interprobe_gap.png",
       interprobe_gap_plot(members, max_join_distance = max_join))

  log_section(sprintf("Done — wrote 9 PNGs to %s", args$output))
}


# ----------------------------------------------------------------------------
# Entry-point guard — only fire main() under `Rscript erv_like_plot_generator.R`.
# ----------------------------------------------------------------------------
if (sys.nframe() == 0L) main()
