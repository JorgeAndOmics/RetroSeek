# =============================================================================
# stage_plot_generator.R — orchestrator
# =============================================================================
# Builds RetroSeek's "middle stage" PNG plots — the homology + LTRharvest/
# LTRdigest integration step that the final-tier plots from plot2sort.R do not
# capture. Input is the per-genome ranges-analysis parquet tables that
# ranges_analysis.R writes to `data/tables/ranges_analysis/` (homology_loci /
# ltr_structure / reduction_multiplicity / counts).
#
# Shared infrastructure (theme, add_titles, auto_dims, empty_plot, save_plot,
# the aggregation-warning helper) is reused from plot2sort/*.R — not
# duplicated — so the new plots match the existing panel's format and
# auto-scaling exactly.
#
# Phases:
#   1. Load the three stage-dataframe parquet types + the counts table.
#   2. Concordance plots (homology ↔ LTR spatial relationship, per-probe yield).
#   3. LTR structure plots (structural completeness, domain composition).
#   4. Refinement funnel plots (per-genome + aggregate cohort).
#   5. Multiplicity plots (M1 by tier, M2 global).
#
# `testthat` sources this file; the `if (sys.nframe() == 0L) main()` guard
# keeps the CLI block from firing during sourcing.

suppressMessages({
  library(argparse)     # Command-line argument parser
  library(arrow)        # Parquet I/O
  library(tidyverse)    # Data manipulation and visualisation
  library(yaml)         # YAML config + manifests
  library(ggsci)        # Scientific colour palettes (planet express)
  library(scales)       # Axis labellers
})


# ----------------------------------------------------------------------------
# Locate sibling scripts + source modules
# ----------------------------------------------------------------------------
.resolve_script_dir <- function() {
  # Walk the call stack: the most recent source() frame carries `ofile`.
  # Works under testthat (source()d several frames deep) and Rscript (falls
  # through to the --file= argument).
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
# Reused from plot2sort: helpers (theme, add_titles, auto_dims, empty_plot,
# order_by_count, collapse_long_tail, futurama_unlimited_palette,
# aggregation_warning) + io (save_plot, verify_required_columns).
source(file.path(.script_dir, "plot2sort", "helpers.R"))
source(file.path(.script_dir, "plot2sort", "io.R"))
# Stage-specific modules.
source(file.path(.script_dir, "stage_plot_generator", "io.R"))
source(file.path(.script_dir, "stage_plot_generator", "plots_concordance.R"))
source(file.path(.script_dir, "stage_plot_generator", "plots_structure.R"))
source(file.path(.script_dir, "stage_plot_generator", "plots_funnel.R"))
source(file.path(.script_dir, "stage_plot_generator", "plots_multiplicity.R"))


# ----------------------------------------------------------------------------
# Pipeline instrumentation — same idiom as ranges_analysis.R / plot2sort.R
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
    description = "Generate RetroSeek middle-stage plots (homology + LTR integration)"
  )
  parser$add_argument("--input", required = TRUE,
                      help = paste("Directory with per-genome ranges-analysis parquet",
                                   "tables (homology_loci / ltr_structure /",
                                   "reduction_multiplicity / counts)."))
  parser$add_argument("--output", required = TRUE,
                      help = "Directory to save output plots.")
  parser$add_argument("--config", required = TRUE,
                      help = "YAML config file with plot parameters.")
  args <- parser$parse_args()

  cfg <- yaml::read_yaml(args$config)
  plot_dpi    <- cfg$plots$dpi    %||% 300
  plot_height <- cfg$plots$height %||% 12
  plot_width  <- cfg$plots$width  %||% 15

  dir.create(args$output, showWarnings = FALSE, recursive = TRUE)
  log_section(sprintf("RetroSeek stage-plot generation (output: %s)", args$output))

  # Multi-value aggregation warning: when virus/label use list/concatenate,
  # plot counts are inflated by entry explosion. Log it and stamp every plot.
  warn_caption <- aggregation_warning(cfg)
  if (!is.null(warn_caption)) {
    warning(warn_caption, call. = FALSE)
    log_section(paste0("  ", warn_caption))
  }

  # Local save wrapper — extracts the per-builder `intended_dims` attribute
  # (set by auto-scaling builders) and threads it to io.R::save_plot.
  emit <- function(name, plot) {
    save_plot(name, plot, args$output,
              dims = attr(plot, "intended_dims"),
              base_w = plot_width, base_h = plot_height, dpi = plot_dpi)
  }

  # ---------- Phase 1: load --------------------------------------------------
  stage     <- load_stage_dataframes(args$input)
  counts_df <- load_counts_table(args$input)

  # ---------- Phase 2: concordance + per-probe yield -------------------------
  log_section("Rendering concordance plots (homology ↔ LTR)")
  emit("homology_ltr_concordance.png",
       concordance_plot(stage$hits, warning_caption = warn_caption))
  emit("probe_yield_funnel.png",
       probe_yield_plot(stage$hits, warning_caption = warn_caption))

  # ---------- Phase 3: LTR structure ----------------------------------------
  log_section("Rendering LTR structure plots")
  emit("ltr_structural_components.png",
       ltr_structure_components_plot(stage$ltr, warning_caption = warn_caption))
  emit("retrotransposon_domain_composition.png",
       domain_composition_plot(stage$ltr, warning_caption = warn_caption))

  # ---------- Phase 4: refinement funnel ------------------------------------
  log_section("Rendering refinement funnel plots")
  emit("refinement_funnel_per_genome.png",
       refinement_funnel_plot(counts_df, warning_caption = warn_caption))
  emit("refinement_funnel_aggregate.png",
       aggregate_funnel_plot(counts_df, warning_caption = warn_caption))

  # ---------- Phase 5: multiplicity -----------------------------------------
  log_section("Rendering multiplicity plots")
  emit("multiplicity_m1_by_tier.png",
       multiplicity_m1_plot(stage$hits, warning_caption = warn_caption))
  emit("multiplicity_m2_global.png",
       multiplicity_m2_plot(stage$reduced, warning_caption = warn_caption))

  log_section(sprintf("Done — wrote 8 PNGs to %s", args$output))
}


# ----------------------------------------------------------------------------
# Entry-point guard — only fire main() under `Rscript stage_plot_generator.R`.
# ----------------------------------------------------------------------------
if (sys.nframe() == 0L) main()
