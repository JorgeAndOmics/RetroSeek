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
source(file.path(.script_dir, "stage_plot_generator", "plots_overlap.R"))
source(file.path(.script_dir, "stage_plot_generator", "plots_ltr_interaction.R"))


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
  # (set by auto-scaling builders) and threads it to io.R::save_plot. `tier`
  # stamps a reduced-state note so each PNG is clear about which range tier it
  # shows; stage plots draw from different tiers, so it's passed per plot. Read
  # intended_dims BEFORE stamping, since `+` drops attributes.
  emit <- function(name, plot, tier = NULL) {
    dims <- attr(plot, "intended_dims")
    plot <- stamp_tier_note(plot, tier)
    save_plot(name, plot, args$output,
              dims = dims, base_w = plot_width, base_h = plot_height, dpi = plot_dpi)
  }

  # ---------- Phase 1: load --------------------------------------------------
  stage        <- load_stage_dataframes(args$input)
  counts_df    <- load_counts_table(args$input)
  coverage_df  <- load_reduction_coverage(args$input)

  # Tier notes: homology_loci is the original tier (first-reduced, non-reduced
  # track); reduction_multiplicity is the globally-reduced tier; ltr_structure
  # describes LTRdigest elements; counts are pipeline-wide.
  tier_original <- "original tier (non-reduced)"
  tier_ltr      <- "LTRdigest elements"
  tier_counts   <- "pipeline counts (all tiers)"

  # ---------- Phase 2: concordance + per-probe yield -------------------------
  log_section("Rendering concordance plots (homology ↔ LTR)")
  emit("homology_ltr_concordance.png",
       concordance_plot(stage$hits, warning_caption = warn_caption), tier_original)
  emit("probe_yield_funnel.png",
       probe_yield_plot(stage$hits, warning_caption = warn_caption), tier_original)

  # ---------- Phase 3: LTR structure ----------------------------------------
  log_section("Rendering LTR structure plots")
  emit("ltr_structural_components.png",
       ltr_structure_components_plot(stage$ltr, warning_caption = warn_caption), tier_ltr)
  emit("retrotransposon_domain_composition.png",
       domain_composition_plot(stage$ltr, warning_caption = warn_caption), tier_ltr)

  # ---------- Phase 4: refinement funnel ------------------------------------
  log_section("Rendering refinement funnel plots")
  emit("refinement_funnel_per_genome.png",
       refinement_funnel_plot(counts_df, warning_caption = warn_caption), tier_counts)
  emit("refinement_funnel_aggregate.png",
       aggregate_funnel_plot(counts_df, warning_caption = warn_caption), tier_counts)

  # ---------- Phase 5: multiplicity -----------------------------------------
  log_section("Rendering multiplicity plots")
  emit("multiplicity_m1_by_tier.png",
       multiplicity_m1_plot(stage$hits, warning_caption = warn_caption), tier_original)
  emit("multiplicity_m2_global.png",
       multiplicity_m2_plot(stage$reduced, warning_caption = warn_caption),
       "globally reduced loci")

  # ---------- Phase 6: pre-reduction overlap / redundancy --------------------
  log_section("Rendering overlap / redundancy plots")
  emit("provirus_overlap_degree.png",
       overlap_degree_plot(stage$overlap, warning_caption = warn_caption), tier_original)
  emit("provirus_reciprocal_overlap.png",
       reciprocal_fraction_plot(stage$overlap, warning_caption = warn_caption), tier_original)
  emit("provirus_reduction_fold.png",
       reduction_fold_plot(stage$hits, stage$reduced, warning_caption = warn_caption),
       "unreduced vs reduced")
  emit("provirus_coverage_before_after.png",
       coverage_before_after_plot(coverage_df, warning_caption = warn_caption),
       "unreduced vs reduced")

  # ---------- Phase 7: LTR feature interactions ------------------------------
  log_section("Rendering LTR-interaction plots")
  emit("ltr_distance_to_retro.png",
       distance_to_retro_plot(stage$ltr_int, warning_caption = warn_caption), tier_original)
  emit("ltr_position_within_provirus.png",
       position_within_provirus_plot(stage$ltr_int, warning_caption = warn_caption), tier_original)
  emit("ltr_feature_breakdown.png",
       ltr_feature_breakdown_plot(stage$ltr_int, warning_caption = warn_caption), tier_original)
  emit("ltr_strand_concordance.png",
       strand_concordance_plot(stage$ltr_int, warning_caption = warn_caption), tier_original)
  emit("ltr_probe_domain_overlap.png",
       probe_domain_heatmap(stage$probe_domain, warning_caption = warn_caption), tier_original)
  emit("ltr_retro_length_vs_hits.png",
       retro_length_vs_hits_plot(stage$ltr, warning_caption = warn_caption), tier_ltr)

  log_section(sprintf("Done — wrote 18 PNGs to %s", args$output))
}


# ----------------------------------------------------------------------------
# Entry-point guard — only fire main() under `Rscript stage_plot_generator.R`.
# ----------------------------------------------------------------------------
if (sys.nframe() == 0L) main()
