# =============================================================================
# stage_plot_generator.R - orchestrator
# =============================================================================
# Builds RetroSeek's integration stage PDF: the homology + LTRharvest/LTRdigest
# integration step that the homology pages from plot2sort.R do not capture.
# Input is the per-genome ranges-analysis parquet tables that ranges_analysis.R
# writes to `data/tables/ranges_analysis/` (homology_loci / ltr_structure /
# reduction_multiplicity / counts / overlap / interaction tables).
#
# Output is ONE PDF (results/plots/ranges/integration/integration.pdf): a key
# page, then the refinement funnel, homology against the LTR elements, the
# elements themselves, and what the reductions collapse. Each page's subtitle
# ends with the range tier it shows.
#
# Shared infrastructure (style, theme, add_titles, empty_plot, the
# aggregation-warning helper, species order) is reused from plot2sort/*.R.
#
# `testthat` sources this file; the `if (sys.nframe() == 0L) main()` guard
# keeps the CLI block from firing during sourcing.

suppressMessages({
  library(argparse)     # Command-line argument parser
  library(arrow)        # Parquet I/O
  library(tidyverse)    # Data manipulation and visualisation
  library(yaml)         # YAML config + manifests
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
# Reused from plot2sort: style (palette, theme, stage PDFs), helpers
# (add_titles, empty_plot, order_by_count, collapse_long_tail,
# aggregation_warning) and tree_axis (species order).
source(file.path(.script_dir, "utils", "log.R"))  # line contract, run_main (ADR-021)
source(file.path(.script_dir, "plot2sort", "style.R"))
source(file.path(.script_dir, "plot2sort", "helpers.R"))
source(file.path(.script_dir, "plot2sort", "tree_axis.R"))
# Stage-specific modules.
source(file.path(.script_dir, "stage_plot_generator", "io.R"))
source(file.path(.script_dir, "stage_plot_generator", "plots_concordance.R"))
source(file.path(.script_dir, "stage_plot_generator", "plots_structure.R"))
source(file.path(.script_dir, "stage_plot_generator", "plots_funnel.R"))
source(file.path(.script_dir, "stage_plot_generator", "plots_multiplicity.R"))
source(file.path(.script_dir, "stage_plot_generator", "plots_overlap.R"))
source(file.path(.script_dir, "stage_plot_generator", "plots_ltr_interaction.R"))


# ----------------------------------------------------------------------------
# main()
# ----------------------------------------------------------------------------
main <- function() {
  parser <- ArgumentParser(
    description = "Generate the RetroSeek integration stage PDF (homology + LTR integration)"
  )
  parser$add_argument("--input", required = TRUE,
                      help = paste("Directory with per-genome ranges-analysis parquet",
                                   "tables (homology_loci / ltr_structure /",
                                   "reduction_multiplicity / counts)."))
  parser$add_argument("--out_pdf", required = TRUE,
                      help = "The stage PDF: a key page, then the integration pages.")
  parser$add_argument("--config", required = TRUE,
                      help = "YAML config file with plot parameters.")
  parser$add_argument("--species_tree_dir", default = "",
                      help = "species_tree_layout.py output: the host tree coordinates")
  parser$add_argument("--log", default = NULL,
                      help = "job log file; the Snakemake log: path")
  args <- parser$parse_args()
  log_job(args$log, "stage_plot_generator")

  use_retroseek_style()
  cfg <- yaml::read_yaml(args$config)
  log_section(sprintf("RetroSeek integration plots (output: %s)", args$out_pdf))

  # Multi-value aggregation warning: when virus/label use list/concatenate,
  # plot counts are inflated by entry explosion. Log it and stamp every page.
  warn <- aggregation_warning(cfg)
  if (!is.null(warn)) warning(warn, call. = FALSE)

  # ---------- Phase 1: load --------------------------------------------------
  stage        <- load_stage_dataframes(args$input)
  counts_df    <- load_counts_table(args$input)
  coverage_df  <- load_reduction_coverage(args$input)
  counts_df$genome <- display_species(counts_df$genome, cfg$species)
  ctx <- panel_ctx(cfg, args$species_tree_dir %||% "")
  # One probe palette for the whole stage, so a probe has one colour on every page.
  probe_cols <- probe_colours(stage$hits$probe)

  # Tier notes, appended to each subtitle: homology_loci is the original tier
  # (first-reduced, not globally reduced); reduction_multiplicity is the globally
  # reduced tier; ltr_structure describes LTRdigest elements; counts span tiers.
  original <- "Original tier, before the global reduction."
  elements <- "LTRdigest elements."
  counts   <- "Pipeline counts, all tiers."
  both     <- "Before and after the global reduction."
  page <- function(plot, tier) stamp_tier_note(plot, tier)

  # ---------- Phase 2: pages -------------------------------------------------
  log_section("Building the pages")
  pages <- list(
    # The refinement funnel.
    page(refinement_funnel_plot(counts_df, warning_caption = warn, ctx = ctx), counts),
    page(aggregate_funnel_plot(counts_df, warning_caption = warn), counts),
    # Homology loci against the LTR elements.
    page(concordance_plot(stage$hits, warning_caption = warn), original),
    page(probe_yield_plot(stage$hits, warning_caption = warn), original),
    page(ltr_feature_breakdown_plot(stage$ltr_int, warning_caption = warn), original),
    page(distance_to_retro_plot(stage$ltr_int, warning_caption = warn), original),
    page(position_within_provirus_plot(stage$ltr_int, warning_caption = warn), original),
    page(strand_concordance_plot(stage$ltr_int, warning_caption = warn), original),
    page(probe_domain_heatmap(stage$probe_domain, warning_caption = warn), original),
    # The LTR elements themselves.
    page(ltr_structure_components_plot(stage$ltr, warning_caption = warn), elements),
    page(domain_composition_plot(stage$ltr, warning_caption = warn), elements),
    page(retro_length_vs_hits_plot(stage$ltr, warning_caption = warn), elements),
    # What the reductions collapse.
    page(multiplicity_m1_plot(stage$hits, warning_caption = warn), original),
    page(multiplicity_m2_plot(stage$reduced, warning_caption = warn),
         "Globally reduced loci."),
    page(overlap_degree_plot(stage$overlap, warning_caption = warn,
                             colours = probe_cols), original),
    page(reciprocal_fraction_plot(stage$overlap, warning_caption = warn), original),
    page(reduction_fold_plot(stage$hits, stage$reduced, warning_caption = warn), both),
    page(coverage_before_after_plot(coverage_df, warning_caption = warn), both)
  )

  key <- key_page(
    "Homology meets LTR structure",
    paste("How the tBLASTn homology loci meet the LTR retrotransposons LTRharvest and",
          "LTRdigest find: where the loci sit relative to the elements, what the",
          "elements carry, and how the reductions collapse overlapping loci. The",
          "last line of each subtitle names the range tier the page shows."),
    colours = c(`LTR-flanked, or overlapping an LTR element` = .TIER_COLOUR[["ltr-flanked"]],
                `Before reduction, or not selected` = .GREY_MID,
                `After reduction, or the flagged subset` = .DATA_COLOUR),
    pages = page_titles(pages))
  n_rows <- max(length(ctx$species_order), length(unique(counts_df$genome)))
  save_stage_pdf(c(list(key), pages), args$out_pdf,
                 height = page_height_for(n_rows, per_species = cfg$plots$per_stratum %||% 0.18))
  log_ok("wrote %s, %s pages", basename(args$out_pdf),
         format(length(pages) + 1L, big.mark = ","))
}


# ----------------------------------------------------------------------------
# Entry-point guard - only fire main() under `Rscript stage_plot_generator.R`.
# run_main() logs how the job ended (ADR-021).
# ----------------------------------------------------------------------------
if (sys.nframe() == 0L) run_main(main)
