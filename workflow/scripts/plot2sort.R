# =============================================================================
# plot2sort.R - orchestrator
# =============================================================================
# Builds RetroSeek's homology stage PDF from the per-genome `final_loci` tables
# that `ranges_analysis.R` writes to `data/tables/ranges_analysis/`. Heavy
# lifting is in `workflow/scripts/plot2sort/*.R`; this file argument-parses,
# sources the modules, prepares the joined per-probe-type frames and assembles
# the pages.
#
# Output is ONE PDF (results/plots/ranges/homology/homology.pdf): a key page,
# two pages over every probe (ranges per host by lineage and by virus), then
# the same ten pages for the main probe set and for the accessory set:
# bitscore density and raincloud, query coverage, ranges per host, per virus
# and host, per probe and host, the virus waffle and three alluvia.
#
# Phases:
#   1. Load + validate the per-genome parquets (verifies probe_type +
#      query_coverage columns are present).
#   2. Readable species names from the YAML species map.
#   3. Aggregate per-probe-type frames + bitscore quartile summaries.
#   4. Build the pages and write the PDF.
#
# `testthat` sources this file. The bottom-of-file `if (sys.nframe() == 0L) main()`
# guard prevents the CLI block from firing during sourcing - testthat then sees
# every helper + builder via the sourced sub-modules.

suppressMessages({
  library(argparse)     # Command-line argument parser
  library(arrow)        # Parquet I/O
  library(tidyverse)    # Data manipulation and visualisation
  library(yaml)         # YAML config
  library(ggalluvial)   # Sankey / alluvial geoms
  library(ggdist)       # Raincloud, halfeye and dots geoms
  library(scales)       # Axis labellers
})


# ----------------------------------------------------------------------------
# Locate sibling scripts + source modules
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
source(file.path(.script_dir, "plot2sort", "tree_axis.R"))  # species rows, host tree
source(file.path(.script_dir, "plot2sort", "plots_distribution.R"))
source(file.path(.script_dir, "plot2sort", "plots_categorical.R"))
source(file.path(.script_dir, "plot2sort", "plots_sankey.R"))


# ----------------------------------------------------------------------------
# Entry point
# ----------------------------------------------------------------------------
main <- function() {
  parser <- ArgumentParser(
    description =
      "Generate the homology stage PDF from RetroSeek per-genome Parquet results"
  )
  parser$add_argument("--input",  required = TRUE,
                      help = paste(
                        "Directory with the per-genome ranges-analysis tables; reads",
                        "{genome}.final_loci.parquet (carries a probe_type column:",
                        "main | accessory)."
                      ))
  parser$add_argument("--out_pdf", required = TRUE,
                      help = "The stage PDF: a key page, then the homology pages.")
  parser$add_argument("--config", required = TRUE,
                      help = "YAML config file with plot parameters")
  parser$add_argument("--species_tree_dir", default = "",
                      help = "species_tree_layout.py output: the host tree coordinates")
  parser$add_argument("--log", default = NULL,
                      help = "job log file; the Snakemake log: path")
  args <- parser$parse_args()
  log_job(args$log, "plot_generator")

  use_retroseek_style()
  cfg <- yaml::read_yaml(args$config)
  x_scale     <- cfg$plots$bitscore_x_scale %||% "linear"
  top_n       <- cfg$plots$sankey_top_n            # NULL = show all
  other_label <- cfg$plots$sankey_other_label %||% "Other"
  unit_hits   <- cfg$plots$waffle_unit_hits        # NULL means derive it automatically
  log_section(sprintf("RetroSeek homology plots (output: %s)", args$out_pdf))

  # Multi-value aggregation warning: when virus/label use list/concatenate,
  # plot counts are inflated by entry explosion. Log it once and stamp every
  # page so the caveat travels with it.
  warn_caption <- aggregation_warning(cfg)
  if (!is.null(warn_caption)) warning(warn_caption, call. = FALSE)

  # ---------- Phase 1: load + validate ---------------------------------------
  all.full <- load_plot_dataframes(args$input)
  log_section("Verifying required columns")
  verify_required_columns(
    all.full,
    c("probe_type", "species", "virus", "probe", "label",
      "abbreviation", "max_bitscore", "query_coverage"),
    source_label = "input parquets"
  )

  # ---------- Phase 2: readable species names --------------------------------
  all.full$species <- display_species(all.full$species, cfg$species)
  all.main      <- all.full %>% dplyr::filter(probe_type == "main")
  all.accessory <- all.full %>% dplyr::filter(probe_type == "accessory")
  log_section(sprintf(
    "Split by probe_type: main = %d rows, accessory = %d rows",
    nrow(all.main), nrow(all.accessory)
  ))

  # ---------- Phase 3: aggregations + bitscore quartiles ---------------------
  log_section("Computing bitscore quartiles + per-probe counts")
  quartiles <- function(df) {
    if (nrow(df) > 0L) q_stats(df) else list(q1 = NA, median = NA, q3 = NA, mean = NA)
  }
  # Sankey input frames: collapse counted_probe to a 2-axis count for each pair.
  pair_counts <- function(df, ax_a, ax_b) {
    df %>%
      dplyr::group_by(.data[[ax_a]], .data[[ax_b]]) %>%
      dplyr::summarise(count = sum(count), .groups = "drop")
  }
  ctx <- panel_ctx(cfg, args$species_tree_dir %||% "")
  # One probe palette and one virus palette for the whole stage, so a probe or a
  # virus has one colour on every page. Main probes first: they take the palette's
  # own colours, the accessory probes what follows.
  probe_cols <- category_colours(c(sort(unique(all.main$probe)),
                                   sort(unique(all.accessory$probe))))
  all_counted <- group_count(all.full)
  virus_cols <- virus_colours(all_counted$virus, all_counted$label, all_counted$count)

  # ---------- Phase 4: pages -------------------------------------------------
  # The same ten pages for each probe set; `lead` opens each subtitle.
  set_pages <- function(ranges, lead) {
    counted <- group_count(ranges)
    q <- quartiles(ranges)
    list(
      density_bitscore_plot(ranges, q$q1, q$median, q$q3, x_scale = x_scale,
                            subset_label = lead, colours = probe_cols),
      raincloud_bitscore_plot(ranges, x_scale = x_scale, subset_label = lead,
                              colours = probe_cols),
      query_coverage_plot(ranges, subset_label = lead, colours = probe_cols),
      bar_plot(counted, subset_label = lead, ctx = ctx),
      balloon_virus_species_plot(counted, subset_label = lead, ctx = ctx),
      heatmap_probe_species_plot(ranges, subset_label = lead, ctx = ctx),
      waffle_virus_plot(ranges, unit_hits = unit_hits, subset_label = lead,
                        colours = virus_cols),
      sankey_species_probe_plot(pair_counts(counted, "species", "probe"),
                                top_n = top_n, other_label = other_label,
                                subset_label = lead, ctx = ctx, colours = probe_cols),
      sankey_label_probe_plot(pair_counts(counted, "label", "probe"),
                              top_n = top_n, other_label = other_label,
                              subset_label = lead, ctx = ctx),
      sankey_species_label_plot(pair_counts(counted, "species", "label"),
                                top_n = top_n, other_label = other_label,
                                subset_label = lead, ctx = ctx)
    )
  }
  log_section("Building the pages")
  pages <- c(
    list(bar_plot(all_counted, subset_label = "All probes", ctx = ctx),
         bar_virus_plot(all_counted, subset_label = "All probes", ctx = ctx,
                        colours = virus_cols)),
    set_pages(all.main, "Main probes"),
    set_pages(all.accessory, "Accessory probes")
  )
  pages <- lapply(pages, stamp_warning_caption, caption = warn_caption)

  genera <- taxon_levels(all.full$label)
  key <- key_page(
    "Homology",
    paste("What the tBLASTn search found, after ranges were merged and globally",
          "reduced: how strong the hits are, which probes and viruses found them,",
          "and in which hosts. Every page shows the valid tier after global",
          "reduction. The main probe set comes first, then the accessory set."),
    colours = stats::setNames(unname(taxon_colours(genera)[genera]), genera),
    pages = page_titles(pages)
  )
  n_rows <- max(length(ctx$species_order), length(unique(all.full$species)))
  save_stage_pdf(
    c(list(key), pages), args$out_pdf,
    height = page_height_for(n_rows, per_species = cfg$plots$per_stratum %||% 0.18)
  )
  log_ok("wrote %s, %s pages", basename(args$out_pdf),
         format(length(pages) + 1L, big.mark = ","))
}


# ----------------------------------------------------------------------------
# Entry-point guard
# ----------------------------------------------------------------------------
# Only invoke main() under `Rscript plot2sort.R ...`. testthat sources this file
# inside test functions where sys.nframe() > 0, so unit tests get the helpers
# and builders without firing the CLI. run_main() logs the ending (ADR-021).
if (sys.nframe() == 0L) run_main(main)
