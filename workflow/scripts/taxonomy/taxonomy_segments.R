# =============================================================================
# taxonomy_segments.R
# =============================================================================
# Splits the authoritative catalog into per-SEGMENT deliverables (ADR-011).
#
# A "segment" is the taxon each locus rolls up to at `classification.segment_rank`
# - genus by default, but any NCBI rank works: the classifier derives the column
# by walking the reference taxonomy hierarchy, so nothing here (or there) hard-
# codes a taxon name. Loci whose call is coarser than the segment rank carry
# `unassigned_at_<rank>` and get their own segment, rather than being dropped or
# given invented precision.
#
# Outputs are split by TYPE, both under a rank-agnostic by_<rank>/ level:
#   <out_dir>/by_<rank>/<segment>.csv    one table per segment
#   segment_summary.csv                  loci / species / HC counts per segment
#   <plots>/by_<rank>/<segment>.pdf      the segment's pages, after a key page
#   <plots>/by_<rank>/overview.pdf       every segment side by side
#
# How much is rendered per segment is set by `plots.segment_panel`:
#   full     (default) every page that means something within one segment: 20 of
#            the 23 taxonomy pages plus all 7 structure pages, so 27.
#   curated  the small legacy subset of 3.
#   none     tables only.
# Cost scales as segments x pages and each page has a row per species, so drop
# to `curated` if a high-genome-count run gets bulky.
# The builders are REUSED from taxonomy_plot_generator.R and
# erv_like_plot_generator.R (sourced for their functions; their
# `if (sys.nframe() == 0L) main()` guards keep the CLIs dormant).

suppressMessages({
  library(argparse)
  library(tidyverse)
  library(yaml)
})


# ----------------------------------------------------------------------------
# Locate sibling scripts + source shared modules
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


# ----------------------------------------------------------------------------
# A filesystem-safe segment name. Segment values come from NCBI taxonomy, so
# they are normally clean, but a name with a space or slash must not create a
# stray directory level or collide with a sibling.
# ----------------------------------------------------------------------------
safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  x[!nzchar(x)] <- "unnamed"
  x
}


# ----------------------------------------------------------------------------
# Per-segment summary: how many loci, across how many species, how many are
# high-confidence, and the tier split. One row per segment, largest first.
# Empty-safe.
# ----------------------------------------------------------------------------
segment_summary <- function(catalog) {
  empty <- tibble(
    segment = character(), n_loci = integer(), n_species = integer(),
    n_hc = integer(), n_ltr_flanked = integer(), n_orphan = integer()
  )
  if (nrow(catalog) == 0L || !"segment" %in% names(catalog)) return(empty)
  catalog %>%
    group_by(segment = as.character(.data$segment)) %>%
    summarise(
      n_loci        = dplyr::n(),
      n_species     = dplyr::n_distinct(.data$species),
      n_hc          = sum(.data$confidence_tag == "HC", na.rm = TRUE),
      n_ltr_flanked = sum(.data$source == "ltr-flanked", na.rm = TRUE),
      n_orphan      = sum(.data$source == "orphan", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(.data$n_loci), .data$segment)
}


# ----------------------------------------------------------------------------
# Load the catalog the way the reused builders expect it.
#
# Two properties, both easy to lose and both failing SILENTLY:
#
#   Everything character. The catalog writes booleans as the strings
#   "True"/"False". readr's default inference turns `resolved` into a logical,
#   and the builders filter on `resolved == "True"` - TRUE coerces to "TRUE",
#   which is not "True", so the filter matches nothing and every segment gets
#   the same "no confident taxon calls" placeholder.
#
#   Numeric companions. `confidence_gradient_plot` guards on `confidence_num`
#   and returns its placeholder when absent. That column is derived, not stored;
#   taxonomy_plot_generator.R's main() makes it, and this path never runs main().
#
# Neither raises, so the symptom was N byte-identical PNGs rather than an error.
# ----------------------------------------------------------------------------
load_catalog <- function(path) {
  catalog <- readr::read_csv(
    path,
    col_types = readr::cols(.default = readr::col_character())
  )
  # BOTH companion sets: the taxonomy builders need confidence_num/n_hits, the
  # structure builders need span_bp and typed completeness/n_main_genes. Each
  # generator derives its own inside a loader this path does not use.
  add_structure_companions(add_numeric_companions(catalog))
}


# ----------------------------------------------------------------------------
# Overview: loci per segment, split by tier. The one page that shows every
# segment at once, so it has its own PDF beside the per-segment ones.
# ----------------------------------------------------------------------------
segment_overview_plot <- function(catalog) {
  if (nrow(catalog) == 0L || !"segment" %in% names(catalog)) {
    return(empty_plot("No segmented loci"))
  }
  counts <- catalog %>%
    count(segment = as.character(.data$segment),
          source = as.character(.data$source), name = "n")
  # Largest segment on top; "unassigned" last whatever its size.
  counts$segment <- factor(counts$segment,
                           levels = rev(taxon_levels(counts$segment, counts$n)))
  p <- ggplot(counts, aes(x = .data$segment, y = .data$n, fill = .data$source)) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
    coord_flip() +
    scale_x_discrete(labels = taxon_labels) +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_manual(values = .TIER_COLOUR, labels = display_label) +
    labs(x = NULL, y = "Loci", fill = NULL) +
    theme(panel.grid.major.y = element_blank())
  add_titles(p, "ERV loci per lineage",
             "Every locus by the lineage it rolls up to at the segment rank, by tier.")
}


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
# NOT called `main`: this script sources taxonomy_plot_generator.R to reuse its
# builders, and that file defines its own `main()`. A shared name would let the
# sourced definition win and this stage would run the plot generator's CLI.
segments_main <- function() {
  parser <- ArgumentParser(description = "Split the ERV catalog by taxonomic segment.")
  parser$add_argument("--catalog", required = TRUE,
                      help = "authoritative catalog.csv from taxonomy_plot_generator")
  parser$add_argument("--output", required = TRUE,
                      help = "directory for the by_<rank>/ tables")
  parser$add_argument("--plots", required = FALSE, default = NULL,
                      help = paste("directory for the by_<rank>/ figures.",
                                   "Defaults to --output for standalone use;",
                                   "the pipeline points it at results/plots/."))
  parser$add_argument("--species_tree_dir", required = FALSE, default = "",
                      help = "species_tree_layout.py output: the host tree coordinates")
  parser$add_argument("--tree_dir", required = FALSE, default = "",
                      help = paste("directory of tree-coordinate CSVs, so the",
                                   "species-tree panels can render per segment"))
  parser$add_argument("--config", required = TRUE, help = "YAML config")
  parser$add_argument("--summary_csv", required = TRUE,
                      help = "per-segment summary CSV")
  args <- parser$parse_args()

  cfg <- yaml::read_yaml(args$config)
  `%||%` <- function(x, y) if (is.null(x)) y else x
  seg_rank    <- cfg$classification$segment_rank %||% "genus"
  # full = every panel entry meaningful for one segment; curated = the small
  # legacy subset; none = tables only. See docs/configuration.md.
  panel_mode  <- cfg$plots$segment_panel %||% "full"
  use_retroseek_style()

  catalog <- load_catalog(args$catalog)
  if (!"segment" %in% names(catalog)) {
    # Catalog predates ADR-011 (or segment_rank is unset): emit empty
    # deliverables so the DAG completes and the cause is visible in the log.
    log_section("catalog has no `segment` column, so nothing to split")
    catalog$segment <- character(nrow(catalog))
  }
  log_section(sprintf("Segmenting %d catalog loci by %s", nrow(catalog), seg_rank))

  # Tables and figures split by TYPE, not by stage: results/tables/ is CSV and
  # results/plots/ is figures. Both keep the rank-agnostic by_<rank>/ level.
  root <- file.path(args$output, paste0("by_", seg_rank))
  dir.create(root, showWarnings = FALSE, recursive = TRUE)
  plot_root <- file.path(args$plots %||% args$output, paste0("by_", seg_rank))
  dir.create(plot_root, showWarnings = FALSE, recursive = TRUE)

  summary_tbl <- segment_summary(catalog)
  dir.create(dirname(args$summary_csv), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(summary_tbl, args$summary_csv)

  tiers <- stats::setNames(unname(.TIER_COLOUR[c("ltr-flanked", "orphan")]),
                           display_label(c("ltr-flanked", "orphan")))
  overview <- segment_overview_plot(catalog)
  save_stage_pdf(
    list(key_page(sprintf("ERV loci by %s", seg_rank),
                  paste("How the study's ERV loci divide among viral lineages at the",
                        "segment rank. Each lineage has its own PDF beside this one."),
                  colours = tiers, pages = page_titles(list(overview))),
         overview),
    file.path(plot_root, "overview.pdf"),
    height = page_height_for(nrow(summary_tbl),
                             per_species = cfg$plots$per_stratum %||% 0.18))

  # Per-segment tables and pages. The registry is resolved once: it is the same
  # declaration the full panels use, filtered to entries that mean something for
  # a single segment (erv_class is constant within a genus, and the taxonomy
  # cladogram collapses to one tip).
  panel <- segment_panel(c(panel_registry(), structure_panel_registry()), panel_mode)
  ctx <- panel_ctx(cfg, args$species_tree_dir %||% "", tree_dir = args$tree_dir %||% "")
  height <- page_height_for(max(length(ctx$species_order), length(unique(catalog$species))),
                            per_species = cfg$plots$per_stratum %||% 0.18)
  log_section(sprintf("Panel mode '%s': %d pages per segment", panel_mode, length(panel)))
  for (seg in summary_tbl$segment) {
    sub <- catalog %>% filter(as.character(.data$segment) == seg)
    stem <- safe_name(seg)
    readr::write_csv(sub, file.path(root, paste0(stem, ".csv")))

    if (length(panel) == 0L) next
    # `data` keeps the tier scope: composition and mosaic pages are LTR-flanked
    # only, so orphans are not silently mixed in.
    sub_loci <- sub %>% filter(as.character(.data$source) == "ltr-flanked")
    pages <- render_panel(panel, sub_loci, sub, ctx)
    key <- key_page(
      sprintf("ERV loci: %s", display_label(seg)),
      sprintf(paste("The taxonomy and structure pages, restricted to the %s loci",
                    "of this lineage. Hosts are rows in the order of the host tree."),
              scales::comma(nrow(sub))),
      colours = tiers, pages = page_titles(pages))
    save_stage_pdf(c(list(key), pages), file.path(plot_root, paste0(stem, ".pdf")),
                   height = height)
  }
  log_section(sprintf("Done: wrote %d segment tables to %s and PDFs to %s",
                      nrow(summary_tbl), root, plot_root))
}


if (sys.nframe() == 0L) {
  .script_dir <- .resolve_script_dir()
  source(file.path(.script_dir, "..", "plot2sort", "style.R"))  # palette, theme, labels, stage PDFs
  source(file.path(.script_dir, "..", "plot2sort", "helpers.R"))
  source(file.path(.script_dir, "..", "plot2sort", "tree_axis.R"))
  .t0 <- Sys.time()
  log_section <- function(name) {
    elapsed <- as.numeric(difftime(Sys.time(), .t0, units = "secs"))
    message(sprintf("[%6.2fs] > %s", elapsed, name))
  }
  # Reuse the taxonomy panel's builders rather than duplicating them. This also
  # pulls in its `main`, hence the distinct name above.
  source(file.path(.script_dir, "taxonomy_plot_generator.R"))
  # Also needed for structure_panel_registry(); CLI-guarded the same way.
  source(file.path(.script_dir, "erv_like_plot_generator.R"))
  segments_main()
}
