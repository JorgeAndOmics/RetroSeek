# =============================================================================
# taxonomy_segments.R
# =============================================================================
# Splits the authoritative catalog into per-SEGMENT deliverables (ADR-011).
#
# A "segment" is the taxon each locus rolls up to at `classification.segment_rank`
# — genus by default, but any NCBI rank works: the classifier derives the column
# by walking the reference taxonomy hierarchy, so nothing here (or there) hard-
# codes a taxon name. Loci whose call is coarser than the segment rank carry
# `unassigned_at_<rank>` and get their own segment, rather than being dropped or
# given invented precision.
#
# Outputs, under <out_dir>/by_<rank>/:
#   <segment>.csv                     one table per segment
#   segment_summary.csv               loci / species / HC counts per segment
#   plots/<segment>/*.png             a small curated panel per segment
#   segment_overview.png              all segments side by side
#
# Only a CURATED subset of the taxonomy panel is rendered per segment: the full
# 20-plot panel times N segments would be hundreds of PNGs for little gain.
# The builders are REUSED from taxonomy_plot_generator.R (sourced for its
# functions — its `if (sys.nframe() == 0L) main()` guard keeps the CLI dormant).

suppressMessages({
  library(argparse)
  library(tidyverse)
  library(yaml)
  library(ggsci)
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
# Overview: loci per segment, split by tier. The one plot that shows every
# segment at once, so it stays at the top level rather than inside a segment dir.
# ----------------------------------------------------------------------------
segment_overview_plot <- function(catalog) {
  if (nrow(catalog) == 0L || !"segment" %in% names(catalog)) {
    return(empty_plot("no segmented loci"))
  }
  counts <- catalog %>%
    count(segment = as.character(.data$segment),
          source = as.character(.data$source), name = "n")
  lvl <- counts %>% group_by(.data$segment) %>%
    summarise(t = sum(.data$n), .groups = "drop") %>%
    arrange(desc(.data$t)) %>% pull(.data$segment)
  counts$segment <- factor(counts$segment, levels = lvl)
  p <- ggplot(counts, aes(x = .data$segment, y = .data$n, fill = .data$source)) +
    geom_col() +
    scale_fill_manual(values = c(`ltr-flanked` = "#1F78B4", orphan = "#33A02C")) +
    labs(x = NULL, y = "loci", fill = "tier") +
    theme_bw()
  add_titles(p, "ERV loci per segment",
             "Per-segment burden, split by LTR-flanked vs orphan tier")
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
                      help = "directory for the by_<rank>/ tables + plots")
  parser$add_argument("--config", required = TRUE, help = "YAML config")
  parser$add_argument("--summary_csv", required = TRUE,
                      help = "per-segment summary CSV")
  args <- parser$parse_args()

  cfg <- yaml::read_yaml(args$config)
  `%||%` <- function(x, y) if (is.null(x)) y else x
  plot_dpi    <- cfg$plots$dpi    %||% 300
  plot_height <- cfg$plots$height %||% 12
  plot_width  <- cfg$plots$width  %||% 15
  per_stratum <- cfg$plots$per_stratum %||% 0.18
  max_dim     <- cfg$plots$max_dim     %||% 60
  seg_rank    <- cfg$classification$segment_rank %||% "genus"

  catalog <- readr::read_csv(args$catalog, show_col_types = FALSE)
  if (!"segment" %in% names(catalog)) {
    # Catalog predates ADR-011 (or segment_rank is unset): emit empty
    # deliverables so the DAG completes and the cause is visible in the log.
    log_section("catalog has no `segment` column - nothing to split")
    catalog$segment <- character(nrow(catalog))
  }
  log_section(sprintf("Segmenting %d catalog loci by %s", nrow(catalog), seg_rank))

  root <- file.path(args$output, paste0("by_", seg_rank))
  dir.create(root, showWarnings = FALSE, recursive = TRUE)

  summary_tbl <- segment_summary(catalog)
  dir.create(dirname(args$summary_csv), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(summary_tbl, args$summary_csv)

  n_species <- length(unique(catalog$species))
  overview <- scale_categorical_axis(
    segment_overview_plot(catalog), max(nrow(summary_tbl), 1L), axis = "x",
    base_w = plot_width, base_h = plot_height,
    per_stratum = per_stratum, cap = max_dim
  )
  save_plot("segment_overview.png", overview, root,
            dims = attr(overview, "intended_dims"),
            base_w = plot_width, base_h = plot_height, dpi = plot_dpi)

  # Per-segment tables + the curated plot subset.
  for (seg in summary_tbl$segment) {
    sub <- catalog %>% filter(as.character(.data$segment) == seg)
    stem <- safe_name(seg)
    readr::write_csv(sub, file.path(root, paste0(stem, ".csv")))

    pdir <- file.path(root, "plots", stem)
    dir.create(pdir, showWarnings = FALSE, recursive = TRUE)
    emit_seg <- function(name, plot) {
      plot <- scale_categorical_axis(plot, n_species, axis = "x",
                                     base_w = plot_width, base_h = plot_height,
                                     per_stratum = per_stratum, cap = max_dim)
      save_plot(name, plot, pdir, dims = attr(plot, "intended_dims"),
                base_w = plot_width, base_h = plot_height, dpi = plot_dpi)
    }
    # Curated subset: who carries it (composition), how sure we are
    # (confidence), and what shape they are in (structure).
    emit_seg("taxon_composition.png",        taxon_composition_plot(sub))
    emit_seg("confidence_gradient.png",      confidence_gradient_plot(sub))
    emit_seg("structure_class_composition.png", structure_class_composition_plot(sub))
  }
  log_section(sprintf("Done - wrote %d segment tables + plots to %s",
                      nrow(summary_tbl), root))
}


if (sys.nframe() == 0L) {
  .script_dir <- .resolve_script_dir()
  source(file.path(.script_dir, "..", "plot2sort", "helpers.R"))
  source(file.path(.script_dir, "..", "plot2sort", "io.R"))
  .t0 <- Sys.time()
  log_section <- function(name) {
    elapsed <- as.numeric(difftime(Sys.time(), .t0, units = "secs"))
    message(sprintf("[%6.2fs] > %s", elapsed, name))
  }
  # Reuse the taxonomy panel's builders rather than duplicating them. This also
  # pulls in its `main`, hence the distinct name above.
  source(file.path(.script_dir, "taxonomy_plot_generator.R"))
  segments_main()
}
