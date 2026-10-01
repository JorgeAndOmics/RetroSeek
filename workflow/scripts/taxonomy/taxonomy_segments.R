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
# Outputs are split by TYPE; only the tables keep a rank-agnostic by_<rank>/ level:
#   <out_dir>/by_<rank>/<segment>.csv    one table per segment
#   segment_summary.csv                  loci / species / HC counts per segment
#   <plots>/<segment>.pdf                the segment's pages, after a key page
#   <overview_pdf>                       every segment side by side
#
# How much is rendered per segment is set by `plots.segment_panel`:
#   full     (default) every page that means something within one segment: 21 of
#            the 24 taxonomy pages plus all 7 structure pages, so 28.
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
# Overview: loci per segment, split by tier.
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
    scale_y_count() +
    scale_fill_manual(values = .TIER_COLOUR, labels = display_label) +
    labs(x = NULL, y = "Loci", fill = NULL) +
    theme(panel.grid.major.y = element_blank())
  add_titles(p, "ERV loci per lineage",
             "Every locus by the lineage it rolls up to at the segment rank, by tier.")
}


# The study-wide overview: a key page, then loci per lineage. Written apart from
# the per-segment PDFs because it describes every segment at once.
#   catalog   the loaded catalog
#   seg_rank  the configured segment rank, named in the key page title
#   tiers     tier label -> colour, for the key page
#   path      where the PDF goes; height  page height in inches
write_overview_pdf <- function(catalog, seg_rank, tiers, path, height) {
  overview <- segment_overview_plot(catalog)
  key <- key_page(sprintf("ERV loci by %s", seg_rank),
                  paste("How the study's ERV loci divide among viral lineages at the",
                        "segment rank. Each lineage has its own PDF in the segments",
                        "folder beside this one."),
                  colours = tiers, pages = page_titles(list(overview)))
  save_stage_pdf(list(key, overview), path, height = height)
}


# Below this many loci a segment gets one page listing them: the full panel
# would draw some 28 pages that each restate the same few loci, most of them
# empty. Measured 2026-10-01: Bovispumavirus, one locus, 16 of 28 pages empty.
SPARSE_SEGMENT <- 5L

# The pages of one segment's PDF, and a sentence for its key page saying what
# was left out (NULL when nothing was).
#   sub       the segment's catalog rows (both tiers)
#   sub_loci  its LTR-flanked rows
#   panel     the registry entries to draw (segment_panel())
segment_pages <- function(sub, sub_loci, panel, ctx) {
  if (nrow(sub) < SPARSE_SEGMENT) {
    return(list(pages = list(segment_loci_page(sub)),
                note = sprintf("With only %s, one page lists them instead.",
                               .count_of(nrow(sub), "locus", "loci"))))
  }
  note <- if (nrow(sub_loci) == 0L) {
    "This lineage has no LTR-flanked loci, so the pages about them are left out."
  }
  list(pages = render_panel(panel, sub_loci, sub, ctx), note = note)
}

# " Nearest virus: X at Y% identity." per locus; "no reference virus found" when
# the classifier found none (a blank cell), nothing at all for a catalog written
# before ADR-024, which has no such columns.
.nearest_phrase <- function(sub) {
  if (!all(c("nearest_virus", "nearest_virus_identity") %in% names(sub))) {
    return(rep("", nrow(sub)))
  }
  found <- !is.na(sub$nearest_virus) & nzchar(sub$nearest_virus)
  ifelse(found,
         sprintf(" Nearest virus: %s at %s%% identity.", sub$nearest_virus,
                 sub$nearest_virus_identity),
         " Nearest virus: no reference virus found.")
}

# One block per locus: the host as a heading, then where it is, what it
# carries and how it was called. For the few loci of a sparse segment.
segment_loci_page <- function(sub) {
  confidence <- as.numeric(sub$confidence)
  confidence <- ifelse(is.na(confidence), "unknown", sprintf("%.2f", confidence))
  template <- "%s, %s:%s to %s (%s). Genes: %s. Called %s at confidence %s.%s"
  details <- sprintf(template, display_label(sub$source), sub$seqname,
                     scales::comma(as.numeric(sub$start)),
                     scales::comma(as.numeric(sub$end)), sub$strand,
                     gsub(",", ", ", sub$genes_present), display_label(sub$taxon_call),
                     confidence, .nearest_phrase(sub))
  details <- vapply(details, function(d) paste(strwrap(d, 160), collapse = "\n"),
                    character(1), USE.NAMES = FALSE)
  y <- 1 - (seq_len(nrow(sub)) - 1) * 0.18
  p <- ggplot() +
    annotate("text", x = 0, y = y, label = sub$species, hjust = 0, vjust = 1,
             size = 4.2, fontface = "bold.italic", family = .FONT, colour = .INK) +
    annotate("text", x = 0, y = y - 0.04, label = details, hjust = 0, vjust = 1,
             size = 3.6, family = .FONT, colour = .INK, lineheight = 1.2) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_retroseek_blank()
  add_titles(p, "The loci of this lineage",
             "Too few loci for the usual pages, so each is listed with its evidence.")
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
  parser$add_argument("--plots", required = TRUE,
                      help = "directory for one PDF per segment")
  parser$add_argument("--overview_pdf", required = TRUE,
                      help = "path of the overview PDF of every segment")
  parser$add_argument("--species_tree_dir", required = FALSE, default = "",
                      help = "species_tree_layout.py output: the host tree coordinates")
  parser$add_argument("--tree_dir", required = FALSE, default = "",
                      help = paste("directory of tree-coordinate CSVs, so the",
                                   "species-tree panels can render per segment"))
  parser$add_argument("--config", required = TRUE, help = "YAML config")
  parser$add_argument("--summary_csv", required = TRUE,
                      help = "per-segment summary CSV")
  parser$add_argument("--log", default = NULL,
                      help = "job log file; the Snakemake log: path")
  args <- parser$parse_args()
  log_job(args$log, "taxonomy_segments")

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
  # results/plots/ is figures. The tables keep a by_<rank>/ level; the figures
  # do not, since results/plots is one level deep and each key page names the
  # rank (docs/visual_style.md, "Output").
  root <- file.path(args$output, paste0("by_", seg_rank))
  dir.create(root, showWarnings = FALSE, recursive = TRUE)
  plot_root <- args$plots
  dir.create(plot_root, showWarnings = FALSE, recursive = TRUE)

  summary_tbl <- segment_summary(catalog)
  dir.create(dirname(args$summary_csv), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(summary_tbl, args$summary_csv)

  tiers <- stats::setNames(unname(.TIER_COLOUR[c("ltr-flanked", "orphan")]),
                           display_label(c("ltr-flanked", "orphan")))
  per_species <- cfg$plots$per_stratum %||% 0.18
  overview_height <- page_height_for(nrow(summary_tbl), per_species = per_species)
  write_overview_pdf(catalog, seg_rank, tiers,
                     path = args$overview_pdf,
                     height = overview_height)

  # Per-segment tables and pages. The registry is resolved once: it is the same
  # declaration the full panels use, filtered to entries that mean something for
  # a single segment (erv_class is constant within a genus, and the taxonomy
  # cladogram collapses to one tip).
  panel <- segment_panel(c(panel_registry(), structure_panel_registry()), panel_mode)
  ctx <- panel_ctx(cfg, args$species_tree_dir %||% "", tree_dir = args$tree_dir %||% "")
  height <- page_height_for(max(length(ctx$species_order),
                                length(unique(catalog$species))),
                            per_species = per_species)
  log_section(sprintf("Panel mode '%s': %d pages per segment",
                      panel_mode, length(panel)))
  for (seg in summary_tbl$segment) {
    sub <- catalog %>% filter(as.character(.data$segment) == seg)
    stem <- safe_name(seg)
    readr::write_csv(sub, file.path(root, paste0(stem, ".csv")))

    if (length(panel) == 0L) next
    # `data` keeps the tier scope: composition and mosaic pages are LTR-flanked
    # only, so orphans are not silently mixed in.
    sub_loci <- sub %>% filter(as.character(.data$source) == "ltr-flanked")
    out <- segment_pages(sub, sub_loci, panel, ctx)
    key <- key_page(
      key_title("ERV loci:", display_label(seg)),
      paste(sprintf(paste("The taxonomy and structure pages, restricted to the %s",
                          "of this lineage. Hosts are rows in the order of the host",
                          "tree."),
                    .count_of(nrow(sub), "locus", "loci")), out$note %||% ""),
      colours = tiers, pages = page_titles(out$pages)
    )
    save_stage_pdf(c(list(key), out$pages), file.path(plot_root, paste0(stem, ".pdf")),
                   height = height)
  }
  log_info("wrote %d segment tables to %s and PDFs to %s",
           nrow(summary_tbl), root, plot_root)
  log_ok("%s loci in %s segments by %s", format(nrow(catalog), big.mark = ","),
         format(nrow(summary_tbl), big.mark = ","), seg_rank)
}


if (sys.nframe() == 0L) {
  .script_dir <- .resolve_script_dir()
  # Line contract, run_main (ADR-021).
  source(file.path(.script_dir, "..", "utils", "log.R"))
  # Palette, theme, labels, stage PDFs.
  source(file.path(.script_dir, "..", "plot2sort", "style.R"))

  source(file.path(.script_dir, "..", "plot2sort", "helpers.R"))
  source(file.path(.script_dir, "..", "plot2sort", "tree_axis.R"))
  # Reuse the taxonomy panel's builders rather than duplicating them. This also
  # pulls in its `main`, hence the distinct name above.
  source(file.path(.script_dir, "taxonomy_plot_generator.R"))
  # Also needed for structure_panel_registry(); CLI-guarded the same way.
  source(file.path(.script_dir, "erv_like_plot_generator.R"))
  # Everything is sourced before log_job() runs inside segments_main(): sourcing
  # log.R again (the two generators do) would reset its state.
  run_main(segments_main)
}
