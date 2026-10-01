# =============================================================================
# highlights.R - the headline findings of a study, in a few slide-shaped pages
# =============================================================================
# Purpose: one short PDF to open first or to show: how many ERV loci each host
# carries, which lineages they belong to, every lineage against every host, and
# the study's focal lineages. Every other figure lives in its stage PDF.
#
# Layout: 16:9 pages with hosts as COLUMNS under the host tree (species_columns()
# and on_columns() in plot2sort/tree_axis.R), the one exception to species on
# rows (docs/visual_style.md): a wide page holds more hosts across than down, and
# the pages drop into a slide deck unchanged.
#
# Inputs:
#   --catalog           results/tables/taxonomy_classification/catalog.csv; only
#                       species, source and segment are read
#   --config            the run config: plots.focal_lineages (may be empty),
#                       classification.segment_rank, species (the column order
#                       without a tree)
#   --species_tree_dir  species_tree_layout.py output, or "" for no tree
# Output:
#   --out_pdf           results/plots/highlights.pdf: a key page, then loci per
#                       host, the lineage mix, the lineage matrix and, when
#                       focal lineages are listed, the focal page. A study with
#                       no loci gets the key page alone.
# Fails when a catalog with rows lacks species, source or segment.
# =============================================================================

suppressMessages({
  library(argparse)
  library(dplyr)
  library(ggplot2)
  library(yaml)
})

.resolve_script_dir <- function() {
  for (frame in rev(sys.frames())) {
    if (!is.null(frame$ofile)) {
      return(dirname(normalizePath(frame$ofile, mustWork = FALSE)))
    }
  }
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) return(dirname(sub("^--file=", "", file_arg[1])))
  "scripts"
}
.script_dir <- .resolve_script_dir()
source(file.path(.script_dir, "utils", "log.R"))
source(file.path(.script_dir, "plot2sort", "style.R"))
source(file.path(.script_dir, "plot2sort", "helpers.R"))
source(file.path(.script_dir, "plot2sort", "tree_axis.R"))

# The two catalog tiers in their house colours, keyed by tier.
.TIERS <- .TIER_COLOUR[c("ltr-flanked", "orphan")]


# Every page reads one small table: loci counted once per host, tier and lineage
# (a few thousand rows at 102 genomes, not the whole catalog).
count_cells <- function(catalog) {
  catalog %>% count(.data$species, .data$source, .data$segment, name = "n")
}

# Loci per host and lineage, both tiers together.
.lineage_counts <- function(cells) {
  cells %>% count(.data$species, .data$segment, wt = .data$n, name = "n")
}


# Page: how many loci each host carries, LTR-flanked and orphan.
loci_by_tier_page <- function(cells, axis) {
  tiers <- cells %>% count(.data$species, .data$source, wt = .data$n, name = "n")
  panel <- column_count_panel(tiers, "source", .TIERS, display_label(names(.TIERS)),
                              length(axis$levels))
  on_columns(list(panel), axis, title = "ERV loci per host",
             subtitle = paste("Every catalogued locus: LTR-flanked proviruses and the",
                              "orphans recovered by their coding sequence."))
}

# Page: each host's lineages, as counts and as shares.
lineage_mix_page <- function(cells, axis) {
  counts <- .lineage_counts(cells)
  order <- taxon_levels(counts$segment, counts$n)
  colours <- taxon_colours(order)
  count_panel <- column_count_panel(counts, "segment", colours, taxon_labels(order),
                                    length(axis$levels))
  # The same fill scale as the count panel, so the two share one legend.
  share_panel <- ggplot(counts, aes(x = .data$species, y = .data$n,
                                    fill = .data$segment)) +
    geom_col(width = 0.72, position = position_fill(reverse = TRUE)) +
    scale_fill_manual(values = colours, breaks = order, labels = taxon_labels(order),
                      limits = order, name = NULL) +
    scale_y_continuous(labels = scales::label_percent(), expand = c(0, 0)) +
    labs(y = "Share of loci") +
    theme(panel.grid.major.x = element_blank())
  on_columns(list(count_panel, share_panel), axis, title = "Viral lineages per host",
             subtitle = paste("Loci by the lineage they roll up to at the segment",
                              "rank, as counts and as each host's shares."))
}


# One row per host and lineage with its locus count, lineages ordered by size.
lineage_matrix_cells <- function(cells) {
  counts <- .lineage_counts(cells)
  counts$segment <- factor(counts$segment,
                           levels = rev(taxon_levels(counts$segment, counts$n)))
  counts
}

# Page: every lineage against every host, one bubble per pair sized by its loci.
lineage_matrix_page <- function(cells, axis) {
  counts <- lineage_matrix_cells(cells)
  panel <- ggplot(counts, aes(x = .data$species, y = .data$segment)) +
    geom_point(aes(size = .data$n, colour = .data$segment), alpha = 0.85) +
    # The count beside each bubble: a rare lineage is a small dot, not nothing.
    geom_text(aes(label = scales::comma(.data$n)), nudge_x = 0.12, hjust = 0,
              size = 2.6, colour = .INK_SOFT, family = .FONT) +
    # A floor on the dot size, so a lineage with one locus still shows.
    scale_size(range = c(1.2, 12), labels = scales::label_comma(), name = "Loci",
               breaks = scales::breaks_extended(n = 4)) +
    scale_colour_manual(values = taxon_colours(levels(counts$segment)),
                        guide = "none") +
    scale_y_discrete(labels = taxon_labels) +
    labs(y = NULL)
  on_columns(list(panel), axis, title = "Which lineages each host carries",
             subtitle = paste("One bubble per host and lineage, sized and labelled by",
                              "its loci. An absent lineage is a gap."))
}


# A focal panel's title: a taxon in bold italics, a non-taxon in bold.
.focal_title <- function(lineage) {
  name <- display_label(lineage)
  if (grepl(.NOT_A_TAXON, name)) name else bquote(bolditalic(.(name)))
}

# Page: the study's focal lineages, one panel each under one tree, or NULL when
# none are configured. A lineage with no loci keeps its panel, saying so, and
# the log names it, since a misspelt name in plots.focal_lineages looks the same.
focal_lineages_page <- function(cells, axis, focal) {
  if (!length(focal)) return(NULL)
  counts <- .lineage_counts(cells)
  n <- length(axis$levels)
  panels <- lapply(focal, function(lineage) {
    d <- counts %>% filter(.data$segment == lineage)
    if (!nrow(d)) {
      log_warn("focal lineage %s has no loci in this catalog: check the name in %s",
               lineage, "plots.focal_lineages")
      return(no_loci_column_panel(n, title = .focal_title(lineage)))
    }
    column_count_panel(d, "segment", stats::setNames(taxon_colours(lineage), lineage),
                       taxon_labels(lineage), n, title = .focal_title(lineage))
  })
  page <- on_columns(panels, axis, title = "The focal lineages",
                     subtitle = paste("Loci per host for each lineage this study looks",
                                      "at (plots.focal_lineages); each panel has its",
                                      "own scale."))
  # Each panel is titled with its lineage, so a legend would only repeat them.
  page & theme(legend.position = "none", plot.title = element_text(size = 11))
}


# The pages after the key page, in reading order; none for a study without loci,
# the focal page only when the config lists focal lineages.
highlights_pages <- function(catalog, axis, focal) {
  if (!nrow(catalog)) return(list())
  cells <- count_cells(catalog)
  pages <- list(loci_by_tier_page(cells, axis), lineage_mix_page(cells, axis),
                lineage_matrix_page(cells, axis),
                focal_lineages_page(cells, axis, focal))
  Filter(Negate(is.null), pages)
}

# The three catalog columns this PDF reads, checked once at the boundary. A blank
# segment becomes "unassigned_at_<rank>", as the classifier writes it.
.read_catalog <- function(path, segment_rank = "genus") {
  wanted <- c("species", "source", "segment")
  catalog <- readr::read_csv(path, col_select = dplyr::any_of(wanted),
                             col_types = readr::cols(.default = "c"), na = "")
  # A study without loci writes a catalog with no columns at all: valid, empty.
  if (!nrow(catalog)) {
    return(tibble::tibble(species = character(), source = character(),
                          segment = character()))
  }
  missing <- setdiff(wanted, names(catalog))
  if (length(missing)) {
    abort_hint(sprintf("%s lacks the column(s) %s", path,
                       paste(missing, collapse = ", ")),
               "rerun --classify to write a current catalog")
  }
  catalog %>% mutate(segment = coalesce(.data$segment,
                                        paste0("unassigned_at_", segment_rank)))
}


main <- function() {
  parser <- ArgumentParser(description = "The headline findings, slide-shaped.")
  parser$add_argument("--catalog", required = TRUE, help = "catalog.csv")
  parser$add_argument("--config", required = TRUE, help = "YAML config file")
  parser$add_argument("--species_tree_dir", default = "",
                      help = "species_tree_layout.py output: the host tree coordinates")
  parser$add_argument("--out_pdf", required = TRUE,
                      help = "results/plots/highlights.pdf")
  parser$add_argument("--log", default = NULL,
                      help = "job log file; the Snakemake log: path")
  args <- parser$parse_args()
  log_job(args$log, "highlights")
  use_retroseek_style()

  cfg <- yaml::read_yaml(args$config)
  rank <- cfg$classification$segment_rank
  catalog <- .read_catalog(args$catalog, if (is.null(rank)) "genus" else rank)
  ctx <- panel_ctx(cfg, args$species_tree_dir)
  axis <- species_columns(catalog$species, ctx$species_tree, ctx$species_order)
  pages <- highlights_pages(catalog, axis, unlist(cfg$plots$focal_lineages))
  key <- key_page("Highlights",
                  paste("The headline findings of this study, one page each. Hosts are",
                        "columns under the host tree. Every other figure is in its",
                        "stage PDF.",
                        if (!nrow(catalog)) "This study has no loci: no pages follow."),
                  colours = tier_key_colours(), pages = page_titles(pages))
  save_stage_pdf(c(list(key), pages), args$out_pdf, width = .SLIDE_WIDTH,
                 height = .SLIDE_HEIGHT)
  log_ok("wrote %s, %d pages, %s loci", basename(args$out_pdf), length(pages) + 1L,
         format(nrow(catalog), big.mark = ","))
}


if (sys.nframe() == 0L) run_main(main)
