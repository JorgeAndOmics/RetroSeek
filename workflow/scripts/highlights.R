# =============================================================================
# highlights.R - the headline findings of a study, in a few slide-shaped pages
# =============================================================================
# Purpose: one short PDF to open first or to show: how many ERV loci each host
# carries, which lineages they belong to, and the focal lineages of the study.
# Every other figure lives in its stage PDF; this one only summarises.
#
# Layout: 16:9 pages with hosts as COLUMNS under the host tree (species_columns()
# and on_columns() in plot2sort/tree_axis.R), the one exception to species on
# rows (docs/visual_style.md): a wide page holds more hosts across than down, and
# the pages drop into a slide deck unchanged.
#
# Inputs:
#   --catalog           results/tables/taxonomy_classification/catalog.csv
#   --config            the run config: plots.focal_lineages (may be empty),
#                       species (the column order without a tree)
#   --species_tree_dir  species_tree_layout.py output, or "" for no tree
# Output:
#   --out_pdf           results/plots/highlights.pdf: a key page, then three
#                       pages, or four when focal lineages are listed
# Fails when the catalog is missing or lacks species, source or segment.
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

# A 16:9 page in inches, the size of a widescreen slide.
.PAGE_16_9 <- c(width = 13.33, height = 7.5)


# Loci per host stacked by `fill`, hosts on a discrete x axis. A total sits above
# each bar while the hosts are few enough for the numbers not to collide.
.count_panel <- function(d, fill, colours, labels, n_species, title = NULL) {
  totals <- d %>% count(.data$species, name = "n")
  p <- ggplot(d, aes(x = .data$species, fill = .data[[fill]])) +
    # show.legend keeps a key for a class a panel lacks, so stacked panels share
    # one legend.
    geom_bar(width = 0.72, position = position_stack(reverse = TRUE),
             show.legend = TRUE) +
    scale_fill_manual(values = colours, breaks = names(colours), labels = labels,
                      limits = names(colours), drop = FALSE, name = NULL) +
    scale_y_count(expand = expansion(mult = c(0, 0.22))) +
    labs(y = "Loci", title = title) +
    theme(panel.grid.major.x = element_blank())
  if (n_species > 30L || !nrow(totals)) return(p)
  p + geom_text(data = totals, aes(x = .data$species, y = .data$n,
                                   label = scales::comma(.data$n)),
                inherit.aes = FALSE, vjust = -0.4, size = 3, colour = .INK_SOFT)
}

# The two tiers in their house colours.
.tier_colours <- function() .TIER_COLOUR[c("ltr-flanked", "orphan")]


# Page: how many loci each host carries, LTR-flanked and orphan.
loci_by_tier_page <- function(catalog, axis) {
  panel <- .count_panel(catalog, "source", .tier_colours(),
                        display_label(names(.tier_colours())), length(axis$levels))
  on_columns(list(panel), axis, .PAGE_16_9[["width"]], tree_share = 0.22,
             title = "ERV loci per host",
             subtitle = paste("Every catalogued locus: LTR-flanked proviruses and the",
                              "orphans recovered by their coding sequence."))
}

# Page: each host's lineages, as counts and as shares.
lineage_mix_page <- function(catalog, axis) {
  weights <- table(catalog$segment)
  order <- taxon_levels(names(weights), as.numeric(weights))
  colours <- taxon_colours(order)
  counts <- .count_panel(catalog, "segment", colours, taxon_labels(order),
                         length(axis$levels))
  shares <- ggplot(catalog, aes(x = .data$species, fill = .data$segment)) +
    geom_bar(width = 0.72, position = position_fill(reverse = TRUE)) +
    scale_fill_manual(values = colours, breaks = order, labels = taxon_labels(order),
                      limits = order, name = NULL) +
    scale_y_continuous(labels = scales::label_percent(), expand = c(0, 0)) +
    labs(y = "Share of loci") +
    theme(panel.grid.major.x = element_blank())
  on_columns(list(counts, shares), axis, .PAGE_16_9[["width"]], tree_share = 0.14,
             title = "Viral lineages per host",
             subtitle = paste("Loci by the lineage they roll up to at the segment",
                              "rank, as counts and as each host's shares."))
}


# One row per host and lineage with its locus count, lineages ordered by size.
genus_matrix_cells <- function(catalog) {
  cells <- catalog %>% count(.data$species, .data$segment, name = "n")
  totals <- tapply(cells$n, cells$segment, sum)
  cells$segment <- factor(cells$segment,
                          levels = rev(taxon_levels(names(totals), as.numeric(totals))))
  cells
}

# Page: every lineage against every host, one bubble per pair sized by its loci.
genus_matrix_page <- function(catalog, axis) {
  cells <- genus_matrix_cells(catalog)
  panel <- ggplot(cells, aes(x = .data$species, y = .data$segment)) +
    geom_point(aes(size = .data$n, colour = .data$segment), alpha = 0.85) +
    # The count beside each bubble: a rare lineage is a small dot, not nothing.
    geom_text(aes(label = scales::comma(.data$n)), nudge_x = 0.12, hjust = 0,
              size = 2.6, colour = .INK_SOFT, family = .FONT) +
    # A floor on the dot size, so a lineage with one locus still shows.
    scale_size(range = c(1.2, 12), labels = scales::label_comma(), name = "Loci",
               breaks = scales::breaks_extended(n = 4)) +
    scale_colour_manual(values = taxon_colours(levels(cells$segment)), guide = "none") +
    scale_y_discrete(labels = taxon_labels) +
    labs(y = NULL)
  on_columns(list(panel), axis, .PAGE_16_9[["width"]], tree_share = 0.14,
             title = "Which lineages each host carries",
             subtitle = paste("One bubble per host and lineage, sized and labelled by",
                              "its loci. An absent lineage is a gap."))
}


# The panel of a focal lineage with no loci: it keeps its place and says so, and
# the log names it, since a misspelt name in plots.focal_lineages looks the same.
.no_loci_panel <- function(lineage, n_columns) {
  log_warn("focal lineage %s has no loci in this catalog: check the name in %s",
           lineage, "plots.focal_lineages")
  ggplot() +
    annotate("text", x = (n_columns + 1) / 2, y = 0.5, label = "No loci",
             colour = .INK_SOFT, size = 3.5, family = .FONT) +
    labs(title = display_label(lineage), y = "Loci") +
    theme(plot.title = element_text(face = "bold.italic", size = 11),
          axis.text.y = element_blank(), panel.grid = element_blank())
}

# Page: the study's focal lineages, one panel each under one tree, or NULL when
# none are configured. A lineage with no loci keeps its panel, saying so, so the
# layout never shifts between studies.
focal_lineages_page <- function(catalog, axis, focal) {
  if (!length(focal)) return(NULL)
  panels <- lapply(focal, function(lineage) {
    d <- catalog %>% filter(.data$segment == lineage)
    if (!nrow(d)) return(.no_loci_panel(lineage, length(axis$levels)))
    colours <- stats::setNames(taxon_colours(lineage), lineage)
    .count_panel(d, "segment", colours, taxon_labels(lineage), length(axis$levels),
                 title = display_label(lineage)) +
      theme(plot.title = element_text(face = "bold.italic", size = 11))
  })
  # Each panel is titled with its lineage, so a legend would only repeat them.
  on_columns(panels, axis, .PAGE_16_9[["width"]], legend = "none",
             title = "The focal lineages",
             subtitle = paste("Loci per host for each lineage this study looks at",
                              "(plots.focal_lineages); each panel has its own scale."))
}


# The pages after the key page, in reading order; the focal page only when the
# config lists focal lineages.
highlights_pages <- function(catalog, axis, focal) {
  if (!nrow(catalog)) return(list())   # a study without loci: the key page only
  pages <- list(loci_by_tier_page(catalog, axis), lineage_mix_page(catalog, axis),
                genus_matrix_page(catalog, axis),
                focal_lineages_page(catalog, axis, focal))
  Filter(Negate(is.null), pages)
}

# The catalog columns this PDF reads, checked once at the boundary.
.read_catalog <- function(path) {
  catalog <- readr::read_csv(path, col_types = readr::cols(.default = "c"), na = "")
  # A study without loci writes a catalog with no columns at all: valid, empty.
  if (!nrow(catalog)) {
    return(tibble::tibble(species = character(), source = character(),
                          segment = character()))
  }
  missing <- setdiff(c("species", "source", "segment"), names(catalog))
  if (length(missing)) {
    abort_hint(sprintf("%s lacks the column(s) %s", path,
                       paste(missing, collapse = ", ")),
               "rerun --classify to write a current catalog")
  }
  catalog %>% mutate(segment = coalesce(.data$segment, "unassigned_at_genus"))
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
  catalog <- .read_catalog(args$catalog)
  ctx <- panel_ctx(cfg, args$species_tree_dir)
  axis <- species_columns(catalog$species, ctx$species_tree, ctx$species_order)
  pages <- highlights_pages(catalog, axis, unlist(cfg$plots$focal_lineages))
  key <- key_page("Highlights",
                  paste("The headline findings of this study, one page each. Hosts are",
                        "columns under the host tree. Every other figure is in its",
                        "stage PDF.",
                        if (!nrow(catalog)) "This study has no loci: no pages follow."),
                  colours = stats::setNames(unname(.tier_colours()),
                                            display_label(names(.tier_colours()))),
                  pages = page_titles(pages))
  save_stage_pdf(c(list(key), pages), args$out_pdf, width = .PAGE_16_9[["width"]],
                 height = .PAGE_16_9[["height"]])
  log_ok("wrote %s, %d pages, %s loci", basename(args$out_pdf), length(pages) + 1L,
         format(nrow(catalog), big.mark = ","))
}


if (sys.nframe() == 0L) run_main(main)
