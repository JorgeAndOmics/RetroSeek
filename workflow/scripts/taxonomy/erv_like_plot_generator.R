# =============================================================================
# erv_like_plot_generator.R
# =============================================================================
# Builds RetroSeek's ERV-like *structural* panel: views of the taxon-founded ERV
# assembly produced by the taxonomy_classify stage. Each valid LTR-element locus
# (one provirus, grouped by its LTR `Parent`) is one row of
# `<genome>.loci.parquet`; this panel summarises their structure - how complete
# they are, whether their genes sit in canonical order, which gene combinations
# occur, their span, and how taxon relates to the genes recovered.
#
# This replaces the retired probe-label-chained erv_like tier: there is now one
# canonical provirus object (the taxon loci table), viewed two ways - the
# taxonomic panel (taxonomy_plot_generator.R) and this structural panel.
#
# Output is ONE PDF (results/plots/classification/structure/structure.pdf): a key
# page, then one page per structure_panel_registry() entry: structural class,
# completeness, gene count, gene combinations, gene order, element length, and
# the lineage by gene heatmap. Pages follow the house style (plot2sort/style.R,
# docs/visual_style.md): hosts on rows beside the host tree, readable names.
#
# Shared infrastructure (style.R, helpers.R, tree_axis.R) is reused from
# plot2sort/. `testthat` and demo_figures.R source this file for its builders;
# the `if (sys.nframe() == 0L) main()` guard keeps the CLI block dormant then.

suppressMessages({
  library(argparse)
  library(arrow)
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
source(file.path(.script_dir, "..", "utils", "log.R"))  # line contract, run_main (ADR-021)
source(file.path(.script_dir, "..", "plot2sort", "style.R"))  # palette, theme, labels, stage PDFs
source(file.path(.script_dir, "..", "plot2sort", "helpers.R"))  # empty_plot, add_titles
source(file.path(.script_dir, "..", "plot2sort", "tree_axis.R"))  # species rows, host tree


# ----------------------------------------------------------------------------
# Load every per-genome taxon-founded loci table, tagging each row with its
# species (filename stem) and coercing the string-typed structural columns the
# classifier emits (completeness / n_main_genes / canonical_order) to numbers.
# Returns one tidy data frame (empty if none).
# ----------------------------------------------------------------------------
load_taxon_loci <- function(input_dir) {
  files <- list.files(input_dir, pattern = "\\.loci\\.parquet$", full.names = TRUE)
  if (length(files) == 0L) return(tibble())
  frames <- lapply(files, function(f) {
    df <- as_tibble(arrow::read_parquet(f))
    if (nrow(df) == 0L) return(NULL)
    df$species <- sub("\\.loci$", "", tools::file_path_sans_ext(basename(f)))
    df
  })
  frames <- Filter(Negate(is.null), frames)
  if (length(frames) == 0L) return(tibble())
  add_structure_companions(bind_rows(frames))
}


# Typed companions the structure builders expect. The loci tables and the
# catalog store every column as a string, so the numeric and logical forms are
# derived rather than parsed, and `span_bp` does not exist on disk at all.
#
# Shared rather than inlined in the loader for the same reason as
# `add_numeric_companions`: taxonomy_segments.R reuses these builders without
# going through this loader, and a missing derived column makes a builder fail
# or fall back to its empty placeholder rather than raising anything useful.
# Empty-safe.
add_structure_companions <- function(df) {
  if (nrow(df) == 0L) return(df)
  df %>%
    mutate(
      completeness    = suppressWarnings(as.numeric(.data$completeness)),  # blank -> NA
      n_main_genes    = suppressWarnings(as.integer(.data$n_main_genes)),  # blank -> NA
      canonical_order = toupper(as.character(.data$canonical_order)) == "TRUE",
      span_bp         = as.numeric(.data$end) - as.numeric(.data$start) + 1
    )
}


# ----------------------------------------------------------------------------
# Structure panel registry: same contract as panel_registry() in
# taxonomy_plot_generator.R, kept here so this file owns its own builders. Both
# main() below and taxonomy_segments.R render it, so the two panels cannot
# drift apart. See that file for the field meanings. Registry order is page
# order in structure.pdf.
#
# All seven read only columns catalog.csv already carries, so the per-segment
# panel can drive them from the catalog slice directly. All are meaningful for
# a single segment: gene content and completeness vary within a genus (measured
# 2026-08-18: genes_present 31 distinct values per segment, n_main_genes 3.1).
# ----------------------------------------------------------------------------
structure_panel_registry <- function() {
  list(
    list(name = "structure_class", build = function(d, ctx) structure_class_plot(d, ctx), data = "loci", segment = TRUE),
    list(name = "completeness", build = function(d, ctx) completeness_plot(d, ctx), data = "loci", segment = TRUE),
    list(name = "n_main_genes", build = function(d, ctx) n_main_genes_plot(d), data = "loci", segment = TRUE),
    list(name = "gene_combinations", build = function(d, ctx) gene_combinations_plot(d), data = "loci", segment = TRUE),
    list(name = "canonical_order", build = function(d, ctx) canonical_order_plot(d, ctx), data = "loci", segment = TRUE),
    list(name = "length_distribution", build = function(d, ctx) length_distribution_plot(d, ctx), data = "loci", segment = TRUE),
    list(name = "composition_heatmap", build = function(d, ctx) composition_heatmap_plot(d), data = "loci", segment = TRUE)
  )
}


# ----------------------------------------------------------------------------
# Plot builders. Each takes the loci frame (and the panel ctx where species are
# rows) and returns a page, or empty_plot() when there is nothing to show.
# Defined at top level so demo_figures.R and the tests can reuse them.
# ----------------------------------------------------------------------------

# Discrete structural class (full / partial / gene) per host: each locus
# committed to one catalogue class (ADR-009). Levels declared in full with
# drop = FALSE so an absent class keeps its slot.
.STRUCTURE_LEVELS <- c("full", "partial", "gene")
structure_class_plot <- function(loci, ctx = NULL) {
  if (nrow(loci) == 0L || !"structure_class" %in% names(loci)) {
    return(empty_plot("No loci"))
  }
  d <- loci %>%
    mutate(structure_class = factor(.data$structure_class, levels = .STRUCTURE_LEVELS)) %>%
    count(.data$species, .data$structure_class, name = "n")
  p <- ggplot(d, aes(x = .data$species, y = .data$n, fill = .data$structure_class)) +
    geom_col(position = position_fill(reverse = TRUE), width = 0.7) +
    scale_fill_manual(values = .STRUCTURE_COLOUR, labels = display_label, drop = FALSE) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "Share of LTR-flanked loci", fill = NULL)
  p <- add_titles(p, "Structural class per host",
                  paste("Full: every main gene present. Partial: some. Single gene:",
                        "one main gene only."))
  on_rows(p, d$species, ctx)
}

# Fraction of main genes present per locus, one row of small multiples per host.
completeness_plot <- function(loci, ctx = NULL) {
  if (nrow(loci) == 0L) return(empty_plot("No loci"))
  d <- loci %>% filter(!is.na(.data$completeness))
  if (nrow(d) == 0L) return(empty_plot("No loci"))
  p <- ggplot(d, aes(x = .data$completeness)) +
    geom_histogram(bins = 20, fill = .DATA_COLOUR, colour = .PAPER, linewidth = 0.2) +
    scale_x_continuous(labels = scales::percent) +
    labs(x = "Main genes present", y = "Loci")
  p <- add_titles(p, "How complete the elements are",
                  "The share of main genes present in each LTR-flanked locus, per host.")
  species_facets(p, ctx)
}

# Number of main genes per locus.
n_main_genes_plot <- function(loci) {
  if (nrow(loci) == 0L) return(empty_plot("No loci"))
  d <- loci %>% filter(!is.na(.data$n_main_genes))
  if (nrow(d) == 0L) return(empty_plot("No loci"))
  counts <- d %>% count(.data$n_main_genes, name = "n")
  p <- ggplot(counts, aes(x = factor(.data$n_main_genes), y = .data$n)) +
    geom_col(fill = .DATA_COLOUR, width = 0.65) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(x = "Main genes per locus", y = "Loci") +
    theme(panel.grid.major.x = element_blank())
  add_titles(p, "Main genes per element",
             "The number of main genes recovered in each LTR-flanked locus.")
}

# Frequency of each gene combination present (e.g. "GAG,POL", "POL"). A genome
# set yields hundreds of rare combinations, so the most common ones get a bar
# each and the rest are pooled into one grey "Other" bar.
.TOP_COMBINATIONS <- 30L
gene_combinations_plot <- function(loci) {
  if (nrow(loci) == 0L) return(empty_plot("No loci"))
  d <- loci %>% filter(nzchar(.data$genes_present))
  if (nrow(d) == 0L) return(empty_plot("No loci"))
  counts <- d %>% count(.data$genes_present, name = "n") %>% arrange(desc(.data$n))
  if (nrow(counts) > .TOP_COMBINATIONS) {
    rest <- counts[-seq_len(.TOP_COMBINATIONS), ]
    counts <- bind_rows(
      counts[seq_len(.TOP_COMBINATIONS), ],
      tibble(genes_present = sprintf("Other (%d combinations)", nrow(rest)),
             n = sum(rest$n)))
  }
  counts <- counts %>%
    mutate(genes_present = factor(.data$genes_present, levels = rev(.data$genes_present)),
           other = grepl("^Other", .data$genes_present))
  p <- ggplot(counts, aes(x = .data$genes_present, y = .data$n, fill = .data$other)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    coord_flip() +
    scale_fill_manual(values = c(`FALSE` = .DATA_COLOUR, `TRUE` = .GREY_OTHER)) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(x = NULL, y = "Loci") +
    # Gene symbols take italics.
    theme(axis.text.y = element_text(face = "italic"),
          panel.grid.major.y = element_blank())
  add_titles(p, "Which genes occur together",
             sprintf(paste("The %d most common combinations of main and diagnostic genes",
                           "found in one locus; the rest pooled."), .TOP_COMBINATIONS))
}

# Canonical versus rearranged main-gene order, per host.
canonical_order_plot <- function(loci, ctx = NULL) {
  if (nrow(loci) == 0L) return(empty_plot("No loci"))
  d <- loci %>%
    mutate(order = ifelse(.data$canonical_order, "canonical", "rearranged")) %>%
    count(.data$species, .data$order, name = "n")
  p <- ggplot(d, aes(x = .data$species, y = .data$n, fill = .data$order)) +
    geom_col(position = position_fill(reverse = TRUE), width = 0.7) +
    scale_fill_manual(values = c(canonical = .GREY_MID, rearranged = .DATA_COLOUR),
                      labels = display_label) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "Share of LTR-flanked loci", fill = NULL)
  p <- add_titles(p, "Gene order",
                  "Main genes in the configured order along the element, or rearranged.")
  on_rows(p, d$species, ctx)
}

# Locus span, one row of small multiples per host.
length_distribution_plot <- function(loci, ctx = NULL) {
  if (nrow(loci) == 0L) return(empty_plot("No loci"))
  d <- loci %>% filter(!is.na(.data$span_bp), .data$span_bp > 0)
  if (nrow(d) == 0L) return(empty_plot("No loci"))
  p <- ggplot(d, aes(x = .data$span_bp)) +
    geom_histogram(bins = 40, fill = .DATA_COLOUR, colour = .PAPER, linewidth = 0.2) +
    scale_x_continuous(labels = scales::label_comma()) +
    labs(x = "Locus span (bp)", y = "Loci")
  p <- add_titles(p, "Element length",
                  "The genomic span of each LTR-flanked locus, per host.")
  species_facets(p, ctx)
}

# Lineage by gene: how often each gene is recovered per lineage.
composition_heatmap_plot <- function(loci) {
  if (nrow(loci) == 0L) return(empty_plot("No taxon-resolved loci"))
  d <- loci %>% filter(.data$resolved == "True", nzchar(.data$genes_present))
  if (nrow(d) == 0L) return(empty_plot("No taxon-resolved loci"))
  long <- d %>%
    separate_rows("genes_present", sep = ",") %>%
    filter(nzchar(.data$genes_present)) %>%
    count(.data$taxon_call, gene = .data$genes_present, name = "n")
  long$ink <- ink_on_ramp(long$n, trans = log10)
  p <- ggplot(long, aes(x = .data$gene, y = .data$taxon_call, fill = .data$n)) +
    geom_tile(colour = .PAPER, linewidth = 0.6) +
    geom_text(aes(label = scales::comma(.data$n), colour = .data$ink), size = 3,
              family = .FONT) +
    scale_colour_identity() +
    scale_fill_ramp(trans = "log10", labels = scales::comma, name = "Loci") +
    scale_y_discrete(labels = taxon_labels) +
    labs(x = NULL, y = NULL) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(face = "italic"))
  add_titles(p, "Genes kept by each lineage",
             "LTR-flanked loci with a confident call, by lineage and by gene recovered.")
}


# ----------------------------------------------------------------------------
# main()
# ----------------------------------------------------------------------------
main <- function() {
  parser <- ArgumentParser(
    description = "Generate RetroSeek ERV-like structural panel (taxon-founded assembly)"
  )
  parser$add_argument("--input", required = TRUE,
                      help = "Directory with per-genome <genome>.loci.parquet tables.")
  parser$add_argument("--out_pdf", required = TRUE,
                      help = "The stage PDF: a key page, then one page per panel entry.")
  parser$add_argument("--config", required = TRUE,
                      help = "YAML config file with plot parameters.")
  parser$add_argument("--species_tree_dir", required = FALSE, default = "",
                      help = "species_tree_layout.py output: the host tree coordinates")
  parser$add_argument("--log", default = NULL,
                      help = "job log file; the Snakemake log: path")
  args <- parser$parse_args()
  log_job(args$log, "erv_like_plot_generator")

  use_retroseek_style()
  cfg <- yaml::read_yaml(args$config)
  log_section(sprintf("RetroSeek structure panel (output: %s)", args$out_pdf))

  loci <- load_taxon_loci(args$input)
  # Readable species names (the config `species:` values) at the single load point.
  if (nrow(loci) > 0L) loci$species <- display_species(loci$species, cfg$species)
  log_section(sprintf("Loaded %d loci across %d species",
                      nrow(loci), length(unique(loci$species))))

  ctx <- panel_ctx(cfg, args$species_tree_dir %||% "")
  pages <- render_panel(structure_panel_registry(), loci, loci, ctx)
  key <- key_page(
    "ERV structure",
    paste("What the LTR-flanked elements are made of: how many of their main genes",
          "survive, in which combinations and order, and how long the elements are.",
          "Hosts are rows in the order of the host tree."),
    colours = stats::setNames(unname(.STRUCTURE_COLOUR), display_label(names(.STRUCTURE_COLOUR))),
    pages = page_titles(pages))
  n_rows <- max(length(ctx$species_order), length(unique(loci$species)))
  save_stage_pdf(c(list(key), pages), args$out_pdf,
                 height = page_height_for(n_rows, per_species = cfg$plots$per_stratum %||% 0.18))
  log_ok("wrote %s, %s pages, %s loci", basename(args$out_pdf),
         format(length(pages) + 1L, big.mark = ","), format(nrow(loci), big.mark = ","))
}


# ----------------------------------------------------------------------------
# Entry-point guard - only fire main() under `Rscript erv_like_plot_generator.R`.
# run_main() logs how the job ended (ADR-021).
# ----------------------------------------------------------------------------
if (sys.nframe() == 0L) run_main(main)
