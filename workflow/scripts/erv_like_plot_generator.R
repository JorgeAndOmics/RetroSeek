# =============================================================================
# erv_like_plot_generator.R
# =============================================================================
# Builds RetroSeek's ERV-like *structural* panel: views of the taxon-founded ERV
# assembly produced by the taxonomy_classify stage. Each valid LTR-element locus
# (one provirus, grouped by its LTR `Parent`) is one row of
# `<genome>.loci.parquet`; this panel summarises their structure — how complete
# they are, whether their genes sit in canonical order, which gene combinations
# occur, their span, and how taxon relates to the genes recovered.
#
# This replaces the retired probe-label-chained erv_like tier: there is now one
# canonical provirus object (the taxon loci table), viewed two ways — the
# taxonomic panel (taxonomy_plot_generator.R) and this structural panel.
#
# Plots (results/plots/erv-like/):
#   1. erv_like_completeness        — fraction of main genes present per locus.
#   2. erv_like_canonical_order     — canonical vs rearranged main-gene order.
#   3. erv_like_gene_combinations   — frequency of each gene set (e.g. GAG,POL).
#   4. erv_like_length_distribution — locus span (bp).
#   5. erv_like_n_main_genes        — number of main genes per locus.
#   6. erv_like_composition_heatmap — taxon x gene (which genes each taxon keeps).
#
# Shared infrastructure (empty_plot, add_titles, save_plot) is reused from
# plot2sort/*.R. `testthat` and demo_figures.R source this file for its builders;
# the `if (sys.nframe() == 0L) main()` guard keeps the CLI block dormant then.

suppressMessages({
  library(argparse)
  library(arrow)
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
source(file.path(.script_dir, "plot2sort", "helpers.R"))  # empty_plot, add_titles
source(file.path(.script_dir, "plot2sort", "io.R"))       # save_plot


# ----------------------------------------------------------------------------
# Pipeline instrumentation (save_plot calls log_section, so define it first).
# ----------------------------------------------------------------------------
.t0 <- Sys.time()
log_section <- function(name) {
  elapsed <- as.numeric(difftime(Sys.time(), .t0, units = "secs"))
  message(sprintf("[%6.2fs] > %s", elapsed, name))
}


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
  df <- bind_rows(frames)
  df %>%
    mutate(
      completeness    = suppressWarnings(as.numeric(.data$completeness)),
      n_main_genes    = suppressWarnings(as.integer(.data$n_main_genes)),
      canonical_order = toupper(as.character(.data$canonical_order)) == "TRUE",
      span_bp         = suppressWarnings(as.numeric(.data$end) - as.numeric(.data$start) + 1)
    )
}


# ----------------------------------------------------------------------------
# Plot builders. Each takes the loci frame and returns a ggplot (or empty_plot()
# when there is nothing to show). Defined at top level so demo_figures.R and the
# tests can reuse them.
# ----------------------------------------------------------------------------

# Fraction of main genes present per locus (0..1), as a per-species histogram.
completeness_plot <- function(loci) {
  if (nrow(loci) == 0L) return(empty_plot("no loci"))
  d <- loci %>% filter(!is.na(.data$completeness))
  if (nrow(d) == 0L) return(empty_plot("no loci"))
  p <- ggplot(d, aes(x = .data$completeness, fill = .data$species)) +
    geom_histogram(bins = 20, colour = NA, alpha = 0.85, position = "stack") +
    scale_fill_igv() +
    labs(x = "main-gene completeness (fraction present)", y = "loci", fill = "species") +
    theme_bw()
  add_titles(p, "ERV-like completeness",
             "Fraction of main genes present per LTR-element locus")
}

# Canonical vs rearranged main-gene order, counts per species.
canonical_order_plot <- function(loci) {
  if (nrow(loci) == 0L) return(empty_plot("no loci"))
  d <- loci %>%
    mutate(order = ifelse(.data$canonical_order, "canonical", "rearranged")) %>%
    count(.data$species, .data$order, name = "n")
  p <- ggplot(d, aes(x = .data$species, y = .data$n, fill = .data$order)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = c(canonical = "#1b9e77", rearranged = "#d95f02")) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "fraction of loci", fill = "gene order") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_titles(p, "ERV-like canonical gene order",
             "Main genes in main_probes order vs rearranged")
}

# Frequency of each gene combination present (e.g. "GAG,POL", "POL").
gene_combinations_plot <- function(loci) {
  if (nrow(loci) == 0L) return(empty_plot("no loci"))
  d <- loci %>% filter(nzchar(.data$genes_present))
  if (nrow(d) == 0L) return(empty_plot("no loci"))
  counts <- d %>% count(.data$genes_present, name = "n") %>% arrange(desc(.data$n))
  counts <- counts %>%
    mutate(genes_present = factor(.data$genes_present, levels = rev(.data$genes_present)))
  p <- ggplot(counts, aes(x = .data$genes_present, y = .data$n)) +
    geom_col(fill = "#386cb0", colour = "black", linewidth = 0.2) +
    coord_flip() +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(x = "genes present", y = "loci") +
    theme_bw()
  add_titles(p, "ERV-like gene combinations",
             "Which main/diagnostic genes co-occur per locus")
}

# Locus span (bp) distribution, per species.
length_distribution_plot <- function(loci) {
  if (nrow(loci) == 0L) return(empty_plot("no loci"))
  d <- loci %>% filter(!is.na(.data$span_bp), .data$span_bp > 0)
  if (nrow(d) == 0L) return(empty_plot("no loci"))
  p <- ggplot(d, aes(x = .data$span_bp, fill = .data$species)) +
    geom_histogram(bins = 40, colour = NA, alpha = 0.85, position = "stack") +
    scale_x_continuous(labels = scales::label_comma()) +
    scale_fill_igv() +
    labs(x = "locus span (bp)", y = "loci", fill = "species") +
    theme_bw()
  add_titles(p, "ERV-like length distribution",
             "Genomic span of each LTR-element locus")
}

# Number of main genes per locus (small-integer bar).
n_main_genes_plot <- function(loci) {
  if (nrow(loci) == 0L) return(empty_plot("no loci"))
  d <- loci %>% filter(!is.na(.data$n_main_genes))
  if (nrow(d) == 0L) return(empty_plot("no loci"))
  counts <- d %>% count(.data$n_main_genes, name = "n")
  p <- ggplot(counts, aes(x = factor(.data$n_main_genes), y = .data$n)) +
    geom_col(fill = "#386cb0", colour = "black", linewidth = 0.2) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(x = "main genes per locus", y = "loci") +
    theme_bw()
  add_titles(p, "ERV-like main-gene count",
             "Number of main genes recovered per LTR-element locus")
}

# Taxon x gene composition heatmap: how often each gene is recovered per taxon.
# Unpacks genes_present into individual genes; counts (taxon_call, gene) pairs.
composition_heatmap_plot <- function(loci) {
  if (nrow(loci) == 0L) return(empty_plot("no taxon-resolved loci"))
  d <- loci %>% filter(.data$resolved == "True", nzchar(.data$genes_present))
  if (nrow(d) == 0L) return(empty_plot("no taxon-resolved loci"))
  long <- d %>%
    separate_rows("genes_present", sep = ",") %>%
    filter(nzchar(.data$genes_present)) %>%
    count(.data$taxon_call, gene = .data$genes_present, name = "n")
  p <- ggplot(long, aes(x = .data$gene, y = .data$taxon_call, fill = .data$n)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = .data$n), size = 3) +
    scale_fill_viridis_c(trans = "log10") +
    labs(x = "gene", y = "taxon", fill = "loci") +
    theme_bw()
  add_titles(p, "ERV-like taxon x gene composition",
             "Genes recovered per confident taxon call")
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
  log_section(sprintf("RetroSeek erv-like structural plots (output: %s)", args$output))

  loci <- load_taxon_loci(args$input)
  log_section(sprintf("Loaded %d loci across %d species",
                      nrow(loci), length(unique(loci$species))))

  emit <- function(name, plot) {
    save_plot(name, plot, args$output,
              base_w = plot_width, base_h = plot_height, dpi = plot_dpi)
  }

  emit("erv_like_completeness.png",        completeness_plot(loci))
  emit("erv_like_canonical_order.png",     canonical_order_plot(loci))
  emit("erv_like_gene_combinations.png",   gene_combinations_plot(loci))
  emit("erv_like_length_distribution.png", length_distribution_plot(loci))
  emit("erv_like_n_main_genes.png",        n_main_genes_plot(loci))
  emit("erv_like_composition_heatmap.png", composition_heatmap_plot(loci))

  log_section(sprintf("Done — wrote 6 PNGs to %s", args$output))
}


# ----------------------------------------------------------------------------
# Entry-point guard — only fire main() under `Rscript erv_like_plot_generator.R`.
# ----------------------------------------------------------------------------
if (sys.nframe() == 0L) main()
