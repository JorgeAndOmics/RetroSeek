# =============================================================================
# taxonomy_plot_generator.R
# =============================================================================
# Builds RetroSeek's taxonomy panel: per-locus ERV genus-call composition and
# resolution plots, derived from the per-genome `<genome>.loci.parquet` tables
# that taxonomy_classify_loci.py writes to data/tables/taxonomy_classification/.
# Because every plot is computed from those same loci tables, the panel is
# concordant with the tables by construction.
#
# Plots:
#   1. genus_composition      — per-species stacked counts of confident genus calls.
#   2. rank_resolution        — per-species resolved-rank distribution (genus /
#                               subfamily / family / unclassified).
#   3. method_mix             — how genus calls were made (placement / lca / presence).
#   4. erv_class_composition  — per-species Class I/II/III composition (the
#                               literature anchor: bats+human Class I, mouse Class II).
#   5. mosaic_alluvial        — gene -> genus flows across mosaic loci (ggalluvial).
#
# Shared infrastructure (empty_plot, add_titles, save_plot) is reused from
# plot2sort/*.R — not duplicated — so the panel matches the existing plots.
# `testthat` sources this file; the `if (sys.nframe() == 0L) main()` guard keeps
# the CLI block from firing during sourcing.

suppressMessages({
  library(argparse)     # Command-line argument parser
  library(arrow)        # Parquet I/O
  library(tidyverse)    # Data manipulation and visualisation
  library(yaml)         # YAML config
  library(ggsci)        # Scientific colour palettes
  library(ggalluvial)   # Alluvial flows for the mosaic plot
})


# ----------------------------------------------------------------------------
# Locate sibling scripts + source shared modules (same idiom as the other
# plot generators, so the panel reuses one set of helpers).
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
# Pipeline instrumentation — same idiom as the other plot generators (save_plot
# calls log_section, so it must exist in the global env at call time).
# ----------------------------------------------------------------------------
.t0 <- Sys.time()
log_section <- function(name) {
  elapsed <- as.numeric(difftime(Sys.time(), .t0, units = "secs"))
  message(sprintf("[%6.2fs] > %s", elapsed, name))
}


# ----------------------------------------------------------------------------
# Load every per-genome loci table in `input_dir`, tagging each row with its
# species (the filename stem). Returns one tidy data frame (empty if none).
# ----------------------------------------------------------------------------
load_loci <- function(input_dir) {
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
  bind_rows(frames)
}


# ----------------------------------------------------------------------------
# Plot builders. Each takes the combined loci frame and returns a ggplot (or
# empty_plot() when there is nothing to show), so the orchestrator stays flat.
# ----------------------------------------------------------------------------

# Confident (genus-rank) calls, stacked per species and coloured by genus.
genus_composition_plot <- function(loci) {
  d <- loci %>% filter(.data$rank == "genus")
  if (nrow(d) == 0L) return(empty_plot("no confident genus calls"))
  counts <- d %>% count(.data$species, .data$genus_call, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$genus_call)) +
    geom_col() +
    scale_fill_igv() +
    labs(x = NULL, y = "loci (genus-resolved)", fill = "genus") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_titles(p, "ERV genus composition", "Confident (genus-rank) calls per species")
}

# Resolved-rank distribution per species (how deep the call goes).
rank_resolution_plot <- function(loci) {
  if (nrow(loci) == 0L) return(empty_plot("no loci"))
  order_lvl <- c("genus", "subfamily", "family", "none")
  d <- loci %>%
    mutate(rank = factor(ifelse(.data$rank %in% order_lvl, .data$rank, "none"),
                         levels = order_lvl)) %>%
    count(.data$species, .data$rank, name = "n")
  p <- ggplot(d, aes(x = .data$species, y = .data$n, fill = .data$rank)) +
    geom_col(position = "fill") +
    scale_fill_npg() +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "fraction of loci", fill = "resolved rank") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_titles(p, "Taxonomic resolution", "Resolved rank per species")
}

# How each genus call was made (placement vs weighted-LCA vs presence).
method_mix_plot <- function(loci) {
  d <- loci %>% filter(.data$rank == "genus")
  if (nrow(d) == 0L) return(empty_plot("no confident genus calls"))
  counts <- d %>% count(.data$species, .data$method, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$method)) +
    geom_col() +
    scale_fill_aaas() +
    labs(x = NULL, y = "genus calls", fill = "method") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_titles(p, "Genus-call method mix", "Placement vs weighted-LCA vs presence")
}

# ERV class (I/II/III) composition per species — the literature anchor.
erv_class_composition_plot <- function(loci) {
  d <- loci %>% filter(nzchar(.data$erv_class))
  if (nrow(d) == 0L) return(empty_plot("no ERV-class assignments"))
  counts <- d %>% count(.data$species, .data$erv_class, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$erv_class)) +
    geom_col(position = "fill") +
    scale_fill_jco() +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "fraction of classified loci", fill = "ERV class") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_titles(p, "ERV class composition",
             "Class I / II / III per species (Jern/Blomberg)")
}

# Gene -> genus flows across mosaic loci (loci whose member genes disagree).
# mosaic_composition is packed as "GENE:Genus;GENE:Genus"; unpack to flows.
mosaic_alluvial_plot <- function(loci) {
  d <- loci %>% filter(.data$is_mosaic == "True", nzchar(.data$mosaic_composition))
  if (nrow(d) == 0L) return(empty_plot("no mosaic loci"))
  flows <- d %>%
    mutate(.locus = row_number()) %>%
    separate_rows("mosaic_composition", sep = ";") %>%
    separate("mosaic_composition", into = c("gene", "genus"),
             sep = ":", fill = "right", extra = "merge") %>%
    filter(nzchar(.data$gene), nzchar(.data$genus))
  if (nrow(flows) == 0L) return(empty_plot("no mosaic loci"))
  counts <- flows %>% count(.data$gene, .data$genus, name = "n")
  p <- ggplot(counts,
              aes(axis1 = .data$gene, axis2 = .data$genus, y = .data$n)) +
    geom_alluvium(aes(fill = .data$genus)) +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("gene", "genus"), expand = c(0.1, 0.1)) +
    scale_fill_igv() +
    labs(y = "mosaic-locus gene calls", fill = "genus") +
    theme_bw()
  add_titles(p, "Mosaic composition", "Per-gene genus calls within mosaic loci")
}


# ----------------------------------------------------------------------------
# main()
# ----------------------------------------------------------------------------
main <- function() {
  parser <- ArgumentParser(
    description = "Generate RetroSeek taxonomy panel (per-locus genus-call plots)"
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
  log_section(sprintf("RetroSeek taxonomy plot generation (output: %s)", args$output))

  loci <- load_loci(args$input)
  log_section(sprintf("Loaded %d classified loci across %d species",
                      nrow(loci), length(unique(loci$species))))

  emit <- function(name, plot) {
    save_plot(name, plot, args$output,
              base_w = plot_width, base_h = plot_height, dpi = plot_dpi)
  }

  emit("genus_composition.png",     genus_composition_plot(loci))
  emit("rank_resolution.png",       rank_resolution_plot(loci))
  emit("method_mix.png",            method_mix_plot(loci))
  emit("erv_class_composition.png", erv_class_composition_plot(loci))
  emit("mosaic_alluvial.png",       mosaic_alluvial_plot(loci))

  log_section(sprintf("Done — wrote 5 PNGs to %s", args$output))
}


# ----------------------------------------------------------------------------
# Entry-point guard — only fire main() under `Rscript taxonomy_plot_generator.R`.
# ----------------------------------------------------------------------------
if (sys.nframe() == 0L) main()
