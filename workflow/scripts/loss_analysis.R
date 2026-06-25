# =============================================================================
# loss_analysis.R
# =============================================================================
# Quantifies how many candidate ERV loci survive each step of the pipeline, in
# one unified funnel that spans BOTH the R ranges-analysis stage and the Python
# blastx-classification stage. The two stages already emit their counts in the
# same tidy (metric, value) shape — ranges_analysis.R writes
# `<genome>.counts.csv`, taxonomy_classify_loci.py writes
# `<genome>.classification_counts.csv` (+ `<genome>.fragments_counts.csv`) — so
# this script just UNIONs them, orders the metrics into a funnel, and computes
# per-step retention.
#
# It also surfaces the investigation deliverable: loci with valid LTR structure
# but ZERO blastx homology (n_blastx_hits == 0) are candidate NOVEL retroviruses,
# exported per genome as `<genome>.novel_candidates.csv`.
#
# Outputs:
#   - loss_analysis.{parquet,csv}      — the per-genome, per-stage funnel.
#   - <genome>.novel_candidates.csv    — no-blastx-hit loci (per genome).
#   - loss_funnel.png                  — stacked per-stage attrition plot.
#
# `testthat` sources this file; the `if (sys.nframe() == 0L) main()` guard keeps
# the CLI block from firing during sourcing.

suppressMessages({
  library(argparse)
  library(arrow)
  library(tidyverse)
  library(yaml)
  library(ggsci)
})


# ----------------------------------------------------------------------------
# Funnel specification — the ordered stages and each stage's PARENT (the stage
# it is measured against for step retention). Reductions (first/global) are
# merges, not data loss; they are kept in the funnel for completeness but their
# step_retained reflects merging, not dropping. Branch tags let the plot and the
# report separate the main chain, the recovered fragments, and the blastx stage.
# ----------------------------------------------------------------------------
.STAGE_SPEC <- tibble::tribble(
  ~metric,                ~stage_order, ~branch,          ~parent,                ~label,
  "raw_blast_hits",        1L,          "main",           NA_character_,          "raw tBLASTn hits",
  "filtered_blast_hits",   2L,          "main",           "raw_blast_hits",       "quality-filtered",
  "first_reduced_ranges",  3L,          "main",           "filtered_blast_hits",  "first reduction",
  "global_reduced_ranges", 4L,          "main",           "first_reduced_ranges", "global reduction",
  "candidate_ranges",      5L,          "main",           "global_reduced_ranges","LTR-overlapping (candidate)",
  "valid_ranges",          6L,          "main",           "candidate_ranges",     "domain-validated (valid)",
  "unanchored_fragments",  7L,          "fragments",      "global_reduced_ranges","non-LTR fragments",
  "fragments_recovered",   8L,          "fragments",      "unanchored_fragments", "fragments recovered",
  "loci_total",            9L,          "classification", "valid_ranges",         "anchored loci",
  "loci_classified",      10L,          "classification", "loci_total",           "loci classified",
  "loci_no_blastx_hit",   11L,          "classification", "loci_total",           "loci w/ no blastx hit"
)


# Build the per-genome funnel from a long (genome, metric, value) counts frame.
# Returns a tibble with one row per (genome, stage) ordered by stage, carrying
# step_retained (value / parent value) and frac_of_input (value / raw hits).
build_loss_funnel <- function(counts_long) {
  if (nrow(counts_long) == 0L) {
    return(tibble::tibble(
      genome = character(), metric = character(), label = character(),
      branch = character(), stage_order = integer(), value = double(),
      parent = character(), step_retained = double(), frac_of_input = double()
    ))
  }
  joined <- dplyr::inner_join(counts_long, .STAGE_SPEC, by = "metric")
  # per-genome lookup: metric -> value, for parent / input references
  lookup <- joined %>%
    dplyr::select("genome", "metric", "value") %>%
    dplyr::distinct()
  value_of <- function(g, m) {
    if (is.na(m)) return(NA_real_)
    v <- lookup$value[lookup$genome == g & lookup$metric == m]
    if (length(v) == 0L) NA_real_ else v[1]
  }
  joined %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      parent_value  = value_of(.data$genome, .data$parent),
      input_value   = value_of(.data$genome, "raw_blast_hits"),
      step_retained = dplyr::if_else(
        is.na(.data$parent_value) | .data$parent_value == 0,
        NA_real_, .data$value / .data$parent_value
      ),
      frac_of_input = dplyr::if_else(
        is.na(.data$input_value) | .data$input_value == 0,
        NA_real_, .data$value / .data$input_value
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-"parent_value", -"input_value") %>%
    dplyr::arrange(.data$genome, .data$stage_order)
}


# Loci with valid LTR structure but zero blastx homology — candidate novel
# retroviruses. Tolerant of an absent column / empty frame.
pick_novel_candidates <- function(loci_df) {
  if (nrow(loci_df) == 0L || !"n_blastx_hits" %in% names(loci_df)) {
    return(loci_df[0, , drop = FALSE])
  }
  loci_df %>% dplyr::filter(as.integer(.data$n_blastx_hits) == 0L)
}


# ----------------------------------------------------------------------------
# I/O helpers — read every `<genome>.<suffix>.csv` in a directory into one long
# (genome, metric, value) frame. Missing dir / no files => empty frame.
# ----------------------------------------------------------------------------
.read_counts <- function(dir, suffix) {
  if (is.null(dir) || !dir.exists(dir)) return(tibble::tibble())
  pat <- paste0("\\.", suffix, "\\.csv$")
  files <- list.files(dir, pattern = pat, full.names = TRUE)
  frames <- lapply(files, function(f) {
    df <- readr::read_csv(f, show_col_types = FALSE)
    if (nrow(df) == 0L || !all(c("metric", "value") %in% names(df))) return(NULL)
    df$genome <- sub(pat, "", basename(f))
    df[, c("genome", "metric", "value")]
  })
  frames <- Filter(Negate(is.null), frames)
  if (length(frames) == 0L) return(tibble::tibble())
  dplyr::bind_rows(frames)
}


# Stacked per-stage attrition plot — value by ordered stage, faceted by genome,
# coloured by branch. empty_plot is reused from plot2sort/helpers.R.
loss_funnel_plot <- function(funnel) {
  if (nrow(funnel) == 0L) return(empty_plot("no counts to plot"))
  d <- funnel %>%
    dplyr::mutate(label = forcats::fct_reorder(.data$label, .data$stage_order))
  p <- ggplot(d, aes(x = .data$label, y = .data$value, fill = .data$branch)) +
    geom_col() +
    facet_wrap(~ .data$genome, scales = "free_y") +
    scale_fill_npg() +
    labs(x = NULL, y = "ranges / loci", fill = "stage") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  add_titles(p, "Pipeline loss funnel",
             "Surviving ranges/loci per stage (incl. recovered fragments)")
}


# ----------------------------------------------------------------------------
# main()
# ----------------------------------------------------------------------------
main <- function() {
  parser <- ArgumentParser(
    description = "Unified per-stage loss funnel + novel-candidate export"
  )
  parser$add_argument("--ranges_counts_dir", required = TRUE,
                      help = "Dir with <genome>.counts.csv (ranges_analysis).")
  parser$add_argument("--classification_counts_dir", required = TRUE,
                      help = "Dir with <genome>.classification_counts.csv + .fragments_counts.csv.")
  parser$add_argument("--loci_dir", required = TRUE,
                      help = "Dir with <genome>.loci.parquet (for novel candidates).")
  parser$add_argument("--out_parquet", required = TRUE)
  parser$add_argument("--out_csv", required = TRUE)
  parser$add_argument("--novel_dir", required = TRUE)
  parser$add_argument("--plot_dir", required = TRUE)
  parser$add_argument("--config", required = TRUE)
  args <- parser$parse_args()

  cfg <- yaml::read_yaml(args$config)
  plot_dpi    <- cfg$plots$dpi    %||% 300
  plot_height <- cfg$plots$height %||% 12
  plot_width  <- cfg$plots$width  %||% 15

  log_section("Loading stage counts (ranges + classification + fragments)")
  counts_long <- dplyr::bind_rows(
    .read_counts(args$ranges_counts_dir, "counts"),
    .read_counts(args$classification_counts_dir, "classification_counts"),
    .read_counts(args$classification_counts_dir, "fragments_counts")
  )
  funnel <- build_loss_funnel(counts_long)
  log_section(sprintf("Funnel: %d stage rows across %d genomes",
                      nrow(funnel), length(unique(funnel$genome))))

  dir.create(dirname(args$out_parquet), showWarnings = FALSE, recursive = TRUE)
  dir.create(dirname(args$out_csv),     showWarnings = FALSE, recursive = TRUE)
  arrow::write_parquet(funnel, args$out_parquet)
  readr::write_csv(funnel, args$out_csv)

  # Per-genome novel candidates from the loci parquet tables.
  dir.create(args$novel_dir, showWarnings = FALSE, recursive = TRUE)
  loci_files <- list.files(args$loci_dir, pattern = "\\.loci\\.parquet$", full.names = TRUE)
  for (f in loci_files) {
    genome <- sub("\\.loci$", "", tools::file_path_sans_ext(basename(f)))
    novel <- pick_novel_candidates(as_tibble(arrow::read_parquet(f)))
    readr::write_csv(novel, file.path(args$novel_dir, paste0(genome, ".novel_candidates.csv")))
  }
  log_section(sprintf("Wrote novel candidates for %d genomes", length(loci_files)))

  save_plot("loss_funnel.png", loss_funnel_plot(funnel), args$plot_dir,
            base_w = plot_width, base_h = plot_height, dpi = plot_dpi)
  log_section("Done")
}


# Source shared plotting helpers (empty_plot, add_titles, save_plot, log_section).
# Placed after the function defs so testthat can source this file without a plot
# environment; main() only runs under Rscript.
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

if (sys.nframe() == 0L) {
  .script_dir <- .resolve_script_dir()
  source(file.path(.script_dir, "plot2sort", "helpers.R"))
  source(file.path(.script_dir, "plot2sort", "io.R"))
  .t0 <- Sys.time()
  log_section <- function(name) {
    elapsed <- as.numeric(difftime(Sys.time(), .t0, units = "secs"))
    message(sprintf("[%6.2fs] > %s", elapsed, name))
  }
  main()
}
