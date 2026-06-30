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
# it is measured against for step retention). The pipeline has TWO reduction
# branches off `first_reduced_ranges` (gr_virus): the anchored spine
# (candidate -> valid -> loci) descends from gr_virus directly, while the
# fragments branch descends from `global_reduced_ranges` (gr_global, a second,
# stronger reduction). So `candidate`'s parent is `first_reduced_ranges`, NOT
# `global_reduced_ranges` — and `global_reduced_ranges` is the head of the
# fragments branch, a SIBLING of `candidate`, not a step in the anchored spine.
# Getting this wrong makes candidate/global > 100% and a non-monotonic funnel.
# `branch` tags let the plot separate the anchored spine, the fragments branch,
# and the classification tier (the last is grouping/quality, not attrition).
# ----------------------------------------------------------------------------
.STAGE_SPEC <- tibble::tribble(
  ~metric,                ~stage_order, ~branch,          ~parent,                 ~label,
  "raw_blast_hits",        1L,          "main",           NA_character_,           "raw tBLASTn hits",
  "filtered_blast_hits",   2L,          "main",           "raw_blast_hits",        "quality-filtered",
  "first_reduced_ranges",  3L,          "main",           "filtered_blast_hits",   "first reduction",
  "candidate_ranges",      4L,          "main",           "first_reduced_ranges",  "LTR-overlapping (candidate)",
  "valid_ranges",          5L,          "main",           "candidate_ranges",      "domain-validated (valid)",
  "global_reduced_ranges", 6L,          "fragments",      "first_reduced_ranges",  "global reduction",
  "unanchored_fragments",  7L,          "fragments",      "global_reduced_ranges", "non-LTR fragments",
  "fragments_recovered",   8L,          "fragments",      "unanchored_fragments",  "fragments recovered",
  "loci_total",            9L,          "classification", "valid_ranges",          "anchored loci (grouped)",
  "loci_classified",      10L,          "classification", "loci_total",            "loci classified",
  "loci_no_blastx_hit",   11L,          "classification", "loci_total",            "loci w/ no blastx hit"
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


# Per-genome novel-candidate burden from the funnel: how many loci had zero
# blastx homology (loci_no_blastx_hit) and what fraction of all loci that is.
# Empty-safe; guards a zero/absent total so frac is 0 (never NaN).
novel_burden_table <- function(funnel) {
  empty <- tibble::tibble(
    genome = character(), n_novel = double(), n_total = double(), frac = double()
  )
  if (nrow(funnel) == 0L) return(empty)
  wide <- funnel %>%
    dplyr::filter(.data$metric %in% c("loci_total", "loci_no_blastx_hit")) %>%
    dplyr::select("genome", "metric", "value") %>%
    tidyr::pivot_wider(names_from = "metric", values_from = "value")
  if (nrow(wide) == 0L) return(empty)
  if (!"loci_no_blastx_hit" %in% names(wide)) wide$loci_no_blastx_hit <- 0
  if (!"loci_total" %in% names(wide)) wide$loci_total <- 0
  wide %>%
    dplyr::transmute(
      genome  = .data$genome,
      n_novel = dplyr::coalesce(.data$loci_no_blastx_hit, 0),
      n_total = dplyr::coalesce(.data$loci_total, 0),
      frac    = dplyr::if_else(dplyr::coalesce(.data$loci_total, 0) > 0,
                               dplyr::coalesce(.data$loci_no_blastx_hit, 0) /
                                 .data$loci_total, 0)
    )
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


# Per-step retention heatmap: genome × stage, fill = fraction of the prior stage
# surviving. Restricted to the genuine attrition/reduction steps (the anchored
# spine + the fragments branch) — the classification tier is grouping/quality,
# not retention, so it is excluded to keep one consistent semantic on the scale.
# Every cell is now a true subset/reduction ratio, so all are <= 100%.
step_retention_plot <- function(funnel) {
  d <- funnel %>%
    dplyr::filter(.data$branch %in% c("main", "fragments"), !is.na(.data$step_retained))
  if (nrow(d) == 0L) return(empty_plot("no step-retention data"))
  d <- d %>% dplyr::mutate(label = forcats::fct_reorder(.data$label, .data$stage_order))
  p <- ggplot(d, aes(x = .data$label, y = .data$genome, fill = .data$step_retained)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.0f%%", 100 * .data$step_retained)), size = 2.8) +
    scale_fill_gradient2(low = "#D73027", mid = "#FEE08B", high = "#1A9850",
                         midpoint = 0.5, limits = c(0, 1), labels = scales::percent) +
    labs(x = NULL, y = NULL, fill = "step\nretained") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  add_titles(p, "Per-step retention",
             "Fraction of the prior stage surviving each attrition/reduction step")
}


# Fragment recovery yield: per genome, unanchored fragments vs the subset that
# earned a taxonomic call (recovered), with the recovered fraction annotated.
fragment_recovery_plot <- function(funnel) {
  d <- funnel %>% dplyr::filter(.data$metric %in% c("unanchored_fragments", "fragments_recovered"))
  if (nrow(d) == 0L) return(empty_plot("no fragments"))
  bars <- d %>% dplyr::mutate(metric = factor(
    .data$metric, levels = c("unanchored_fragments", "fragments_recovered"),
    labels = c("unanchored", "recovered")))
  fr <- d %>%
    dplyr::select("genome", "metric", "value") %>%
    tidyr::pivot_wider(names_from = "metric", values_from = "value") %>%
    dplyr::mutate(frac = dplyr::if_else(
      dplyr::coalesce(.data$unanchored_fragments, 0) > 0,
      dplyr::coalesce(.data$fragments_recovered, 0) / .data$unanchored_fragments, 0))
  p <- ggplot(bars, aes(x = .data$genome, y = .data$value, fill = .data$metric)) +
    geom_col(position = position_dodge(width = 0.8)) +
    geom_text(data = fr, inherit.aes = FALSE,
              aes(x = .data$genome, y = .data$fragments_recovered,
                  label = scales::percent(.data$frac, accuracy = 1)),
              vjust = -0.4, size = 2.8) +
    scale_fill_manual(values = c(unanchored = "#9E9E9E", recovered = "#1A9850")) +
    labs(x = NULL, y = "fragments", fill = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_titles(p, "Fragment recovery", "Non-LTR fragments recovered as classified loci")
}


# Novel-candidate burden: per-genome SHARE of loci with zero blastx homology
# (candidate novel retroviruses). Plotted as a fraction (not a raw count) so a
# single novel locus among thousands reads as ~0 rather than filling the panel;
# the absolute count is kept on the bar label. The y-axis is the per-locus
# novelty RATE, which also exposes which genome is relatively more novel-rich.
novel_burden_plot <- function(funnel) {
  bt <- novel_burden_table(funnel)
  if (nrow(bt) == 0L) return(empty_plot("no loci"))
  p <- ggplot(bt, aes(x = .data$genome, y = .data$frac)) +
    geom_col(fill = "#762A83") +
    geom_text(aes(label = sprintf("n=%d (%.2f%%)", as.integer(.data$n_novel),
                                  100 * .data$frac)),
              vjust = -0.4, size = 2.8) +
    scale_y_continuous(labels = scales::percent,
                       expand = expansion(mult = c(0, 0.18))) +
    labs(x = NULL, y = "share of loci with no blastx hit") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_titles(p, "Novel-candidate burden",
             "Anchored loci with zero blastx homology (label = count; bar = share of all loci)")
}


# Loss waterfall: the main-chain survival counts per genome, ordered by stage,
# with the absolute count on each bar. Cleaner single-genome funnel than the
# all-branch overview in loss_funnel.png.
loss_waterfall_plot <- function(funnel) {
  d <- funnel %>% dplyr::filter(.data$branch == "main")
  if (nrow(d) == 0L) return(empty_plot("no funnel"))
  d <- d %>% dplyr::mutate(label = forcats::fct_reorder(.data$label, .data$stage_order))
  p <- ggplot(d, aes(x = .data$label, y = .data$value)) +
    geom_col(fill = "#4575B4") +
    geom_text(aes(label = .data$value), vjust = -0.3, size = 2.6) +
    facet_wrap(~ .data$genome, scales = "free_y") +
    labs(x = NULL, y = "surviving ranges / loci") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  add_titles(p, "Loss waterfall", "Surviving ranges along the main chain (per genome)")
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

  # Relabel genome stems to config display names for the PLOTS only; the written
  # table + novel-candidate files keep the stem key for downstream joins.
  funnel_disp <- funnel
  funnel_disp$genome <- relabel_species(funnel_disp$genome, cfg$species)

  emit <- function(name, plot) {
    save_plot(name, plot, args$plot_dir,
              base_w = plot_width, base_h = plot_height, dpi = plot_dpi)
  }
  emit("loss_funnel.png",      loss_funnel_plot(funnel_disp))
  emit("step_retention.png",   step_retention_plot(funnel_disp))
  emit("fragment_recovery.png", fragment_recovery_plot(funnel_disp))
  emit("novel_burden.png",     novel_burden_plot(funnel_disp))
  emit("loss_waterfall.png",   loss_waterfall_plot(funnel_disp))
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
