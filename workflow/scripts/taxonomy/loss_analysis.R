# =============================================================================
# loss_analysis.R
# =============================================================================
# Quantifies how many candidate ERV loci survive each step of the pipeline, in
# one unified funnel that spans BOTH the R ranges-analysis stage and the Python
# blastx-classification stage. The two stages already emit their counts in the
# same tidy (metric, value) shape - ranges_analysis.R writes
# `<genome>.counts.csv`, taxonomy_classify_loci.py writes
# `<genome>.classification_counts.csv` (+ `<genome>.orphans_counts.csv`) - so
# this script just UNIONs them, orders the metrics into a funnel, and computes
# per-step retention.
#
# It also surfaces the investigation deliverable: loci with valid LTR structure
# but ZERO blastx homology (n_blastx_hits == 0) are candidate NOVEL retroviruses,
# exported per genome as `<genome>.novel_candidates.csv`.
#
# Outputs:
#   - loss_analysis.{parquet,csv}      - the per-genome, per-stage funnel.
#   - <genome>.novel_candidates.csv    - no-blastx-hit loci (per genome).
#   - loss.pdf                         - a key page and five loss pages.
#
# `testthat` sources this file; the `if (sys.nframe() == 0L) main()` guard keeps
# the CLI block from firing during sourcing.

suppressMessages({
  library(argparse)
  library(arrow)
  library(tidyverse)
  library(yaml)
})


# ----------------------------------------------------------------------------
# Funnel specification - the ordered stages and each stage's PARENT (the stage
# it is measured against for step retention). The pipeline has TWO reduction
# branches off `first_reduced_ranges` (gr_virus): the LTR-flanked spine
# (candidate -> valid -> loci) descends from gr_virus directly, while the
# orphan branch descends from `global_reduced_ranges` (gr_global, a second,
# stronger reduction). So `candidate`'s parent is `first_reduced_ranges`, NOT
# `global_reduced_ranges` - and `global_reduced_ranges` is the head of the
# orphan branch, a SIBLING of the ltr-flanked spine, not a step in it.
# Getting this wrong makes the ratios exceed 100% and the funnel non-monotonic.
# `branch` tags let the plot separate the LTR-flanked spine, the orphan branch,
# and the classification tier (the last is grouping/quality, not attrition).
# ----------------------------------------------------------------------------
.STAGE_SPEC <- tibble::tribble(
  ~metric,                ~stage_order, ~branch,          ~parent,                 ~label,
  "raw_blast_hits",        1L,          "main",           NA_character_,           "Raw tBLASTn hits",
  "filtered_blast_hits",   2L,          "main",           "raw_blast_hits",        "Quality-filtered hits",
  "first_reduced_ranges",  3L,          "main",           "filtered_blast_hits",   "First reduction",
  "element_hits_ranges",          4L,          "main",           "first_reduced_ranges",  "Overlapping an LTR element",
  "global_reduced_ranges", 5L,          "orphan",         "first_reduced_ranges",  "Global reduction",
  "orphans",               6L,          "orphan",         "global_reduced_ranges", "Orphan hits",
  "orphans_total",         7L,          "orphan",         "orphans",               "Orphan loci",
  "orphans_recovered",     8L,          "orphan",         "orphans_total",         "Orphans recovered",
  "loci_total",            9L,          "classification", "element_hits_ranges",          "LTR-flanked loci",
  "loci_classified",      10L,          "classification", "loci_total",            "Loci classified",
  "loci_no_blastx_hit",   11L,          "classification", "loci_total",            "Loci with no blastx hit"
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


# Loci with valid LTR structure but zero blastx homology - candidate novel
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
# I/O helpers - read every `<genome>.<suffix>.csv` in a directory into one long
# (genome, metric, value) frame. Missing dir / no files => empty frame.
# ----------------------------------------------------------------------------
.read_counts <- function(dir, suffix) {
  if (is.null(dir) || !dir.exists(dir)) return(tibble::tibble())
  pat <- paste0("\\.", suffix, "\\.csv$")
  files <- list.files(dir, pattern = pat, full.names = TRUE)
  frames <- Filter(Negate(is.null), lapply(files, .read_count_file, pat = pat))
  if (length(frames) == 0L) return(tibble::tibble())
  dplyr::bind_rows(frames)
}

# One genome's (genome, metric, value) rows, the genome taken from the file
# name; NULL for an empty file or one without metric/value columns.
.read_count_file <- function(f, pat) {
  df <- readr::read_csv(f, show_col_types = FALSE)
  if (nrow(df) == 0L || !all(c("metric", "value") %in% names(df))) return(NULL)
  df$genome <- sub(pat, "", basename(f))
  df[, c("genome", "metric", "value")]
}


# ----------------------------------------------------------------------------
# The loss pages. Genomes follow the canonical species order (host tree, else
# config): as rows where genome is an axis, as small multiples in that order
# where each genome gets its own funnel. The two branches wear the tier colours
# they lead to: the main chain ends in LTR-flanked loci, the other in orphans.
# ----------------------------------------------------------------------------
# A function, not a constant: this file is sourced without style.R by the table
# tests, and the colours only exist once style.R is loaded.
.branch_colours <- function() {
  c(main           = .TIER_COLOUR[["ltr-flanked"]],
    orphan         = .TIER_COLOUR[["orphan"]],
    classification = .GREY_MID)
}
.BRANCH_LABELS <- c(main = "Main chain, to LTR-flanked loci",
                    orphan = "Orphan branch",
                    classification = "Classification of loci")

# Genomes as a factor in canonical order, first on top (small multiples read
# from the top left).
.genome_factor <- function(genome, ctx) {
  factor(genome, levels = rev(species_order(genome, ctx$species_tree, ctx$species_order)))
}

# Surviving ranges or loci per stage, one small multiple per genome, stages as
# rows in pipeline order.
loss_funnel_plot <- function(funnel, ctx = NULL) {
  if (nrow(funnel) == 0L) return(empty_plot("No counts to plot"))
  d <- funnel %>%
    dplyr::mutate(label = forcats::fct_rev(forcats::fct_reorder(.data$label, .data$stage_order)),
                  genome = .genome_factor(.data$genome, ctx))
  p <- ggplot(d, aes(x = .data$value, y = .data$label, fill = .data$branch)) +
    geom_col(width = 0.7) +
    facet_wrap(~ .data$genome, scales = "free_x") +
    scale_fill_manual(values = .branch_colours(), labels = .BRANCH_LABELS,
                      breaks = names(.BRANCH_LABELS)) +
    scale_x_continuous(labels = scales::label_comma()) +
    labs(x = "Ranges or loci", y = NULL, fill = NULL) +
    theme(panel.grid.major.y = element_blank(),
          strip.text = element_text(face = "bold.italic"))
  add_titles(p, "What survives each step",
             "Ranges and loci left after every step of the pipeline, recovered orphans included.")
}

# Per-step retention: genome by step, the fraction of the previous stage that
# survives. Only true reduction steps (the main chain and the orphan branch);
# classification groups loci rather than removing them, so it would put a second
# meaning on the scale.
step_retention_plot <- function(funnel, ctx = NULL) {
  d <- funnel %>%
    dplyr::filter(.data$branch %in% c("main", "orphan"), !is.na(.data$step_retained))
  if (nrow(d) == 0L) return(empty_plot("No step retention data"))
  d <- d %>% dplyr::mutate(label = forcats::fct_reorder(.data$label, .data$stage_order))
  d$ink <- ink_on_ramp(d$step_retained)
  p <- ggplot(d, aes(x = .data$label, y = .data$genome, fill = .data$step_retained)) +
    geom_tile(colour = .PAPER, linewidth = 0.6) +
    geom_text(aes(label = scales::percent(.data$step_retained, accuracy = 1),
                  colour = .data$ink), size = 3.2, family = .FONT) +
    scale_colour_identity() +
    scale_fill_ramp(limits = c(0, 1), labels = scales::percent, name = "Retained") +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, 14)) +
    labs(x = NULL, y = NULL) +
    theme(panel.grid = element_blank())
  p <- add_titles(p, "How much each step keeps",
                  "The share of the previous stage that survives each reduction step.")
  on_rows(p, d$genome, ctx, axis = "y")
}

# Orphan recovery: clustered orphan loci per genome, and the share that earned a
# lineage call. Both in loci, so the share is a clean gate.
orphan_recovery_plot <- function(funnel, ctx = NULL) {
  d <- funnel %>% dplyr::filter(.data$metric %in% c("orphans_total", "orphans_recovered"))
  if (nrow(d) == 0L) return(empty_plot("No orphans"))
  wide <- d %>%
    dplyr::select("genome", "metric", "value") %>%
    tidyr::pivot_wider(names_from = "metric", values_from = "value") %>%
    dplyr::mutate(
      orphans_total = dplyr::coalesce(.data$orphans_total, 0),
      orphans_recovered = dplyr::coalesce(.data$orphans_recovered, 0),
      share = dplyr::if_else(.data$orphans_total > 0,
                             .data$orphans_recovered / .data$orphans_total, 0))
  p <- ggplot(wide, aes(x = .data$genome)) +
    geom_col(aes(y = .data$orphans_total), fill = .GREY_OTHER, width = 0.7) +
    geom_col(aes(y = .data$orphans_recovered), fill = .TIER_COLOUR[["orphan"]], width = 0.7) +
    geom_text(aes(y = .data$orphans_total,
                  label = sprintf("%s recovered", scales::percent(.data$share, accuracy = 1))),
              hjust = -0.1, size = 3.2, family = .FONT) +
    scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.2))) +
    labs(x = NULL, y = "Orphan loci")
  p <- add_titles(p, "Orphans recovered",
                  "Clustered orphan loci (grey) and those that earned a lineage call (teal).")
  on_rows(p, wide$genome, ctx)
}

# Novel-candidate burden: the share of loci per genome with no blastx homology
# (candidate novel retroviruses). A share rather than a count, so one novel
# locus among thousands reads as close to zero; the count is on the label.
novel_burden_plot <- function(funnel, ctx = NULL) {
  bt <- novel_burden_table(funnel)
  if (nrow(bt) == 0L) return(empty_plot("No loci"))
  p <- ggplot(bt, aes(x = .data$genome, y = .data$frac)) +
    geom_col(fill = .DATA_COLOUR, width = 0.7) +
    geom_text(aes(label = sprintf("%s loci (%.2f%%)", scales::comma(.data$n_novel),
                                  100 * .data$frac)),
              hjust = -0.1, size = 3.2, family = .FONT) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.25))) +
    labs(x = NULL, y = "Share of loci with no blastx hit")
  p <- add_titles(p, "Candidate novel retroviruses",
                  "LTR-flanked loci with no blastx homology at all, as a share of each genome's loci.")
  on_rows(p, bt$genome, ctx)
}

# The main chain alone, one small multiple per genome, with each count printed:
# a cleaner single-genome funnel than the all-branch overview.
loss_waterfall_plot <- function(funnel, ctx = NULL) {
  d <- funnel %>% dplyr::filter(.data$branch == "main")
  if (nrow(d) == 0L) return(empty_plot("No funnel"))
  d <- d %>%
    dplyr::mutate(label = forcats::fct_rev(forcats::fct_reorder(.data$label, .data$stage_order)),
                  genome = .genome_factor(.data$genome, ctx))
  p <- ggplot(d, aes(x = .data$value, y = .data$label)) +
    geom_col(fill = .branch_colours()[["main"]], width = 0.7) +
    geom_text(aes(label = scales::comma(.data$value)), hjust = -0.1, size = 3,
              family = .FONT) +
    facet_wrap(~ .data$genome, scales = "free_x") +
    scale_x_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.3))) +
    labs(x = "Surviving ranges or loci", y = NULL) +
    theme(panel.grid.major.y = element_blank(),
          strip.text = element_text(face = "bold.italic"))
  add_titles(p, "The main chain, step by step",
             "Surviving ranges along the chain that ends in LTR-flanked loci, per genome.")
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
                      help = "Dir with <genome>.classification_counts.csv + .orphans_counts.csv.")
  parser$add_argument("--loci_dir", required = TRUE,
                      help = "Dir with <genome>.loci.parquet (for novel candidates).")
  parser$add_argument("--out_parquet", required = TRUE)
  parser$add_argument("--out_csv", required = TRUE)
  parser$add_argument("--novel_dir", required = TRUE)
  parser$add_argument("--out_pdf", required = TRUE,
                      help = "The stage PDF: a key page, then the loss pages.")
  parser$add_argument("--config", required = TRUE)
  parser$add_argument("--species_tree_dir", default = "",
                      help = "species_tree_layout.py output: the host tree coordinates")
  parser$add_argument("--log", default = NULL,
                      help = "job log file; the Snakemake log: path")
  args <- parser$parse_args()
  log_job(args$log, "loss_analysis")

  use_retroseek_style()
  cfg <- yaml::read_yaml(args$config)

  log_section("Loading stage counts (ranges + classification + orphans)")
  counts_long <- dplyr::bind_rows(
    .read_counts(args$ranges_counts_dir, "counts"),
    .read_counts(args$classification_counts_dir, "classification_counts"),
    .read_counts(args$classification_counts_dir, "orphans_counts")
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
  n_novel <- 0L
  for (f in loci_files) {
    genome <- sub("\\.loci$", "", tools::file_path_sans_ext(basename(f)))
    novel <- pick_novel_candidates(as_tibble(arrow::read_parquet(f)))
    n_novel <- n_novel + nrow(novel)
    readr::write_csv(novel, file.path(args$novel_dir, paste0(genome, ".novel_candidates.csv")))
  }
  log_section(sprintf("Wrote novel candidates for %d genomes", length(loci_files)))

  # Readable species names for the PAGES only; the written table and the
  # novel-candidate files keep the stem key for downstream joins.
  funnel_disp <- funnel
  funnel_disp$genome <- display_species(funnel_disp$genome, cfg$species)
  ctx <- panel_ctx(cfg, args$species_tree_dir %||% "")
  pages <- list(loss_funnel_plot(funnel_disp, ctx), loss_waterfall_plot(funnel_disp, ctx),
                step_retention_plot(funnel_disp, ctx), orphan_recovery_plot(funnel_disp, ctx),
                novel_burden_plot(funnel_disp, ctx))
  key <- key_page(
    "Where candidates are lost",
    paste("Every step from raw tBLASTn hits to classified ERV loci, with how many",
          "ranges survive each one, for the main chain that ends in LTR-flanked",
          "loci and for the orphan branch. Also the loci with no blastx homology,",
          "exported as candidate novel retroviruses."),
    colours = stats::setNames(unname(.branch_colours()), .BRANCH_LABELS[names(.branch_colours())]),
    pages = page_titles(pages))
  n_rows <- max(length(ctx$species_order), length(unique(funnel_disp$genome)))
  save_stage_pdf(c(list(key), pages), args$out_pdf,
                 height = page_height_for(n_rows, per_species = cfg$plots$per_stratum %||% 0.18))
  log_ok("loss funnel over %s genomes, %s novel candidates",
         format(length(unique(funnel$genome)), big.mark = ","), format(n_novel, big.mark = ","))
}


# Source shared plotting helpers (style, empty_plot, add_titles, species order).
# Placed after the function defs so testthat can source this file without a plot
# environment; main() only runs under Rscript.
# Where this script lives, so it can source its siblings. The file name travels
# in `ofile` when the script is source()d (testthat does, several frames deep)
# and in `--file=` when Rscript runs it. The scripts that tests source carry
# this copy; a script cannot source a shared helper before it knows where it
# lives.
.resolve_script_dir <- function() {
  for (frame in rev(sys.frames())) {
    if (!is.null(frame$ofile)) return(dirname(normalizePath(frame$ofile, mustWork = FALSE)))
  }
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) return(dirname(sub("^--file=", "", file_arg[1])))
  "scripts"
}

if (sys.nframe() == 0L) {
  .script_dir <- .resolve_script_dir()
  source(file.path(.script_dir, "..", "utils", "log.R"))  # line contract, run_main (ADR-021)
  source(file.path(.script_dir, "..", "plot2sort", "style.R"))  # palette, theme, labels, stage PDFs
  source(file.path(.script_dir, "..", "plot2sort", "helpers.R"))
  source(file.path(.script_dir, "..", "plot2sort", "tree_axis.R"))
  run_main(main)
}
