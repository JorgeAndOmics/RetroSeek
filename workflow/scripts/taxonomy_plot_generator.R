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


# Load every per-genome fragments table (the recovered non-LTR tier). Same schema
# as the loci tables (source == "fragment"); empty if none.
load_fragments <- function(input_dir) {
  files <- list.files(input_dir, pattern = "\\.fragments\\.parquet$", full.names = TRUE)
  if (length(files) == 0L) return(tibble())
  frames <- lapply(files, function(f) {
    df <- as_tibble(arrow::read_parquet(f))
    if (nrow(df) == 0L) return(NULL)
    df$species <- sub("\\.fragments$", "", tools::file_path_sans_ext(basename(f)))
    df
  })
  frames <- Filter(Negate(is.null), frames)
  if (length(frames) == 0L) return(tibble())
  bind_rows(frames)
}


# Build the tidy classification report from the combined (anchored + fragment)
# loci frame: counts by genus, by confidence, by method, plus mosaic and
# integration totals — split by `source` so the two tiers stay distinguishable.
# Returns a long tibble (source, dimension, level, count); empty-safe.
build_report <- function(combined) {
  if (nrow(combined) == 0L) {
    return(tibble(
      source = character(), dimension = character(),
      level = character(), count = integer()
    ))
  }
  cols <- c("source", "dimension", "level", "count")
  if (!"source" %in% names(combined)) combined$source <- "anchored"
  per_source <- function(df, src) {
    genus <- df %>% filter(.data$rank == "genus") %>%
      count(level = .data$genus_call, name = "count") %>%
      mutate(dimension = "genus")
    conf <- df %>%
      count(level = .data$confidence_tag, name = "count") %>%
      mutate(dimension = "confidence")
    method <- df %>% filter(.data$rank == "genus") %>%
      count(level = .data$method, name = "count") %>%
      mutate(dimension = "method")
    summary <- tibble(
      dimension = "summary",
      level     = c("mosaic", "integrations"),
      count     = c(sum(df$is_mosaic == "True"), nrow(df))
    )
    bind_rows(genus, conf, method, summary) %>%
      mutate(source = src) %>%
      select(all_of(cols))
  }
  combined %>%
    group_split(.data$source) %>%
    lapply(function(df) per_source(df, df$source[1])) %>%
    bind_rows()
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

# Confidence-tag composition per species (HC/LC), faceted by tier so anchored
# loci and recovered fragments are both visible. Reads the combined frame.
confidence_plot <- function(combined) {
  if (nrow(combined) == 0L) return(empty_plot("no classified loci"))
  if (!"source" %in% names(combined)) combined$source <- "anchored"
  counts <- combined %>%
    count(.data$species, .data$source, .data$confidence_tag, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$confidence_tag)) +
    geom_col(position = "fill") +
    facet_wrap(~ .data$source) +
    scale_fill_manual(values = c(HC = "#1B9E77", LC = "#D95F02")) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "fraction of loci", fill = "confidence") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_titles(p, "Call confidence", "High vs low confidence (< confidence_min) per species")
}

# Bucket a per-locus blastx hit count into ordered evidence bands. Pure helper
# (unit-tested): 0 / 1 / 2–5 / 6+. Robust to numeric (non-integer) input.
bucket_evidence <- function(n) {
  n <- as.integer(n)
  out <- dplyr::case_when(
    is.na(n) ~ NA_character_,
    n <= 0L  ~ "0",
    n == 1L  ~ "1",
    n <= 5L  ~ "2–5",
    TRUE     ~ "6+"
  )
  factor(out, levels = c("0", "1", "2–5", "6+"))
}

# Per-locus blastx evidence depth; the zero bin is the candidate-novel-retrovirus
# pile. Pseudo-log x so the long tail of well-supported loci stays readable.
evidence_depth_plot <- function(combined) {
  if (nrow(combined) == 0L) return(empty_plot("no classified loci"))
  d <- combined %>% mutate(is_novel = .data$n_hits == 0L)
  n_novel <- sum(d$is_novel, na.rm = TRUE)
  p <- ggplot(d, aes(x = .data$n_hits, fill = .data$is_novel)) +
    geom_histogram(binwidth = 1, colour = "grey30", linewidth = 0.2) +
    scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
    scale_fill_manual(values = c("FALSE" = "#4575B4", "TRUE" = "#D73027"),
                      labels = c("FALSE" = "has homology", "TRUE" = "novel (0 hits)")) +
    facet_wrap(~ .data$source, scales = "free_y") +
    labs(x = "blastx hits per locus (pseudo-log)", y = "loci", fill = NULL) +
    theme_bw()
  add_titles(p, "Blastx evidence depth",
             sprintf("Per-locus homology — %d novel candidates (0 hits)", n_novel))
}

# Raw confidence distribution split by method, with the HC/LC threshold line.
# A histogram (not a KDE): confidence values pile up at discrete points — most
# at exactly 1.000 — so a density estimate smears mass past the [0,1] domain and
# drops single-value groups (e.g. presence, always 1.0). Faceting by method keeps
# the very different scales legible; the dashed line marks confidence_min.
confidence_density_plot <- function(combined, confidence_min = 0.5) {
  d <- combined %>% filter(!is.na(.data$confidence_num))
  if (nrow(d) == 0L) return(empty_plot("no confidence values"))
  p <- ggplot(d, aes(x = .data$confidence_num, fill = .data$method)) +
    geom_histogram(binwidth = 0.05, boundary = 0, colour = "grey30", linewidth = 0.15) +
    geom_vline(xintercept = confidence_min, linetype = "dashed", colour = "grey20") +
    facet_wrap(~ .data$method, scales = "free_y", ncol = 1) +
    scale_fill_aaas(guide = "none") +
    scale_x_continuous(limits = c(-0.02, 1.02)) +
    labs(x = "call confidence", y = "loci") +
    theme_bw()
  add_titles(p, "Confidence calibration",
             sprintf("Confidence by method; dashed = confidence_min (%.2f)", confidence_min))
}

# Call confidence across blastx evidence-depth buckets — does more homology mean
# a more confident call?
confidence_vs_evidence_plot <- function(combined) {
  d <- combined %>%
    filter(!is.na(.data$confidence_num), !is.na(.data$n_hits)) %>%
    mutate(bucket = bucket_evidence(.data$n_hits))
  if (nrow(d) == 0L) return(empty_plot("no data"))
  p <- ggplot(d, aes(x = .data$bucket, y = .data$confidence_num)) +
    geom_boxplot(fill = "#74ADD1", outlier.size = 0.6) +
    labs(x = "blastx hits per locus", y = "call confidence") +
    theme_bw()
  add_titles(p, "Confidence vs evidence", "Call confidence across blastx hit-count buckets")
}

# Structural completeness by tier — how many main genes each locus carries,
# anchored proviruses vs recovered fragments. (Replaces a novel-vs-classified
# view: 0-hit "novel" loci are essentially absent here — they are domain-validated
# so they have homology — so that comparison was empty. This populated view is
# the useful one: it shows fragments are structurally simpler, mostly single
# markers, while anchored loci carry more of the gag/pol/env complement.)
# Counts are shown as a fraction within each tier so the two tiers' very
# different sizes don't swamp the comparison.
structure_by_tier_plot <- function(combined) {
  if (nrow(combined) == 0L || !"completeness_num" %in% names(combined)) {
    return(empty_plot("no loci"))
  }
  d <- combined %>% filter(!is.na(.data$completeness_num))
  if (nrow(d) == 0L) return(empty_plot("no completeness data"))
  p <- ggplot(d, aes(x = .data$completeness_num, y = after_stat(.data$density),
                     fill = .data$source)) +
    geom_histogram(binwidth = 0.1, boundary = 0, position = "identity",
                   alpha = 0.5, colour = "grey40", linewidth = 0.15) +
    scale_fill_manual(values = c("anchored" = "#1F78B4", "fragment" = "#33A02C")) +
    scale_x_continuous(labels = scales::percent) +
    labs(x = "completeness (fraction of main genes present)",
         y = "within-tier density", fill = "tier") +
    theme_bw()
  add_titles(p, "Structural completeness by tier",
             "Anchored proviruses carry more genes; fragments are mostly single markers")
}

# Yield boost from the fragments tier — loci recovered per tier, per species.
source_yield_plot <- function(combined) {
  if (nrow(combined) == 0L) return(empty_plot("no loci"))
  counts <- combined %>% count(.data$species, .data$source, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$source)) +
    geom_col(position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = c("anchored" = "#1F78B4", "fragment" = "#33A02C")) +
    labs(x = NULL, y = "loci", fill = "tier") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_titles(p, "Anchored vs fragment yield", "Loci recovered per tier, per species")
}

# Genus composition split by tier — surfaces genera present only in the fragment
# tier (the novel-lineage check). Long genus tail folded via collapse_long_tail.
genus_by_source_plot <- function(combined) {
  d <- combined %>% filter(.data$rank == "genus")
  if (nrow(d) == 0L) return(empty_plot("no confident genus calls"))
  d <- collapse_long_tail(d, "genus_call", top_n = 20)
  counts <- d %>% count(.data$species, .data$source, .data$genus_call, name = "n")
  p <- ggplot(counts, aes(x = .data$species, y = .data$n, fill = .data$genus_call)) +
    geom_col() +
    facet_wrap(~ .data$source, scales = "free_x") +
    scale_fill_igv() +
    labs(x = NULL, y = "genus-resolved loci", fill = "genus") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_titles(p, "Genus composition by tier",
             "Genera per species, anchored vs fragment (new-lineage check)")
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
  parser$add_argument("--report_csv", required = TRUE,
                      help = "Output path for the tidy classification report CSV.")
  args <- parser$parse_args()

  cfg <- yaml::read_yaml(args$config)
  plot_dpi       <- cfg$plots$dpi    %||% 300
  plot_height    <- cfg$plots$height %||% 12
  plot_width     <- cfg$plots$width  %||% 15
  confidence_min <- cfg$classification$confidence_min %||% 0.5

  dir.create(args$output, showWarnings = FALSE, recursive = TRUE)
  log_section(sprintf("RetroSeek taxonomy plot generation (output: %s)", args$output))

  loci <- load_loci(args$input)
  fragments <- load_fragments(args$input)
  # Relabel genome stems to the config display names (species: map) for all axes.
  if (nrow(loci) > 0L) loci$species <- relabel_species(loci$species, cfg$species)
  if (nrow(fragments) > 0L) {
    fragments$species <- relabel_species(fragments$species, cfg$species)
  }
  combined <- bind_rows(loci, fragments)
  # numeric companions for the evidence/confidence plots (the loci tables store
  # every column as a string). Guarded so an all-empty input stays well-formed.
  if (nrow(combined) > 0L) {
    combined <- combined %>%
      mutate(
        confidence_num   = suppressWarnings(as.numeric(.data$confidence)),
        n_hits           = suppressWarnings(as.integer(.data$n_blastx_hits)),
        completeness_num = suppressWarnings(as.numeric(.data$completeness))
      )
  }
  log_section(sprintf("Loaded %d anchored loci + %d recovered fragments across %d species",
                      nrow(loci), nrow(fragments), length(unique(combined$species))))

  emit <- function(name, plot) {
    save_plot(name, plot, args$output,
              base_w = plot_width, base_h = plot_height, dpi = plot_dpi)
  }

  # Taxonomy composition plots stay on the anchored loci (the genus-founded
  # assembly); the confidence plot spans both tiers.
  emit("genus_composition.png",     genus_composition_plot(loci))
  emit("rank_resolution.png",       rank_resolution_plot(loci))
  emit("method_mix.png",            method_mix_plot(loci))
  emit("erv_class_composition.png", erv_class_composition_plot(loci))
  emit("mosaic_alluvial.png",       mosaic_alluvial_plot(loci))
  emit("confidence.png",            confidence_plot(combined))

  # Evidence / confidence / fragments panel (spans both tiers).
  emit("evidence_depth.png",        evidence_depth_plot(combined))
  emit("confidence_density.png",    confidence_density_plot(combined, confidence_min))
  emit("confidence_vs_evidence.png", confidence_vs_evidence_plot(combined))
  emit("structure_by_tier.png",     structure_by_tier_plot(combined))
  emit("source_yield.png",          source_yield_plot(combined))
  emit("genus_by_source.png",       genus_by_source_plot(combined))

  # Tidy report: counts by genus / confidence / method + mosaic + integrations,
  # split by tier. Concordant with the plots (same combined frame).
  dir.create(dirname(args$report_csv), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(build_report(combined), args$report_csv)

  log_section(sprintf("Done — wrote 12 PNGs to %s + report %s",
                      args$output, args$report_csv))
}


# ----------------------------------------------------------------------------
# Entry-point guard — only fire main() under `Rscript taxonomy_plot_generator.R`.
# ----------------------------------------------------------------------------
if (sys.nframe() == 0L) main()
