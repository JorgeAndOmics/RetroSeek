# =============================================================================
# solo_plots.R
# =============================================================================
# The figure panel for native solo-LTR detection (ADR-017).
#
# A solo LTR is what a provirus leaves behind when its two LTRs recombine and
# excise everything between them. The detector finds them by subtraction: blast
# the LTR arms of ERV-bearing elements against the genome, then remove the hits
# that overlap an intact element (those are flanking arms) and the hits sitting
# close to an orphan locus (those are monoLTRs beside surviving coding sequence,
# where the provirus is damaged rather than excised).
#
# Because the method IS a subtraction, the panel's job is to make each subtraction
# and each threshold visible, so a reader can see where candidates died and judge
# whether the thresholds sit in sensible places for this genome.
#
# Per genome (results/plots/classification/solo_ltr/):
#   1. funnel            - the subtraction as a waterfall: raw hits to solos.
#   2. identity_by_class - identity to bait per class, with the cut drawn.
#   3. length_scatter    - bait length vs hit length, with both cuts drawn.
#   4. orphan_distance   - distance to nearest orphan, with the pad drawn.
#   5. chromosome_density- solos per chromosome beside intact elements.
#   6. family_abundance  - solos per seeding element, rank-ordered.
#   7. divergence_age    - divergence from the bait exemplar, as time.
#   8. tree_enrichment   - same-class-sister observed vs permutation null.
#
# Across genomes:
#   9. all_species_solo_intact_ratio - the headline biological number.
#  10. all_species_class_composition - the three fates per genome.
#
# Shared infrastructure (empty_plot, add_titles, save_plot, palettes) is reused
# from plot2sort/*.R. The `if (sys.nframe() == 0L) main()` guard keeps the CLI
# dormant when testthat sources this file for its builders.

suppressMessages({
  library(argparse)
  library(data.table)
  library(ggplot2)
  library(yaml)
})

# Two one-line helpers every R stage script in this project defines for itself
# rather than importing, following the existing convention (see hotspot_detector.R
# and taxonomy_segments.R).
`%||%` <- function(x, y) if (is.null(x)) y else x

log_section <- function(name) {
  message(sprintf("[solo_plots] %s", name))
}


# The three fates, as solo_finder writes them, with one palette used everywhere so
# a colour means the same thing in every panel.
.FATE_LEVELS <- c("solo", "mono_ltr_at_orphan", "intact_flank")
.FATE_LABELS <- c(
  solo = "solo LTR",
  mono_ltr_at_orphan = "monoLTR at an orphan",
  intact_flank = "flank of an intact element"
)
.FATE_FILL <- c(
  solo = "#B03A2E",
  mono_ltr_at_orphan = "#D4AC0D",
  intact_flank = "#5499C7"
)

# Neutral substitution rate for mammals, the same constant the ADR uses to turn
# divergence into time. Two LTRs are identical the day an element inserts, so
# divergence accumulates at twice the per-site rate.
.NEUTRAL_RATE <- 2.2e-9

age_from_divergence <- function(pct) pct / 100 / (2 * .NEUTRAL_RATE)


# ---------------------------------------------------------------------------
# Builders. Each takes plain data.tables and returns a ggplot, so testthat can
# call them without touching the filesystem.
# ---------------------------------------------------------------------------

#' The detection funnel, as a waterfall.
#'
#' The most informative single panel: it shows how 4.5 million raw hits become
#' eleven thousand solos, and which criterion removed each order of magnitude.
funnel_plot <- function(funnel, genome) {
  stages <- c("raw_hits", "accepted_hits", "merged_candidates",
              "intact_flank", "mono_ltr_at_orphan", "solo")
  labels <- c(
    raw_hits = "raw blastn hits",
    accepted_hits = "pass all criteria",
    merged_candidates = "merged into loci",
    intact_flank = "- intact element flanks",
    mono_ltr_at_orphan = "- monoLTRs at orphans",
    solo = "= solo LTRs"
  )
  d <- funnel[stage %in% stages]
  if (!nrow(d)) return(empty_plot("no funnel data"))
  d[, label := factor(labels[stage], levels = labels[stages])]

  p <- ggplot(d, aes(x = .data$label, y = .data$count)) +
    geom_col(fill = "#5499C7", width = 0.65) +
    geom_text(aes(label = format(.data$count, big.mark = ",")),
              hjust = -0.1, size = 3.2) +
    scale_y_log10(expand = expansion(mult = c(0, 0.25))) +
    coord_flip() +
    theme_minimal()
  add_titles(
    p,
    title = "How raw LTR matches become solo-LTR calls",
    subtitle = paste("Log scale. Every LTR began as one of a pair flanking a",
                     "provirus, so removing the pairable ones and the ones beside",
                     "surviving coding sequence leaves the solos."),
    subset_label = genome
  ) + labs(x = NULL, y = "loci (log scale)")
}


#' Identity to the bait arm, per class, with the acceptance threshold drawn.
#'
#' Identity is an age proxy, so this shows the age structure of each class and
#' whether the threshold sits on a real feature of the distribution or cuts
#' arbitrarily through it.
identity_by_class_plot <- function(candidates, genome, min_identity) {
  if (!nrow(candidates)) return(empty_plot("no candidates"))
  d <- copy(candidates)
  d[, fate := factor(fate, levels = .FATE_LEVELS)]

  p <- ggplot(d, aes(x = .data$best_identity, fill = .data$fate)) +
    geom_histogram(bins = 40, alpha = 0.75, position = "identity") +
    geom_vline(xintercept = min_identity, linetype = "dashed", colour = "grey30") +
    scale_fill_manual(values = .FATE_FILL, labels = .FATE_LABELS, drop = FALSE) +
    theme_minimal() +
    theme(legend.position = "bottom")
  add_titles(
    p,
    title = "Identity to the bait LTR, by fate",
    subtitle = sprintf(paste("Dashed line: the %.0f%% acceptance threshold.",
                             "Identity is an age filter, so the left tail is the",
                             "ancient material this method gives up."), min_identity),
    subset_label = genome
  ) + labs(x = "percent identity to the bait arm", y = "candidate loci", fill = NULL)
}


#' Bait length against hit length, with both cuts drawn.
#'
#' The near-full-length requirement is the criterion that took the solo/intact
#' ratio from 492:1 to 27.6:1, so it earns a panel showing exactly what it removes.
length_scatter_plot <- function(candidates, genome, min_hit_length) {
  if (!nrow(candidates)) return(empty_plot("no candidates"))
  d <- copy(candidates)
  d[, fate := factor(fate, levels = .FATE_LEVELS)]

  p <- ggplot(d, aes(x = .data$length, y = .data$best_identity, colour = .data$fate)) +
    geom_point(alpha = 0.35, size = 0.8) +
    geom_vline(xintercept = min_hit_length, linetype = "dashed", colour = "grey30") +
    scale_x_log10() +
    scale_colour_manual(values = .FATE_FILL, labels = .FATE_LABELS, drop = FALSE) +
    theme_minimal() +
    theme(legend.position = "bottom")
  add_titles(
    p,
    title = "Candidate length against identity",
    subtitle = sprintf(paste("Dashed line: the %d bp minimum. Without it the method",
                             "accepts short partial matches and over-reports solos",
                             "by an order of magnitude."), min_hit_length),
    subset_label = genome
  ) + labs(x = "candidate length (bp, log scale)", y = "percent identity", colour = NULL)
}


#' Distance to the nearest orphan locus, with the proviral pad drawn.
#'
#' This is the criterion only RetroSeek can apply, so it should be visible: the
#' panel shows how many candidates the pad reclassifies and whether the chosen
#' distance sits on a feature of the distribution.
orphan_distance_plot <- function(candidates, genome, orphan_pad) {
  d <- candidates[!is.na(orphan_distance) & fate != "intact_flank"]
  if (!nrow(d)) return(empty_plot("no orphan distances"))

  p <- ggplot(d, aes(x = pmax(.data$orphan_distance, 1))) +
    geom_histogram(bins = 50, fill = "#7D3C98", alpha = 0.8) +
    geom_vline(xintercept = orphan_pad, linetype = "dashed", colour = "grey30") +
    scale_x_log10(labels = scales::comma) +
    theme_minimal()
  add_titles(
    p,
    title = "Distance from each candidate to the nearest orphan locus",
    subtitle = sprintf(paste("Dashed line: the %s bp proviral distance. Candidates",
                             "inside it are monoLTRs beside surviving coding",
                             "sequence, not solos."), format(orphan_pad, big.mark = ",")),
    subset_label = genome
  ) + labs(x = "distance to nearest orphan (bp, log scale)", y = "candidate loci")
}


#' Solos per chromosome, beside intact element flanks on the same axis.
#'
#' Answers whether solos sit where intact elements sit. They should broadly: a
#' solo marks an integration into the same kind of genomic neighbourhood.
chromosome_density_plot <- function(candidates, genome, top_n = 25) {
  if (!nrow(candidates)) return(empty_plot("no candidates"))
  counts <- candidates[, .N, by = .(seqname, fate)]
  keep <- counts[, .(total = sum(N)), by = seqname][order(-total)][seq_len(min(top_n, .N))]
  d <- counts[seqname %in% keep$seqname]
  d[, seqname := factor(seqname, levels = keep$seqname)]
  d[, fate := factor(fate, levels = .FATE_LEVELS)]

  p <- ggplot(d, aes(x = .data$seqname, y = .data$N, fill = .data$fate)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = .FATE_FILL, labels = .FATE_LABELS, drop = FALSE) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1),
          legend.position = "bottom")
  add_titles(
    p,
    title = "Candidate loci per sequence",
    subtitle = sprintf("The %d sequences carrying the most candidates.", top_n),
    subset_label = genome
  ) + labs(x = NULL, y = "candidate loci", fill = NULL)
}


#' Solos per seeding element, rank-ordered.
#'
#' Shows whether recombination is concentrated in a few prolific families or
#' spread evenly. A steep curve means a handful of families account for most of
#' the solo burden.
family_abundance_plot <- function(candidates, genome) {
  d <- candidates[fate == "solo", .N, by = parent][order(-N)]
  if (!nrow(d)) return(empty_plot("no solos"))
  d[, rank := seq_len(.N)]

  p <- ggplot(d, aes(x = .data$rank, y = .data$N)) +
    geom_col(fill = "#B03A2E", width = 1) +
    theme_minimal()
  add_titles(
    p,
    title = "Solo LTRs per seeding element",
    subtitle = paste("Elements ranked by how many solos their LTR caught.",
                     "A steep curve means recombination is concentrated in a few",
                     "families."),
    subset_label = genome
  ) + labs(x = "seeding element, ranked", y = "solo LTRs")
}


#' Divergence from the bait exemplar, expressed as time.
#'
#' Deliberately NOT called insertion age. The two-arm molecular clock does not
#' apply to a solo: it has one arm, so there is no internal pair to date. What this
#' measures is how far the solo has drifted from a surviving modern relative, which
#' is a lower bound on its age and is labelled as such.
divergence_age_plot <- function(candidates, genome) {
  d <- candidates[fate == "solo"]
  if (!nrow(d)) return(empty_plot("no solos"))
  d[, age_my := age_from_divergence(100 - best_identity) / 1e6]

  p <- ggplot(d, aes(x = .data$age_my)) +
    geom_histogram(bins = 40, fill = "#1F618D", alpha = 0.85) +
    theme_minimal()
  add_titles(
    p,
    title = "Divergence from the bait exemplar, as time",
    subtitle = paste("At 2.2e-9 substitutions/site/year. This is divergence from a",
                     "surviving relative, NOT insertion age: a solo has one arm, so",
                     "the two-arm clock cannot be applied to it."),
    subset_label = genome
  ) + labs(x = "divergence from the bait exemplar (My equivalent)", y = "solo LTRs")
}


#' The tree's clustering statistic against its permutation null.
#'
#' Without the null this number is uninterpretable, because any structured tree
#' shows some clustering. Both bars are therefore always drawn together.
tree_enrichment_plot <- function(summary_dt, genome) {
  get <- function(key) {
    value <- summary_dt[metric == key, value][1]
    if (length(value) == 0 || is.na(value)) return(NA_real_)
    suppressWarnings(as.numeric(value))
  }
  observed <- get("same_class_sister_observed")
  null_mean <- get("same_class_sister_null_mean")
  null_sd <- get("same_class_sister_null_sd")
  control <- get("arm_sisterhood_fraction")
  if (is.na(observed) || is.na(null_mean)) return(empty_plot("no tree statistics"))

  d <- data.table(
    what = factor(c("observed", "label-permuted null"),
                  levels = c("observed", "label-permuted null")),
    value = c(observed, null_mean),
    lower = c(observed, null_mean - null_sd),
    upper = c(observed, null_mean + null_sd)
  )
  p <- ggplot(d, aes(x = .data$what, y = .data$value, fill = .data$what)) +
    geom_col(width = 0.55) +
    geom_errorbar(aes(ymin = .data$lower, ymax = .data$upper), width = 0.15) +
    scale_fill_manual(values = c("observed" = "#B03A2E",
                                 "label-permuted null" = "grey70")) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    theme_minimal() +
    theme(legend.position = "none")
  add_titles(
    p,
    title = "Do the three fates form their own clades?",
    subtitle = sprintf(paste("Tips whose sister group shares their class, against a",
                             "null that permutes the labels on a fixed topology.",
                             "\nPositive control: %.0f%% of elements have their two",
                             "arms recovered as sister tips, which they must be."),
                       100 * ifelse(is.na(control), 0, control)),
    subset_label = genome
  ) + labs(x = NULL, y = "tips with a same-class sister")
}


#' The headline biological number: solos per intact element, per genome.
solo_intact_ratio_plot <- function(report) {
  if (!nrow(report)) return(empty_plot("no per-genome summary"))
  d <- report[order(-solo_to_intact_ratio)]
  d[, species := factor(species, levels = species)]

  p <- ggplot(d, aes(x = .data$species, y = .data$solo_to_intact_ratio)) +
    # The published mammalian range: solos outnumber intact proviruses by one to
    # two orders of magnitude. A bar far outside it is a red flag, not a finding.
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = 100,
             alpha = 0.12, fill = "#2E86C1") +
    geom_col(fill = "#B03A2E", width = 0.6) +
    geom_text(aes(label = sprintf("%.1f:1", .data$solo_to_intact_ratio)),
              vjust = -0.5, size = 3.2) +
    scale_y_log10(expand = expansion(mult = c(0, 0.15))) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  add_titles(
    p,
    title = "Solo LTRs per intact ERV locus",
    subtitle = paste("Shaded band: the published mammalian range of one to two",
                     "orders of magnitude. A high ratio means most of a lineage's",
                     "integrations have had time to recombine away.")
  ) + labs(x = NULL, y = "solo : intact (log scale)")
}


#' The three fates as composition per genome.
class_composition_plot <- function(all_candidates) {
  if (!nrow(all_candidates)) return(empty_plot("no candidates"))
  d <- all_candidates[, .N, by = .(species, fate)]
  d[, fate := factor(fate, levels = .FATE_LEVELS)]

  p <- ggplot(d, aes(x = .data$species, y = .data$N, fill = .data$fate)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = .FATE_FILL, labels = .FATE_LABELS, drop = FALSE) +
    scale_y_continuous(labels = scales::percent) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "bottom")
  add_titles(
    p,
    title = "What the LTR matches turn out to be",
    subtitle = paste("Every candidate locus, by fate. The monoLTR class is the one",
                     "LTR_retriever structurally cannot separate out.")
  ) + labs(x = NULL, y = "share of candidate loci", fill = NULL)
}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

.resolve_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  getwd()
}

main <- function() {
  script_dir <- .resolve_script_dir()
  source(file.path(script_dir, "..", "plot2sort", "helpers.R"))
  source(file.path(script_dir, "..", "plot2sort", "io.R"))

  parser <- ArgumentParser(description = "Solo-LTR figure panel (ADR-017)")
  parser$add_argument("--table_dir", required = TRUE)
  parser$add_argument("--plot_dir", required = TRUE)
  parser$add_argument("--tree_dir", required = TRUE)
  parser$add_argument("--config", required = TRUE)
  parser$add_argument("--report", required = TRUE)
  args <- parser$parse_args()

  cfg <- yaml::read_yaml(args$config)
  solo <- cfg$solo_ltr
  species_map <- cfg$species
  plot_cfg <- cfg$plots %||% list()
  base_w <- plot_cfg$width %||% 11
  base_h <- plot_cfg$height %||% 7
  dpi <- plot_cfg$dpi %||% 150

  dir.create(args$plot_dir, recursive = TRUE, showWarnings = FALSE)

  funnels <- list.files(args$table_dir, pattern = "\\.funnel\\.csv$", full.names = TRUE)
  if (!length(funnels)) {
    log_section("No solo-LTR tables found; nothing to plot.")
    fwrite(data.table(), args$report)
    return(invisible(NULL))
  }
  genomes <- sub("\\.funnel\\.csv$", "", basename(funnels))

  emit <- function(name, plot) {
    save_plot(name, plot, args$plot_dir, base_w = base_w, base_h = base_h, dpi = dpi)
  }

  report_rows <- list()
  all_candidates <- list()

  for (genome in genomes) {
    funnel <- fread(file.path(args$table_dir, paste0(genome, ".funnel.csv")))
    candidates_path <- file.path(args$table_dir, paste0(genome, ".candidates.csv"))
    candidates <- if (file.exists(candidates_path)) fread(candidates_path) else data.table()

    emit(sprintf("%s_funnel.pdf", genome), funnel_plot(funnel, genome))
    emit(sprintf("%s_identity_by_class.pdf", genome),
         identity_by_class_plot(candidates, genome, solo$min_identity))
    emit(sprintf("%s_length_scatter.pdf", genome),
         length_scatter_plot(candidates, genome, solo$min_hit_length))
    emit(sprintf("%s_orphan_distance.pdf", genome),
         orphan_distance_plot(candidates, genome, solo$orphan_pad))
    emit(sprintf("%s_chromosome_density.pdf", genome),
         chromosome_density_plot(candidates, genome))
    emit(sprintf("%s_family_abundance.pdf", genome),
         family_abundance_plot(candidates, genome))
    emit(sprintf("%s_divergence_age.pdf", genome),
         divergence_age_plot(candidates, genome))

    tree_summary_path <- file.path(args$table_dir, paste0(genome, ".tree_summary.csv"))
    if (file.exists(tree_summary_path)) {
      emit(sprintf("%s_tree_enrichment.pdf", genome),
           tree_enrichment_plot(fread(tree_summary_path), genome))
    }

    solos <- funnel[stage == "solo", count][1]
    intact <- funnel[stage == "intact_loci", count][1]
    report_rows[[genome]] <- data.table(
      genome = genome,
      species = relabel_species(genome, species_map),
      solo = solos,
      mono_ltr_at_orphan = funnel[stage == "mono_ltr_at_orphan", count][1],
      intact_flank = funnel[stage == "intact_flank", count][1],
      intact_loci = intact,
      solo_to_intact_ratio = if (isTRUE(intact > 0)) solos / intact else NA_real_
    )
    if (nrow(candidates)) {
      candidates[, species := relabel_species(genome, species_map)]
      all_candidates[[genome]] <- candidates[, .(species, fate)]
    }
  }

  report <- rbindlist(report_rows, fill = TRUE)
  emit("all_species_solo_intact_ratio.pdf", solo_intact_ratio_plot(report))
  emit("all_species_class_composition.pdf",
       class_composition_plot(rbindlist(all_candidates, fill = TRUE)))

  fwrite(report, args$report)
  # Count what is actually on disk rather than predicting it: the per-genome
  # tree panel is conditional, so an arithmetic guess would drift.
  log_section(sprintf("Done - %d genomes, %d PDFs in %s", nrow(report),
                      length(list.files(args$plot_dir, pattern = "\\.pdf$")),
                      args$plot_dir))
}


# ----------------------------------------------------------------------------
# Entry-point guard - only fire main() under `Rscript solo_plots.R`.
# ----------------------------------------------------------------------------
if (sys.nframe() == 0L) main()
