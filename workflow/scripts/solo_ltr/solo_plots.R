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
# Output is ONE multi-page PDF per genome, and one for the cross-genome view, so a
# reader opens a single file per species rather than hunting through a folder:
#
# Both open with a key page (what the stage does, the fate colours, the page list)
# and follow the house style in plot2sort/style.R (docs/visual_style.md).
#
#   {genome}.solo_ltr.pdf (results/plots/classification/solo_ltr/), one page each:
#     1. funnel             - the subtraction as a waterfall: raw hits to solos.
#     2. identity_by_class  - identity to bait per class, with the cut drawn.
#     3. length_scatter     - candidate length vs identity, with the cut drawn.
#     4. orphan_distance    - distance to nearest orphan, with the pad drawn.
#     5. chromosome_density - candidates per sequence, by fate.
#     6. family_abundance   - solos per seeding element, rank-ordered.
#     7. divergence_age     - divergence from the bait exemplar, as time.
#     8. ltr_tree           - the evidence tree itself, tips coloured by fate.
#     9. tree_enrichment    - same-class-sister observed vs permutation null.
#    10. family_census      - LTR families on the tree, with and without an intact
#                             member, and what each is made of.
#    11. solo_tree          - the tree pruned to its solos, coloured by family.
#    12. family_subtrees    - the largest families of each kind, side by side.
#
#   all_species.solo_ltr.pdf, species on rows beside the host tree:
#     1. solo_intact_ratio  - the headline biological number, per genome.
#     2. class_composition  - the three fates per genome.
#
# Pages 8 to 12 appear only when the tree stage ran (solo_ltr.tree.enable).
#
# Shared infrastructure (style.R, empty_plot, add_titles) is reused
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


# The three fates, as solo_finder writes them. Their colours and words come from
# the shared style (plot2sort/style.R): each fate shares its colour with the
# pipeline tier describing the same situation, so a reader learns them once.
.FATE_LEVELS <- c("solo", "mono_ltr_at_orphan", "intact_flank")

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
funnel_plot <- function(funnel, species) {
  stages <- c("raw_hits", "accepted_hits", "merged_candidates",
              "intact_flank", "mono_ltr_at_orphan", "solo")
  labels <- c(
    raw_hits = "Raw blastn hits",
    accepted_hits = "Pass every criterion",
    merged_candidates = "Merged into loci",
    intact_flank = "Removed: flanks of intact elements",
    mono_ltr_at_orphan = "Removed: monoLTRs at orphans",
    solo = "Solo LTRs"
  )
  d <- funnel[stage %in% stages]
  if (!nrow(d)) return(empty_plot("No funnel data"))
  d[, label := factor(labels[stage], levels = rev(labels[stages]))]
  # The three outcomes wear their fate colours; the filtering steps are grey.
  d[, fill := ifelse(stage %in% names(.FATE_COLOUR), .FATE_COLOUR[stage], .GREY_MID)]

  p <- ggplot(d, aes(x = .data$label, y = .data$count, fill = .data$fill)) +
    geom_col(width = 0.66) +
    geom_text(aes(label = scales::comma(.data$count)), hjust = -0.12, size = 3.3,
              family = .FONT) +
    scale_fill_identity() +
    scale_y_log10(labels = scales::comma, expand = expansion(mult = c(0, 0.2))) +
    coord_flip() +
    theme_retroseek() +
    theme(panel.grid.major.y = element_blank())
  add_titles(
    p,
    title = "From raw LTR matches to solo LTRs",
    subtitle = paste("Every LTR began as one of a pair flanking a provirus, so",
                     "removing the pairable ones and those beside surviving coding",
                     "sequence leaves the solos. Log scale."),
    subset_label = species
  ) + labs(x = NULL, y = "Count (log scale)")
}


#' Identity to the bait arm, per class, with the acceptance threshold drawn.
#'
#' Identity is an age proxy, so this shows the age structure of each class and
#' whether the threshold sits on a real feature of the distribution or cuts
#' arbitrarily through it.
identity_by_class_plot <- function(candidates, species, min_identity) {
  if (!nrow(candidates)) return(empty_plot("No candidates"))
  d <- copy(candidates)
  d[, fate := factor(fate, levels = .FATE_LEVELS)]

  p <- ggplot(d, aes(x = .data$best_identity, fill = .data$fate)) +
    geom_histogram(bins = 40, alpha = 0.75, position = "identity") +
    geom_vline(xintercept = min_identity, linetype = "dashed", colour = .INK_SOFT) +
    scale_fill_manual(values = .FATE_COLOUR, labels = display_label, drop = FALSE) +
    theme_retroseek() +
    theme(legend.position = "bottom")
  add_titles(
    p,
    title = "Identity to the bait LTR, by fate",
    subtitle = sprintf(paste("Dashed line: the %.0f%% acceptance threshold.",
                             "Identity is an age filter, so the left tail is the",
                             "ancient material this method gives up."), min_identity),
    subset_label = species
  ) + labs(x = "Percent identity to the bait arm", y = "Candidate loci", fill = NULL)
}


#' Bait length against hit length, with both cuts drawn.
#'
#' The near-full-length requirement is the criterion that took the solo/intact
#' ratio from 492:1 to 27.6:1, so it earns a panel showing exactly what it removes.
length_scatter_plot <- function(candidates, species, min_hit_length) {
  if (!nrow(candidates)) return(empty_plot("No candidates"))
  d <- copy(candidates)
  d[, fate := factor(fate, levels = .FATE_LEVELS)]

  p <- ggplot(d, aes(x = .data$length, y = .data$best_identity, colour = .data$fate)) +
    geom_point(alpha = 0.35, size = 0.8) +
    geom_vline(xintercept = min_hit_length, linetype = "dashed", colour = .INK_SOFT) +
    scale_x_log10() +
    scale_colour_manual(values = .FATE_COLOUR, labels = display_label, drop = FALSE) +
    # Legend keys at full size and opacity; the points themselves are faint.
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme_retroseek()
  add_titles(
    p,
    title = "Candidate length against identity",
    subtitle = sprintf(paste("Dashed line: the %d bp minimum. Without it the method",
                             "accepts short partial matches and over-reports solos",
                             "by an order of magnitude."), min_hit_length),
    subset_label = species
  ) + labs(x = "Candidate length (bp, log scale)", y = "Percent identity", colour = NULL)
}


#' Distance to the nearest orphan locus, with the proviral pad drawn.
#'
#' This is the criterion only RetroSeek can apply, so it should be visible: the
#' panel shows how many candidates the pad reclassifies and whether the chosen
#' distance sits on a feature of the distribution.
orphan_distance_plot <- function(candidates, species, orphan_pad) {
  d <- candidates[!is.na(orphan_distance) & fate != "intact_flank"]
  if (!nrow(d)) return(empty_plot("No orphan distances"))

  p <- ggplot(d, aes(x = pmax(.data$orphan_distance, 1))) +
    geom_histogram(bins = 50, fill = .TIER_COLOUR[["orphan"]]) +
    geom_vline(xintercept = orphan_pad, linetype = "dashed", colour = .INK_SOFT) +
    scale_x_log10(labels = scales::comma) +
    theme_retroseek()
  add_titles(
    p,
    title = "Distance from each candidate to the nearest orphan locus",
    subtitle = sprintf(paste("Dashed line: the %s bp proviral distance. Candidates",
                             "inside it are monoLTRs beside surviving coding",
                             "sequence, not solos."), format(orphan_pad, big.mark = ",")),
    subset_label = species
  ) + labs(x = "Distance to nearest orphan (bp, log scale)", y = "Candidate loci")
}


#' Solos per chromosome, beside intact element flanks on the same axis.
#'
#' Answers whether solos sit where intact elements sit. They should broadly: a
#' solo marks an integration into the same kind of genomic neighbourhood.
chromosome_density_plot <- function(candidates, species, top_n = 25) {
  if (!nrow(candidates)) return(empty_plot("No candidates"))
  counts <- candidates[, .N, by = .(seqname, fate)]
  keep <- counts[, .(total = sum(N)), by = seqname][order(-total)][seq_len(min(top_n, .N))]
  d <- counts[seqname %in% keep$seqname]
  d[, seqname := factor(seqname, levels = keep$seqname)]
  d[, fate := factor(fate, levels = .FATE_LEVELS)]

  p <- ggplot(d, aes(x = .data$seqname, y = .data$N, fill = .data$fate)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = .FATE_COLOUR, labels = display_label, drop = FALSE) +
    theme_retroseek() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1),
          legend.position = "bottom")
  add_titles(
    p,
    title = "Candidate loci per sequence",
    subtitle = sprintf("The %d sequences carrying the most candidates.", top_n),
    subset_label = species
  ) + labs(x = NULL, y = "Candidate loci", fill = NULL)
}


#' Solos per seeding element, rank-ordered.
#'
#' Shows whether recombination is concentrated in a few prolific families or
#' spread evenly. A steep curve means a handful of families account for most of
#' the solo burden.
family_abundance_plot <- function(candidates, species) {
  d <- candidates[fate == "solo", .N, by = parent][order(-N)]
  if (!nrow(d)) return(empty_plot("No solos"))
  d[, rank := seq_len(.N)]

  p <- ggplot(d, aes(x = .data$rank, y = .data$N)) +
    geom_col(fill = .FATE_COLOUR[["solo"]], width = 1) +
    theme_retroseek()
  add_titles(
    p,
    title = "Solo LTRs per seeding element",
    subtitle = paste("Elements ranked by how many solos their LTR caught.",
                     "A steep curve means recombination is concentrated in a few",
                     "families."),
    subset_label = species
  ) + labs(x = "Seeding element, ranked", y = "Solo LTRs")
}


#' Divergence from the bait exemplar, expressed as time.
#'
#' Deliberately NOT called insertion age. The two-arm molecular clock does not
#' apply to a solo: it has one arm, so there is no internal pair to date. What this
#' measures is how far the solo has drifted from a surviving modern relative, which
#' is a lower bound on its age and is labelled as such.
divergence_age_plot <- function(candidates, species) {
  d <- candidates[fate == "solo"]
  if (!nrow(d)) return(empty_plot("No solos"))
  d[, age_my := age_from_divergence(100 - best_identity) / 1e6]

  p <- ggplot(d, aes(x = .data$age_my)) +
    geom_histogram(bins = 40, fill = .FATE_COLOUR[["solo"]]) +
    theme_retroseek()
  add_titles(
    p,
    title = "Divergence from the bait exemplar, as time",
    subtitle = paste("At 2.2e-9 substitutions/site/year. This is divergence from a",
                     "surviving relative, NOT insertion age: a solo has one arm, so",
                     "the two-arm clock cannot be applied to it."),
    subset_label = species
  ) + labs(x = "Divergence from the bait exemplar (My equivalent)", y = "Solo LTRs")
}


#' The evidence tree itself, drawn from solo_tree_layout.py's coordinates.
#'
#' Tips are points, not labels: with around a thousand tips a label per tip is
#' unreadable, and what the eye needs is where the colours cluster: solo-rich and
#' flank-rich regions of the tree are families with different solo histories.
ltr_tree_plot <- function(tips, segs, summary_dt, species) {
  if (is.null(tips) || !nrow(tips)) return(empty_plot("No tree"))
  class_to_fate <- c(FLANK = "intact_flank", SOLO = "solo", MONO = "mono_ltr_at_orphan")
  d <- copy(tips)
  d[, fate := factor(class_to_fate[class], levels = .FATE_LEVELS)]

  get <- function(key) {
    value <- summary_dt[metric == key, value][1]
    if (length(value) == 0 || is.na(value)) return(NA_real_)
    suppressWarnings(as.numeric(value))
  }
  subtitle <- sprintf(paste(
    "%d tips: both LTR arms of sampled ERV-bearing elements (every sampled solo's",
    "seed among them), sampled solos and monoLTRs.\nSeed control: %.0f%% of solos",
    "sit within 0.1 substitutions/site of the arm that caught them. Arm control:",
    "%.0f%%. Same-class sisters %.0f%% against a %.0f%% permutation null (%.2fx)."),
    nrow(d), 100 * get("solos_near_seed_fraction"),
    100 * get("arm_sisterhood_fraction"),
    100 * get("same_class_sister_observed"),
    100 * get("same_class_sister_null_mean"), get("enrichment"))

  p <- ggplot() +
    { if (!is.null(segs) && nrow(segs)) {
        geom_segment(data = segs, aes(x = .data$x, y = .data$y,
                                      xend = .data$xend, yend = .data$yend),
                     colour = .GREY_MID, linewidth = 0.15)
      } } +
    geom_point(data = d, aes(x = .data$x, y = .data$y, colour = .data$fate),
               size = 0.7) +
    scale_colour_manual(values = .FATE_COLOUR, labels = display_label, drop = FALSE) +
    theme_retroseek_blank()
  add_titles(p, title = "The LTR evidence tree",
             subtitle = subtitle, subset_label = species) +
    labs(colour = NULL)
}


#' The tree's clustering statistic against its permutation null.
#'
#' Without the null this number is uninterpretable, because any structured tree
#' shows some clustering. Both bars are therefore always drawn together.
tree_enrichment_plot <- function(summary_dt, species) {
  get <- function(key) {
    value <- summary_dt[metric == key, value][1]
    if (length(value) == 0 || is.na(value)) return(NA_real_)
    suppressWarnings(as.numeric(value))
  }
  observed <- get("same_class_sister_observed")
  null_mean <- get("same_class_sister_null_mean")
  null_sd <- get("same_class_sister_null_sd")
  control <- get("arm_sisterhood_fraction")
  seed_control <- get("solos_near_seed_fraction")
  if (is.na(observed) || is.na(null_mean)) return(empty_plot("No tree statistics"))

  d <- data.table(
    what = factor(c("Observed", "Label-permuted null"),
                  levels = c("Observed", "Label-permuted null")),
    value = c(observed, null_mean),
    lower = c(observed, null_mean - null_sd),
    upper = c(observed, null_mean + null_sd)
  )
  p <- ggplot(d, aes(x = .data$what, y = .data$value, fill = .data$what)) +
    geom_col(width = 0.55) +
    geom_errorbar(aes(ymin = .data$lower, ymax = .data$upper), width = 0.15) +
    scale_fill_manual(values = c("Observed" = .FATE_COLOUR[["solo"]],
                                 "Label-permuted null" = .GREY_MID)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    theme_retroseek() +
    theme(legend.position = "none")
  add_titles(
    p,
    title = "Do the three fates cluster on the tree?",
    subtitle = sprintf(paste("Tips whose sister group shares their class, against a",
                             "null that permutes the labels on a fixed topology.",
                             "\nControls: %.0f%% of solos sit beside the arm that caught",
                             "them; %.0f%% of elements have their two arms as sisters",
                             "(young bursts of near-identical copies blur this one)."),
                       100 * ifelse(is.na(seed_control), 0, seed_control),
                       100 * ifelse(is.na(control), 0, control)),
    subset_label = species
  ) + labs(x = NULL, y = "Tips with a same-class sister")
}


#' The headline biological number: solos per intact element, per genome.
#'
#' Species on rows in the canonical order (host tree, else config), like every
#' cross-genome figure, so this row sits where the same genome sits elsewhere.
solo_intact_ratio_plot <- function(report, tree = NULL, order = NULL) {
  if (!nrow(report)) return(empty_plot("No per-genome summary"))
  p <- ggplot(report, aes(x = .data$species, y = .data$solo_to_intact_ratio)) +
    # The published mammalian range: solos outnumber intact proviruses by one to
    # two orders of magnitude. A bar far outside it is a red flag, not a finding.
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = 100,
             fill = .TOL_MUTED[["cyan"]], alpha = 0.18) +
    geom_col(fill = .FATE_COLOUR[["solo"]], width = 0.62) +
    geom_text(aes(label = sprintf("%.1f to 1", .data$solo_to_intact_ratio)),
              hjust = -0.15, size = 3.3, family = .FONT) +
    scale_y_log10(expand = expansion(mult = c(0, 0.18))) +
    labs(x = NULL, y = "Solo LTRs per intact locus (log scale)")
  p <- add_titles(
    p,
    title = "Solo LTRs per intact ERV locus",
    subtitle = paste("Shaded band: the published mammalian range, one to two",
                     "orders of magnitude. A high ratio means most of a lineage's",
                     "integrations have had time to recombine away.")
  )
  species_rows(p, report$species, tree = tree, fallback_order = order)
}


#' The three fates as composition per genome, species on rows.
class_composition_plot <- function(all_candidates, tree = NULL, order = NULL) {
  if (!nrow(all_candidates)) return(empty_plot("No candidates"))
  d <- all_candidates[, .N, by = .(species, fate)]
  d[, fate := factor(fate, levels = .FATE_LEVELS)]
  p <- ggplot(d, aes(x = .data$species, y = .data$N, fill = .data$fate)) +
    geom_col(position = "fill", width = 0.66) +
    # Reversed so the legend reads left to right in the same order as the bars.
    scale_fill_manual(values = .FATE_COLOUR, labels = display_label,
                      guide = guide_legend(reverse = TRUE)) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "Share of candidate loci", fill = NULL)
  p <- add_titles(
    p,
    title = "What the LTR matches turn out to be",
    subtitle = paste("Every candidate locus, by fate. The monoLTR class is the one",
                     "LTR_retriever structurally could not separate out.")
  )
  species_rows(p, unique(d$species), tree = tree, fallback_order = order)
}


# Family kinds as tree_families.py writes them, with reader-facing labels.
.KIND_LABELS <- c(no_intact = "No intact member", with_intact = "Has an intact member")


#' LTR families on the tree: which have an intact member, and what they hold.
#'
#' A family is a maximal clade within solo_ltr.tree.family_max_distance (see
#' tree_families.py). The families that matter here are the ones containing
#' solos; those WITHOUT an intact member are LTR families the solo detector
#' reaches and LTRharvest could not.
family_census_plot <- function(families, species, top_n = 30) {
  if (is.null(families)) return(empty_plot("No families with solos"))
  # An unrecognised kind would otherwise render silently as NA, which is exactly
  # how a stale family table from an older run went unnoticed once.
  unknown <- setdiff(unique(families$kind), c(names(.KIND_LABELS), "no_solo"))
  if (length(unknown)) {
    stop(sprintf("unknown family kind(s) %s; regenerate the tree views",
                 paste(unknown, collapse = ", ")))
  }
  d <- families[n_solo > 0]
  if (!nrow(d)) return(empty_plot("No families with solos"))
  n_no_intact <- sum(d$kind == "no_intact")
  solos_no_intact <- sum(d[kind == "no_intact", n_solo])
  subtitle <- sprintf(paste(
    "%d families contain solos; %d of them have no intact member and hold %d of",
    "the %d sampled solos.\nShowing the %d largest. A family is a clade of the",
    "evidence tree whose members are all within the configured distance."),
    nrow(d), n_no_intact, solos_no_intact, sum(d$n_solo), min(top_n, nrow(d)))

  d <- d[order(-n_tips)][seq_len(min(top_n, .N))]
  long <- melt(d, id.vars = c("family", "kind"),
               measure.vars = c("n_solo", "n_mono", "n_flank"),
               variable.name = "what", value.name = "n")
  long[, fate := factor(c(n_solo = "solo", n_mono = "mono_ltr_at_orphan",
                          n_flank = "intact_flank")[as.character(what)],
                        levels = .FATE_LEVELS)]
  long[, family := factor(family, levels = d$family)]
  long[, kind := factor(.KIND_LABELS[kind], levels = .KIND_LABELS)]

  p <- ggplot(long, aes(x = .data$family, y = .data$n, fill = .data$fate)) +
    geom_col() +
    facet_grid(~ kind, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = .FATE_COLOUR, labels = display_label, drop = FALSE) +
    theme_retroseek() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7),
          legend.position = "bottom")
  add_titles(p, title = "LTR families on the evidence tree",
             subtitle = subtitle, subset_label = species) +
    labs(x = "Family", y = "Tips on the tree", fill = NULL)
}


#' The tree pruned to its solos, coloured by family.
#'
#' Pruning keeps the relationships the full tree inferred, so this is the same
#' evidence with everything but the solos removed. Only the largest families get
#' their own colour; a colour per family would be unreadable past a handful.
solo_tree_plot <- function(tips, segs, species, n_colours = 8) {
  if (is.null(tips) || nrow(tips) < 3) return(empty_plot("Too few solos for a tree"))
  d <- copy(tips)
  top <- d[, .N, by = family][order(-N)][seq_len(min(n_colours, .N)), family]
  d[, colour := ifelse(family %in% top, family, "Other")]
  d[, colour := factor(colour, levels = c(top, "Other"))]
  palette <- c(category_colours(top), Other = .GREY_OTHER)

  p <- ggplot() +
    { if (!is.null(segs) && nrow(segs)) {
        geom_segment(data = segs, aes(x = .data$x, y = .data$y,
                                      xend = .data$xend, yend = .data$yend),
                     colour = .GREY_MID, linewidth = 0.2)
      } } +
    geom_point(data = d, aes(x = .data$x, y = .data$y, colour = .data$colour),
               size = 1.1) +
    scale_colour_manual(values = palette) +
    theme_retroseek_blank() +
    theme(legend.position = "right")
  add_titles(p, title = "The solo-only tree",
             subtitle = sprintf(paste(
               "The evidence tree pruned to its %d sampled solos; the largest %d",
               "families coloured.\nSolos sharing a colour belong to one LTR family."),
               nrow(d), length(top)),
             subset_label = species) +
    labs(colour = "Family")
}


#' The largest families of each kind, each as its own small tree.
#'
#' The contrast is the point: a family with no intact member beside one that
#' still has one. Branch lengths are kept, so a long branch is a diverged copy.
family_subtrees_plot <- function(tips, segs, families, species) {
  if (is.null(tips) || !nrow(tips)) return(empty_plot("No families to show"))
  labels <- families[, .(family, panel = sprintf(
    "%s, %s\n%d solos, %d monoLTRs, %d intact flanks",
    family, tolower(.KIND_LABELS[kind]), n_solo, n_mono, n_flank), kind)]
  order_ <- labels[order(kind == "with_intact", family), panel]
  d <- merge(tips, labels, by = "family")
  d[, fate := factor(c(FLANK = "intact_flank", SOLO = "solo",
                       MONO = "mono_ltr_at_orphan")[class], levels = .FATE_LEVELS)]
  d[, panel := factor(panel, levels = order_)]
  sg <- merge(segs, labels, by = "family")
  sg[, panel := factor(panel, levels = order_)]

  p <- ggplot() +
    geom_segment(data = sg, aes(x = .data$x, y = .data$y,
                                xend = .data$xend, yend = .data$yend),
                 colour = .GREY_MID, linewidth = 0.3) +
    geom_point(data = d, aes(x = .data$x, y = .data$y, colour = .data$fate),
               size = 1.6) +
    facet_wrap(~ panel, scales = "free") +
    scale_colour_manual(values = .FATE_COLOUR, labels = display_label, drop = FALSE) +
    theme_retroseek_blank()
  add_titles(p, title = "The largest LTR families, with and without an intact member",
             subtitle = paste("Each panel is one family cut from the evidence tree. \"No",
                              "intact member\" means none among the sampled elements:",
                              "every solo has an intact relative at 95% identity or more."),
             subset_label = species) +
    labs(colour = NULL)
}


# Read an optional table: NULL when the file is absent (tree stage disabled).
read_optional <- function(path) {
  if (file.exists(path)) fread(path) else NULL
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

# The key pages' colour list: the three fates, in words.
.fate_key <- function() stats::setNames(unname(.FATE_COLOUR[.FATE_LEVELS]),
                                        display_label(.FATE_LEVELS))


# Two modes, because Snakemake renders one PDF per genome in parallel and the
# cross-genome PDF once, after all of them:
#   --genome G      write G's multi-page PDF
#   (no --genome)   write the all-species PDF and the summary report
main <- function() {
  script_dir <- .resolve_script_dir()
  source(file.path(script_dir, "..", "plot2sort", "style.R"))  # palette, theme, labels, stage PDFs
  source(file.path(script_dir, "..", "plot2sort", "helpers.R"))
  source(file.path(script_dir, "..", "plot2sort", "io.R"))
  source(file.path(script_dir, "..", "plot2sort", "tree_axis.R"))  # species rows, host tree
  use_retroseek_style()

  parser <- ArgumentParser(description = "Solo-LTR figure panel (ADR-017)")
  parser$add_argument("--table_dir", required = TRUE)
  parser$add_argument("--config", required = TRUE)
  parser$add_argument("--out_pdf", required = TRUE)
  parser$add_argument("--genome", default = NULL,
                      help = "Render this genome's PDF; omit for the summary.")
  parser$add_argument("--report", default = NULL,
                      help = "Summary mode only: where to write the report CSV.")
  parser$add_argument("--species_tree_dir", default = "",
                      help = "Summary mode: the host tree coordinates (species_tree_layout.py).")
  args <- parser$parse_args()

  cfg <- yaml::read_yaml(args$config)
  solo <- cfg$solo_ltr
  species_map <- cfg$species
  table <- function(genome, suffix) file.path(args$table_dir, paste0(genome, suffix))

  if (!is.null(args$genome)) {
    genome <- args$genome
    species <- display_species(genome, species_map)
    funnel <- fread(table(genome, ".funnel.csv"))
    candidates <- fread(table(genome, ".candidates.csv"))
    tree_summary <- read_optional(table(genome, ".tree_summary.csv"))

    plots <- list(
      funnel_plot(funnel, species),
      identity_by_class_plot(candidates, species, solo$min_identity),
      length_scatter_plot(candidates, species, solo$min_hit_length),
      orphan_distance_plot(candidates, species, solo$orphan_pad),
      chromosome_density_plot(candidates, species),
      family_abundance_plot(candidates, species),
      divergence_age_plot(candidates, species)
    )
    if (!is.null(tree_summary)) {
      families <- read_optional(table(genome, ".tree_families.csv"))
      plots <- c(plots, list(
        ltr_tree_plot(read_optional(table(genome, ".tree_tips.csv")),
                      read_optional(table(genome, ".tree_segments.csv")),
                      tree_summary, species),
        tree_enrichment_plot(tree_summary, species),
        family_census_plot(families, species),
        solo_tree_plot(read_optional(table(genome, ".solo_tree_tips.csv")),
                       read_optional(table(genome, ".solo_tree_segments.csv")),
                       species),
        family_subtrees_plot(read_optional(table(genome, ".family_tree_tips.csv")),
                             read_optional(table(genome, ".family_tree_segments.csv")),
                             families, species)
      ))
    }
    key <- key_page(
      sprintf("Solo LTRs in %s", species),
      paste("A solo LTR is what a provirus leaves behind when its two LTRs recombine",
            "and excise everything between them. This stage finds every copy of a",
            "known retroviral LTR and removes the copies that are something else:",
            "flanks of intact elements, and lone LTRs beside surviving coding sequence.",
            "What remains are the solos. The later pages show the LTR evidence tree",
            "and the LTR families cut from it."),
      colours = .fate_key(), pages = page_titles(plots))
    save_stage_pdf(c(list(key), plots), args$out_pdf)
    log_section(sprintf("wrote %s (%d pages)", args$out_pdf, length(plots) + 1L))
    return(invisible(NULL))
  }

  # Summary mode: every genome with a funnel table.
  funnels <- list.files(args$table_dir, pattern = "\\.funnel\\.csv$", full.names = TRUE)
  genomes <- sub("\\.funnel\\.csv$", "", basename(funnels))
  report_rows <- list()
  all_candidates <- list()
  for (genome in genomes) {
    funnel <- fread(table(genome, ".funnel.csv"))
    solos <- funnel[stage == "solo", count][1]
    intact <- funnel[stage == "intact_loci", count][1]
    report_rows[[genome]] <- data.table(
      genome = genome,
      species = display_species(genome, species_map),
      solo = solos,
      mono_ltr_at_orphan = funnel[stage == "mono_ltr_at_orphan", count][1],
      intact_flank = funnel[stage == "intact_flank", count][1],
      intact_loci = intact,
      solo_to_intact_ratio = if (isTRUE(intact > 0)) solos / intact else NA_real_
    )
    candidates <- fread(table(genome, ".candidates.csv"), select = "fate")
    candidates[, species := display_species(genome, species_map)]
    all_candidates[[genome]] <- candidates
  }
  report <- rbindlist(report_rows, fill = TRUE)
  tree <- read_species_tree(args$species_tree_dir)
  order <- display_species(names(species_map), species_map)
  pages <- list(solo_intact_ratio_plot(report, tree, order),
                class_composition_plot(rbindlist(all_candidates, fill = TRUE), tree, order))
  key <- key_page(
    "Solo LTRs across genomes",
    paste("The solo-LTR stage for every genome side by side, genomes on rows in the",
          "host-tree order used by every RetroSeek figure. Per-genome detail is in",
          "each genome's own PDF."),
    colours = .fate_key(), pages = page_titles(pages))
  save_stage_pdf(c(list(key), pages), args$out_pdf,
                 height = page_height_for(nrow(report)))
  fwrite(report, args$report)
  log_section(sprintf("wrote %s and %s (%d genomes)", args$out_pdf, args$report,
                      nrow(report)))
}


# ----------------------------------------------------------------------------
# Entry-point guard - only fire main() under `Rscript solo_plots.R`.
# ----------------------------------------------------------------------------
if (sys.nframe() == 0L) main()
