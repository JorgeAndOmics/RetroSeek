# =============================================================================
# stage_plot_generator/plots_ltr_interaction.R — hit ↔ LTR feature interplay
# =============================================================================
# How unreduced tBLASTn loci interact with the LTRharvest/LTRdigest features:
# distance to the nearest retrotransposon, position within the enclosing
# element, probe ↔ Pfam-domain agreement, the per-feature overlap breakdown,
# element length vs recovered hits, and strand concordance. Same builder
# contract as the other plot modules. Reuses plot2sort helpers. All describe the
# NON-reduced original tier; the orchestrator stamps the reduced-state note.

.FEATURE_FILL <- c(domain_overlap = "#1b9e77", inside_retro = "#66a61e",
                   flanking_ltr = "#d95f02", disjoint = "#7570b3")
.CONC_FILL    <- c(concordant = "#1b9e77", discordant = "#d95f02")


# Histogram of distance from each unreduced locus to its nearest LTR
# retrotransposon (log10 of distance+1; 0 → inside). The left spike at 0 is the
# candidate pool; the bulk to the right are homology-only (disjoint) loci.
distance_to_retro_plot <- function(ltr_int_df, subset_label = NULL,
                                   warning_caption = NULL) {
  if (nrow(ltr_int_df) == 0L) return(empty_plot())
  d <- ltr_int_df %>% dplyr::filter(!is.na(distance_to_nearest_retro))
  if (nrow(d) == 0L) return(empty_plot("no LTR retrotransposons"))

  p <- ggplot(d, aes(x = distance_to_nearest_retro + 1)) +
    geom_histogram(bins = 40, fill = "#386cb0", colour = NA, alpha = 0.85) +
    scale_x_log10(labels = scales::label_comma()) +
    scale_y_continuous(labels = scales::label_comma()) +
    theme_minimal() +
    labs(x = "Distance to nearest LTR retrotransposon (bp + 1, log10)",
         y = "Unreduced loci") +
    theme(text = element_text(face = "bold"))
  add_titles(
    p,
    title    = "Locus distance to nearest LTR retrotransposon",
    subtitle = "1 (= distance 0) means the locus lies inside a retrotransposon",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}


# Density of where inside-retrotransposon loci sit along the enclosing element
# (strand-aware 5'→3', 0–1), per probe — a metagene of probe placement.
position_within_provirus_plot <- function(ltr_int_df, subset_label = NULL,
                                          warning_caption = NULL) {
  if (nrow(ltr_int_df) == 0L) return(empty_plot())
  d <- ltr_int_df %>% dplyr::filter(!is.na(relative_position_in_retro))
  if (nrow(d) < 2L) return(empty_plot("too few inside-retrotransposon loci"))
  probes <- sort(unique(d$probe))

  p <- ggplot(d, aes(x = relative_position_in_retro, colour = probe, fill = probe)) +
    geom_density(alpha = 0.25, adjust = 1.2) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_colour_manual(values = futurama_unlimited_palette(12, length(probes))) +
    scale_fill_manual(values = futurama_unlimited_palette(12, length(probes))) +
    theme_minimal() +
    labs(x = "Relative position within provirus (5'→3')", y = "Density",
         colour = "Probe", fill = "Probe") +
    theme(text = element_text(face = "bold"))
  add_titles(
    p,
    title    = "Probe position within enclosing provirus",
    subtitle = "Strand-aware metagene of where each probe lands inside the retrotransposon",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}


# Stacked bars, per probe, of the per-LTR-feature overlap class:
# domain_overlap > inside_retro > flanking_ltr > disjoint.
ltr_feature_breakdown_plot <- function(ltr_int_df, subset_label = NULL,
                                       warning_caption = NULL) {
  if (nrow(ltr_int_df) == 0L) return(empty_plot())
  feat_levels <- c("domain_overlap", "inside_retro", "flanking_ltr", "disjoint")
  d <- ltr_int_df %>%
    dplyr::mutate(feature_class = factor(feature_class, levels = feat_levels)) %>%
    dplyr::count(probe, feature_class, name = "count")
  ordered_probe <- order_by_count(d, "probe", weight = "count")
  d <- d %>% dplyr::mutate(probe = factor(probe, levels = ordered_probe))

  p <- ggplot(d, aes(x = probe, y = count, fill = feature_class)) +
    geom_col(colour = "black", linewidth = 0.2) +
    scale_fill_manual(values = .FEATURE_FILL, drop = FALSE) +
    scale_y_continuous(labels = scales::label_comma()) +
    theme_minimal() +
    labs(x = "Probe", y = "Unreduced loci", fill = "LTR feature") +
    theme(text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  out <- add_titles(
    p,
    title    = "Locus overlap with LTR features, per probe",
    subtitle = "domain_overlap = overlaps a probe-assigned Pfam domain (the validation signal)",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
  attr(out, "intended_dims") <- auto_dims(length(ordered_probe), axis = "x")
  out
}


# Strand concordance of inside-retrotransposon loci vs the enclosing element,
# per probe. Discordant hits (opposite strand) flag likely spurious overlaps.
strand_concordance_plot <- function(ltr_int_df, subset_label = NULL,
                                    warning_caption = NULL) {
  if (nrow(ltr_int_df) == 0L) return(empty_plot())
  d <- ltr_int_df %>% dplyr::filter(!is.na(strand_concordant))
  if (nrow(d) == 0L) return(empty_plot("no inside-retrotransposon loci"))
  d <- d %>%
    dplyr::mutate(concordance = factor(ifelse(strand_concordant,
                                              "concordant", "discordant"),
                                       levels = c("concordant", "discordant"))) %>%
    dplyr::count(probe, concordance, name = "count")
  ordered_probe <- order_by_count(d, "probe", weight = "count")
  d <- d %>% dplyr::mutate(probe = factor(probe, levels = ordered_probe))

  p <- ggplot(d, aes(x = probe, y = count, fill = concordance)) +
    geom_col(colour = "black", linewidth = 0.2) +
    scale_fill_manual(values = .CONC_FILL) +
    scale_y_continuous(labels = scales::label_comma()) +
    theme_minimal() +
    labs(x = "Probe", y = "Inside-retrotransposon loci", fill = "Strand") +
    theme(text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  out <- add_titles(
    p,
    title    = "Hit-vs-element strand concordance",
    subtitle = "Among loci inside a retrotransposon; discordant = opposite strand",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
  attr(out, "intended_dims") <- auto_dims(length(ordered_probe), axis = "x")
  out
}


# Tile heatmap of tBLASTn-locus probe (hit_probe) vs the probe assigned to the
# overlapped Pfam domain (domain_probe). The diagonal is agreement (the
# validation signal); off-diagonal cells are cross-probe overlaps.
probe_domain_heatmap <- function(probe_domain_df, subset_label = NULL,
                                 warning_caption = NULL) {
  if (nrow(probe_domain_df) == 0L) return(empty_plot())
  d <- probe_domain_df %>% dplyr::count(hit_probe, domain_probe, name = "count")
  ordered_hit <- order_by_count(d, "hit_probe", weight = "count")
  ordered_dom <- order_by_count(d, "domain_probe", weight = "count")
  d <- d %>% dplyr::mutate(
    hit_probe    = factor(hit_probe, levels = ordered_hit),
    domain_probe = factor(domain_probe, levels = ordered_dom)
  )

  p <- ggplot(d, aes(x = domain_probe, y = hit_probe, fill = count)) +
    geom_tile(colour = "grey90") +
    geom_text(aes(label = count), size = 3, fontface = "bold") +
    scale_fill_gradient(low = "#deebf7", high = "#08519c") +
    theme_minimal() +
    labs(x = "Pfam-domain probe", y = "tBLASTn-locus probe", fill = "Overlaps") +
    theme(text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  out <- add_titles(
    p,
    title    = "Probe × Pfam-domain overlap",
    subtitle = "Diagonal = locus probe matches the overlapped domain's probe (validation agreement)",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
  attr(out, "intended_dims") <- auto_dims(length(ordered_hit), axis = "y")
  out
}


# Scatter of retrotransposon length vs the number of unreduced loci overlapping
# it (from the ltr_structure table). Longer, hit-rich elements are the
# best-supported candidate proviruses.
retro_length_vs_hits_plot <- function(ltr_df, subset_label = NULL,
                                      warning_caption = NULL) {
  if (nrow(ltr_df) == 0L) return(empty_plot())
  d <- ltr_df %>% dplyr::filter(n_overlapping_hits > 0L)
  if (nrow(d) == 0L) return(empty_plot("no hit-overlapping retrotransposons"))

  p <- ggplot(d, aes(x = width, y = n_overlapping_hits)) +
    geom_point(alpha = 0.4, colour = "#386cb0") +
    scale_x_continuous(labels = scales::label_comma()) +
    scale_y_continuous(labels = scales::label_comma()) +
    theme_minimal() +
    labs(x = "Retrotransposon length (bp)", y = "Overlapping unreduced loci") +
    theme(text = element_text(face = "bold"))
  add_titles(
    p,
    title    = "Retrotransposon length vs overlapping hits",
    subtitle = "Per LTRdigest element with at least one overlapping tBLASTn locus",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}
