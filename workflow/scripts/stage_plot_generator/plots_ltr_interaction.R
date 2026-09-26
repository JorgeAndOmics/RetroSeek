# =============================================================================
# stage_plot_generator/plots_ltr_interaction.R - hit <-> LTR feature interplay
# =============================================================================
# How unreduced tBLASTn loci interact with the LTRharvest/LTRdigest features:
# distance to the nearest retrotransposon, position within the enclosing
# element, probe <-> Pfam-domain agreement, the per-feature overlap breakdown,
# element length vs recovered hits, and strand concordance. Same builder
# contract as the other plot modules. All describe the NON-reduced original
# tier; the orchestrator stamps the tier note.


# Distance from each unreduced locus to its nearest LTR retrotransposon (log10
# of distance + 1; 0 -> inside). The spike at 1 is the candidate pool; the bulk
# to the right are homology-only (disjoint) loci.
distance_to_retro_plot <- function(ltr_int_df, subset_label = NULL,
                                   warning_caption = NULL) {
  if (nrow(ltr_int_df) == 0L) return(empty_plot())
  d <- ltr_int_df %>% dplyr::filter(!is.na(distance_to_nearest_retro))
  if (nrow(d) == 0L) return(empty_plot("no LTR retrotransposons"))

  p <- ggplot(d, aes(x = distance_to_nearest_retro + 1)) +
    geom_histogram(bins = 40, fill = .DATA_COLOUR, colour = NA) +
    scale_x_log10(labels = scales::label_comma()) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(x = "Distance to the nearest LTR retrotransposon (bp + 1, log scale)",
         y = "Unreduced loci")
  add_titles(
    p,
    title    = "How far each locus lies from an LTR element",
    subtitle = "A distance of 1 means the locus lies inside a retrotransposon.",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}


# Where inside-retrotransposon loci sit along the enclosing element (strand-aware
# 5' to 3', 0 to 1): a metagene of probe placement, one small multiple per probe,
# the probe with the most loci first. Overlaid densities of a dozen probes would
# be unreadable.
position_within_provirus_plot <- function(ltr_int_df, subset_label = NULL,
                                          warning_caption = NULL) {
  if (nrow(ltr_int_df) == 0L) return(empty_plot())
  d <- ltr_int_df %>% dplyr::filter(!is.na(relative_position_in_retro))
  if (nrow(d) < 2L) return(empty_plot("too few inside-retrotransposon loci"))
  # A density needs a few points; probes with fewer are left out, and said so.
  counts <- d %>% dplyr::count(probe, name = "n")
  kept <- counts$probe[counts$n >= 5L]
  d <- d %>% dplyr::filter(probe %in% kept)
  if (nrow(d) == 0L) return(empty_plot("too few inside-retrotransposon loci"))
  d$probe <- factor(d$probe, levels = order_by_count(d, "probe"))

  p <- ggplot(d, aes(x = relative_position_in_retro)) +
    geom_density(fill = .DATA_COLOUR, colour = NA, alpha = 0.85, adjust = 1.2) +
    facet_wrap(~ probe, scales = "free_y") +
    coord_cartesian(xlim = c(0, 1)) +
    scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("5'", "Middle", "3'")) +
    labs(x = "Position along the enclosing element", y = "Density") +
    theme(strip.text = element_text(face = "bold.italic"),
          axis.text.y = element_blank())
  add_titles(
    p,
    title    = "Where each probe lands inside an element",
    subtitle = sprintf(
      paste("Position of each locus along the retrotransposon that contains",
            "it, strand-aware. Probes with fewer than 5 such loci left out (%d)."),
      sum(counts$n < 5L)
    ),
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}


# Per probe, each locus's closest relationship to an LTR element: overlapping a
# domain assigned to the probe, inside an element, flanking an LTR, or none.
# Closer means darker on the ramp.
ltr_feature_breakdown_plot <- function(ltr_int_df, subset_label = NULL,
                                       warning_caption = NULL) {
  if (nrow(ltr_int_df) == 0L) return(empty_plot())
  feat_levels <- c("domain_overlap", "inside_retro", "flanking_ltr", "disjoint")
  d <- ltr_int_df %>%
    dplyr::mutate(feature_class = factor(feature_class, levels = feat_levels)) %>%
    dplyr::count(probe, feature_class, name = "count")
  ordered_probe <- rev(order_by_count(d, "probe", weight = "count"))
  d <- d %>% dplyr::mutate(probe = factor(probe, levels = ordered_probe))

  p <- ggplot(d, aes(x = probe, y = count, fill = feature_class)) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
    scale_fill_manual(values = c(domain_overlap = seq_colours(4)[4],
                                 inside_retro   = seq_colours(4)[3],
                                 flanking_ltr   = seq_colours(4)[2],
                                 disjoint       = .GREY_OTHER),
                      labels = c(domain_overlap = "On a matching Pfam domain",
                                 inside_retro   = "Inside an element",
                                 flanking_ltr   = "Flanking an LTR",
                                 disjoint       = "None"),
                      drop = FALSE) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(x = NULL, y = "Unreduced loci", fill = NULL)
  add_titles(
    .probe_rows(p),
    title    = "How close each probe's loci come to an LTR element",
    subtitle =
      "Overlapping a Pfam domain assigned to the probe is the validation signal.",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}


# Strand agreement of inside-retrotransposon loci with the enclosing element,
# per probe. Opposite-strand loci flag likely spurious overlaps, so they take
# the data colour and agreeing loci the grey.
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
  ordered_probe <- rev(order_by_count(d, "probe", weight = "count"))
  d <- d %>% dplyr::mutate(probe = factor(probe, levels = ordered_probe))

  p <- ggplot(d, aes(x = probe, y = count, fill = concordance)) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
    scale_fill_manual(values = c(concordant = .GREY_MID, discordant = .DATA_COLOUR),
                      labels = c(concordant = "Same strand",
                                 discordant = "Opposite strand")) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(x = NULL, y = "Loci inside an element", fill = NULL)
  add_titles(
    .probe_rows(p),
    title    = "Strand agreement with the enclosing element",
    subtitle = "Among loci inside a retrotransposon; an opposite-strand hit is likely a chance overlap.",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}


# The tBLASTn locus probe (hit_probe) against the Pfam domain it overlaps
# (domain_name, raw). See ADR-016: agreement between the two is the validation
# signal. A genome set overlaps a hundred or so domains, so the most frequent
# get a column each and the rest are pooled; counts span orders of magnitude,
# so the ramp is logarithmic, and cell numbers are drawn only while they fit.
.TOP_DOMAINS <- 25L
probe_domain_heatmap <- function(probe_domain_df, subset_label = NULL,
                                 warning_caption = NULL) {
  if (nrow(probe_domain_df) == 0L) return(empty_plot())
  d <- probe_domain_df %>%
    dplyr::count(hit_probe, domain_name, name = "count") %>%
    collapse_long_tail("domain_name", top_n = .TOP_DOMAINS, weight = "count") %>%
    dplyr::group_by(hit_probe, domain_name) %>%
    dplyr::summarise(count = sum(count), .groups = "drop")
  ordered_hit <- rev(order_by_count(d, "hit_probe", weight = "count"))
  ordered_dom <- order_by_count(d, "domain_name", weight = "count")
  ordered_dom <- c(ordered_dom[!grepl("^Other", ordered_dom)],
                   ordered_dom[grepl("^Other", ordered_dom)])
  d <- d %>% dplyr::mutate(
    hit_probe   = factor(hit_probe, levels = ordered_hit),
    domain_name = factor(domain_name, levels = ordered_dom)
  )
  d$ink <- ink_on_ramp(d$count, trans = log10)

  p <- ggplot(d, aes(x = domain_name, y = hit_probe, fill = count)) +
    geom_tile(colour = .PAPER, linewidth = 0.6) +
    scale_colour_identity() +
    scale_fill_ramp(trans = "log10", labels = scales::label_comma(),
                    name = "Overlaps") +
    labs(x = "Pfam domain", y = "Probe of the tBLASTn locus") +
    # Pfam names are long: the one place a tilted axis label is the lesser evil.
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 40, hjust = 1),
          axis.text.y = element_text(face = "italic"))
  if (length(ordered_dom) <= 15L) {
    p <- p + geom_text(aes(label = scales::comma(count), colour = ink), size = 3,
                       family = .FONT)
  }
  add_titles(
    p,
    title    = "Probes against the Pfam domains they overlap",
    subtitle = sprintf(
      paste("The %d domains overlapped most often, the rest pooled, log colour",
            "scale. A probe meeting its own domains is the validation agreement."),
      .TOP_DOMAINS
    ),
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}


# Retrotransposon length against the number of unreduced loci overlapping it
# (from the ltr_structure table). Longer, hit-rich elements are the
# best-supported candidate proviruses.
retro_length_vs_hits_plot <- function(ltr_df, subset_label = NULL,
                                      warning_caption = NULL) {
  if (nrow(ltr_df) == 0L) return(empty_plot())
  d <- ltr_df %>% dplyr::filter(n_overlapping_hits > 0L)
  if (nrow(d) == 0L) return(empty_plot("no hit-overlapping retrotransposons"))

  p <- ggplot(d, aes(x = width, y = n_overlapping_hits)) +
    geom_point(alpha = 0.35, size = 1, colour = .DATA_COLOUR) +
    scale_x_continuous(labels = scales::label_comma()) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(x = "Retrotransposon length (bp)", y = "Unreduced loci overlapping it")
  add_titles(
    p,
    title    = "Element length against the loci it holds",
    subtitle = "Each LTRdigest element with at least one overlapping tBLASTn locus.",
    subset_label    = subset_label,
    warning_caption = warning_caption
  )
}
