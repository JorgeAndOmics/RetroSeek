# -----------------------------------------------------------------------------
# range_analysis / stage_dataframe.R
# -----------------------------------------------------------------------------
# Build the per-genome "middle stage" dataframes consumed by
# stage_plot_generator.R. Where plot_dataframe.R characterises the final
# (valid-reduced) tier, these tables characterise the homology + LTRharvest/
# LTRdigest integration stage: how tBLASTn hits sit relative to LTR elements,
# how structurally complete the LTR calls are, and how many primitive
# sequences collapse into each reduced locus (the multiplicity metrics).
#
# Three single-grain tables (one parquet each), mirroring the one-parquet-per-
# genome pattern of plot_dataframe.R:
#   * hits     — one row per gr_virus locus (homology tier): spatial
#                concordance vs retrotransposons, tier flags, M1 (n_hits).
#   * ltr      — one row per retrotransposon: structural completeness +
#                probe-domain composition.
#   * reduced  — one row per gr_global locus: M2 (n_loci) + M1 roll-up.

suppressMessages({
  library(GenomicRanges)
  library(S4Vectors)
})


# Window (bp) within which a non-overlapping hit counts as "flanking" rather
# than "disjoint" relative to the nearest retrotransposon.
.STAGE_FLANKING_WINDOW <- 1000L


# One row per gr_virus locus. `concordance` classifies each homology locus by
# its spatial relationship to the LTRdigest retrotransposons; `is_candidate` /
# `is_valid` flag whether the locus survives each refinement step. `n_hits` is
# M1 — the count of threshold-passing raw tBLASTn hits collapsed into the
# locus by reduce_first.
build_stage_hits_df <- function(gr_virus, retrotransposons, candidate_hits,
                                valid_hits,
                                flanking_window = .STAGE_FLANKING_WINDOW) {
  empty <- tibble::tibble(
    seqnames = character(0), start = integer(0), end = integer(0),
    width = integer(0), strand = character(0),
    probe = character(0), virus = character(0), label = character(0),
    species = character(0), n_hits = integer(0),
    max_bitscore = numeric(0), max_identity = numeric(0),
    query_coverage = numeric(0), concordance = character(0),
    is_candidate = logical(0), is_valid = logical(0)
  )
  if (length(gr_virus) == 0L) return(empty)

  # Spatial concordance: inside (overlaps a retrotransposon) / flanking (within
  # `flanking_window` bp of one) / disjoint (farther, or no LTR calls at all).
  concordance <- rep("disjoint", length(gr_virus))
  if (length(retrotransposons) > 0L) {
    d2n <- GenomicRanges::distanceToNearest(gr_virus, retrotransposons,
                                            ignore.strand = TRUE)
    if (length(d2n) > 0L) {
      qh   <- S4Vectors::queryHits(d2n)
      dist <- S4Vectors::mcols(d2n)$distance
      concordance[qh[dist == 0L]] <- "inside"
      concordance[qh[dist > 0L & dist <= flanking_window]] <- "flanking"
    }
  }

  ids <- as.character(S4Vectors::mcols(gr_virus)$ID)
  cand_ids  <- as.character(S4Vectors::mcols(candidate_hits)$ID)
  valid_ids <- as.character(S4Vectors::mcols(valid_hits)$ID)

  df <- as.data.frame(gr_virus, stringsAsFactors = FALSE)
  tibble::tibble(
    seqnames       = as.character(df$seqnames),
    start          = as.integer(df$start),
    end            = as.integer(df$end),
    width          = as.integer(df$width),
    strand         = as.character(df$strand),
    probe          = as.character(df$probe),
    virus          = as.character(df$virus),
    label          = as.character(df$label),
    species        = as.character(df$species),
    n_hits         = as.integer(df$n_hits),          # M1
    max_bitscore   = as.numeric(df$max_bitscore),
    max_identity   = as.numeric(df$max_identity),
    query_coverage = as.numeric(df$query_coverage),
    concordance    = concordance,
    is_candidate   = ids %in% cand_ids,
    is_valid       = ids %in% valid_ids
  )
}


# One row per retrotransposon. Structural-completeness flags are derived by
# grouping the LTRdigest child features by `Parent` ID. Two `Parent` namespaces
# are in play: flanking LTRs, Pfam domains and RR-tracts (PPT) are children of
# the `LTR_retrotransposon` (so they join on `retro_ids`), while target-site
# duplications are children of the enclosing `repeat_region` (so they join on
# each retrotransposon's own `Parent`, i.e. `retro_parent`). `domain_probes` is
# a "; "-joined sorted set for the domain-composition plot.
build_stage_ltr_df <- function(retrotransposons, flanking_ltrs,
                               domains_w_probes, ltr_data, gr_virus) {
  empty <- tibble::tibble(
    seqnames = character(0), start = integer(0), end = integer(0),
    width = integer(0), strand = character(0), ID = character(0),
    n_flanking_ltrs = integer(0), has_both_ltrs = logical(0),
    n_probe_domains = integer(0), n_domains_total = integer(0),
    domain_probes = character(0),
    has_tsd = logical(0), n_tsd = integer(0),
    has_ppt = logical(0), n_ppt = integer(0),
    n_overlapping_hits = integer(0)
  )
  if (length(retrotransposons) == 0L) return(empty)

  retro_ids <- as.character(S4Vectors::mcols(retrotransposons)$ID)
  retro_lvl <- factor(retro_ids, levels = retro_ids)
  # Each retrotransposon's own `Parent` is the enclosing `repeat_region` ID —
  # the namespace that target-site-duplication features are parented to.
  retro_parent <- as.character(S4Vectors::mcols(retrotransposons)$Parent)

  # Flanking LTR arms per parent retrotransposon.
  flank_parent <- if (length(flanking_ltrs) > 0L) {
    as.character(S4Vectors::mcols(flanking_ltrs)$Parent)
  } else character(0)
  n_flank <- as.integer(table(factor(flank_parent, levels = retro_ids)))

  # Probe-assigned Pfam domains per parent (the config-regex-matched subset).
  dom_parent <- if (length(domains_w_probes) > 0L) {
    as.character(S4Vectors::mcols(domains_w_probes)$Parent)
  } else character(0)
  dom_probe <- if (length(domains_w_probes) > 0L) {
    as.character(S4Vectors::mcols(domains_w_probes)$probe)
  } else character(0)
  n_probe_domains <- as.integer(table(factor(dom_parent, levels = retro_ids)))
  probe_by_parent <- split(dom_probe, factor(dom_parent, levels = retro_ids))
  domain_probes <- vapply(probe_by_parent, function(p) {
    if (length(p) == 0L) return(NA_character_)
    paste(sort(unique(p)), collapse = "; ")
  }, character(1))

  # All Pfam `protein_match` features per parent, regardless of probe
  # assignment — LTRdigest's view of coding capacity, independent of the probe
  # panel. `protein_match` is a child of the `LTR_retrotransposon`.
  pm <- ltr_data[ltr_data$type == "protein_match"]
  pm_parent <- if (length(pm) > 0L) {
    as.character(S4Vectors::mcols(pm)$Parent)
  } else character(0)
  n_domains_total <- as.integer(table(factor(pm_parent, levels = retro_ids)))

  # Target-site duplications per parent. TSDs are children of the enclosing
  # `repeat_region`, not the `LTR_retrotransposon` — so they join on
  # `retro_parent` (the retrotransposon's own `Parent`), not `retro_ids`.
  tsd <- ltr_data[ltr_data$type == "target_site_duplication"]
  tsd_parent <- if (length(tsd) > 0L) {
    as.character(S4Vectors::mcols(tsd)$Parent)
  } else character(0)
  n_tsd <- as.integer(table(factor(tsd_parent, levels = retro_parent)))
  has_tsd <- n_tsd > 0L

  # Polypurine tracts (RR_tract / PPT) per parent — a marker of structural
  # intactness. Children of the `LTR_retrotransposon`, so they join on
  # `retro_ids`. LTRdigest emits these for only a fraction of elements.
  ppt <- ltr_data[ltr_data$type == "RR_tract"]
  ppt_parent <- if (length(ppt) > 0L) {
    as.character(S4Vectors::mcols(ppt)$Parent)
  } else character(0)
  n_ppt <- as.integer(table(factor(ppt_parent, levels = retro_ids)))
  has_ppt <- n_ppt > 0L

  # tBLASTn (gr_virus) loci landing inside each retrotransposon.
  n_hits_in <- if (length(gr_virus) > 0L) {
    as.integer(GenomicRanges::countOverlaps(retrotransposons, gr_virus,
                                            ignore.strand = TRUE))
  } else {
    rep(0L, length(retrotransposons))
  }

  df <- as.data.frame(retrotransposons, stringsAsFactors = FALSE)
  tibble::tibble(
    seqnames           = as.character(df$seqnames),
    start              = as.integer(df$start),
    end                = as.integer(df$end),
    width              = as.integer(df$width),
    strand             = as.character(df$strand),
    ID                 = retro_ids,
    n_flanking_ltrs    = n_flank,
    has_both_ltrs      = n_flank >= 2L,
    n_probe_domains    = n_probe_domains,
    n_domains_total    = n_domains_total,
    domain_probes      = unname(domain_probes),
    has_tsd            = has_tsd,
    n_tsd              = n_tsd,
    has_ppt            = has_ppt,
    n_ppt              = n_ppt,
    n_overlapping_hits = n_hits_in
  )
}


# One row per gr_global locus. `n_loci` is M2 — the count of gr_virus loci
# collapsed into this globally-reduced locus by reduce_global; `n_hits` is the
# M1 roll-up (total raw hits beneath it).
build_stage_reduced_df <- function(gr_global) {
  empty <- tibble::tibble(
    seqnames = character(0), start = integer(0), end = integer(0),
    width = integer(0), strand = character(0), ID = character(0),
    probe = character(0), virus = character(0), label = character(0),
    n_hits = integer(0), n_loci = integer(0), max_bitscore = numeric(0)
  )
  if (length(gr_global) == 0L) return(empty)

  df <- as.data.frame(gr_global, stringsAsFactors = FALSE)
  tibble::tibble(
    seqnames     = as.character(df$seqnames),
    start        = as.integer(df$start),
    end          = as.integer(df$end),
    width        = as.integer(df$width),
    strand       = as.character(df$strand),
    ID           = as.character(df$ID),
    probe        = as.character(df$probe),
    virus        = as.character(df$virus),
    label        = as.character(df$label),
    n_hits       = as.integer(df$n_hits),    # M1 roll-up
    n_loci       = as.integer(df$n_loci),    # M2
    max_bitscore = as.numeric(df$max_bitscore)
  )
}
