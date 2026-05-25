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


# One row per gr_virus locus, characterising PRE-REDUCTION spatial redundancy —
# how much the unreduced original tier overlaps itself (the redundancy that
# reduction collapses). `overlap_degree` is the number of OTHER gr_virus loci
# this one overlaps (strand-agnostic — redundancy is positional);
# `max_reciprocal_fraction` is the largest intersection / min(width) over those
# partners (1.0 = fully contained in / containing a neighbour).
build_stage_overlap_df <- function(gr_virus) {
  empty <- tibble::tibble(
    seqnames = character(0), start = integer(0), end = integer(0),
    width = integer(0), strand = character(0), probe = character(0),
    virus = character(0), label = character(0),
    overlap_degree = integer(0), max_reciprocal_fraction = numeric(0)
  )
  if (length(gr_virus) == 0L) return(empty)

  hits <- GenomicRanges::findOverlaps(gr_virus, gr_virus, ignore.strand = TRUE)
  q <- S4Vectors::queryHits(hits); s <- S4Vectors::subjectHits(hits)
  keep <- q != s
  q <- q[keep]; s <- s[keep]
  overlap_degree <- tabulate(q, nbins = length(gr_virus))
  max_rec <- numeric(length(gr_virus))
  if (length(q) > 0L) {
    inter <- BiocGenerics::width(
      GenomicRanges::pintersect(gr_virus[q], gr_virus[s], ignore.strand = TRUE)
    )
    frac <- inter / pmin(BiocGenerics::width(gr_virus)[q],
                         BiocGenerics::width(gr_virus)[s])
    per_q <- tapply(frac, q, max)
    max_rec[as.integer(names(per_q))] <- as.numeric(per_q)
  }

  df <- as.data.frame(gr_virus, stringsAsFactors = FALSE)
  tibble::tibble(
    seqnames = as.character(df$seqnames), start = as.integer(df$start),
    end = as.integer(df$end), width = as.integer(df$width),
    strand = as.character(df$strand), probe = as.character(df$probe),
    virus = as.character(df$virus), label = as.character(df$label),
    overlap_degree = as.integer(overlap_degree),
    max_reciprocal_fraction = as.numeric(max_rec)
  )
}


# One row per gr_virus locus, characterising its interaction with the nearest
# LTRdigest retrotransposon. `feature_class` refines the homology concordance
# with a `domain_overlap` level (overlaps a probe-assigned Pfam domain — the
# validation signal); `relative_position_in_retro` is the strand-aware 5'→3'
# position [0,1] of the locus midpoint within its enclosing retrotransposon
# (NA when not inside); `strand_concordant` compares the hit strand to the
# enclosing element's (NA when not inside; `*` retro strand counts as concordant).
build_stage_ltr_interaction_df <- function(gr_virus, retrotransposons,
                                           domains_with_probes,
                                           flanking_window = .STAGE_FLANKING_WINDOW) {
  empty <- tibble::tibble(
    seqnames = character(0), start = integer(0), end = integer(0),
    strand = character(0), probe = character(0), virus = character(0),
    distance_to_nearest_retro = numeric(0), feature_class = character(0),
    relative_position_in_retro = numeric(0), strand_concordant = logical(0),
    enclosing_retro_width = integer(0)
  )
  if (length(gr_virus) == 0L) return(empty)

  n <- length(gr_virus)
  dist  <- rep(NA_real_, n)
  nearest_subj <- rep(NA_integer_, n)
  if (length(retrotransposons) > 0L) {
    d2n <- GenomicRanges::distanceToNearest(gr_virus, retrotransposons,
                                            ignore.strand = TRUE)
    if (length(d2n) > 0L) {
      qh <- S4Vectors::queryHits(d2n)
      dist[qh]         <- as.numeric(S4Vectors::mcols(d2n)$distance)
      nearest_subj[qh] <- S4Vectors::subjectHits(d2n)
    }
  }
  inside <- !is.na(dist) & dist == 0

  # Relative position + strand concordance within the enclosing (nearest, dist 0)
  # retrotransposon.
  rel_pos <- rep(NA_real_, n)
  strand_conc <- rep(NA, n)
  retro_w <- rep(NA_integer_, n)
  if (any(inside)) {
    idx   <- which(inside)
    rsub  <- nearest_subj[idx]
    rstart <- BiocGenerics::start(retrotransposons)[rsub]
    rwidth <- BiocGenerics::width(retrotransposons)[rsub]
    rstr   <- as.character(BiocGenerics::strand(retrotransposons))[rsub]
    mid    <- (BiocGenerics::start(gr_virus)[idx] + BiocGenerics::end(gr_virus)[idx]) / 2
    frac   <- pmin(1, pmax(0, (mid - rstart) / rwidth))
    rel_pos[idx] <- ifelse(rstr == "-", 1 - frac, frac)
    hstr <- as.character(BiocGenerics::strand(gr_virus))[idx]
    strand_conc[idx] <- (rstr == "*") | (hstr == rstr)
    retro_w[idx] <- as.integer(rwidth)
  }

  # feature_class: domain_overlap > inside_retro > flanking_ltr > disjoint.
  feature_class <- rep("disjoint", n)
  feature_class[!is.na(dist) & dist > 0 & dist <= flanking_window] <- "flanking_ltr"
  feature_class[inside] <- "inside_retro"
  if (length(domains_with_probes) > 0L) {
    dom_ov <- S4Vectors::queryHits(
      GenomicRanges::findOverlaps(gr_virus, domains_with_probes, ignore.strand = TRUE)
    )
    feature_class[unique(dom_ov)] <- "domain_overlap"
  }

  tibble::tibble(
    seqnames = as.character(GenomicRanges::seqnames(gr_virus)),
    start    = as.integer(BiocGenerics::start(gr_virus)),
    end      = as.integer(BiocGenerics::end(gr_virus)),
    strand   = as.character(BiocGenerics::strand(gr_virus)),
    probe    = as.character(S4Vectors::mcols(gr_virus)$probe),
    virus    = as.character(S4Vectors::mcols(gr_virus)$virus),
    distance_to_nearest_retro  = dist,
    feature_class              = feature_class,
    relative_position_in_retro = rel_pos,
    strand_concordant          = strand_conc,
    enclosing_retro_width      = retro_w
  )
}


# Long table (one row per gr_virus-locus × overlapped probe-domain pair): the
# raw co-occurrence used by the probe × Pfam-domain heatmap. `hit_probe` is the
# tBLASTn locus's probe; `domain_probe` is the probe assigned to the overlapped
# Pfam domain. Their agreement is exactly the validation signal made visible.
build_stage_probe_domain_df <- function(gr_virus, domains_with_probes) {
  empty <- tibble::tibble(hit_probe = character(0), domain_probe = character(0))
  if (length(gr_virus) == 0L || length(domains_with_probes) == 0L) return(empty)
  ov <- GenomicRanges::findOverlaps(gr_virus, domains_with_probes,
                                    ignore.strand = TRUE)
  if (length(ov) == 0L) return(empty)
  tibble::tibble(
    hit_probe    = as.character(S4Vectors::mcols(gr_virus)$probe)[S4Vectors::queryHits(ov)],
    domain_probe = as.character(S4Vectors::mcols(domains_with_probes)$probe)[S4Vectors::subjectHits(ov)]
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
