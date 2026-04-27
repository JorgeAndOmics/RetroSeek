# -----------------------------------------------------------------------------
# range_analysis / validation.R
# -----------------------------------------------------------------------------
# Two refinement steps applied to the reduced BLAST GRanges:
#
#   - candidate hits  = reduced BLAST hits that overlap an LTR retrotransposon
#                       (ERV); these are the spatial candidates for ERV identity.
#   - valid hits      = candidates whose probe label matches a probe label of
#                       at least one Pfam-annotated domain inside the *same*
#                       enclosing retrotransposon. This is the semantic check
#                       that turns "near-an-LTR" into "consistent with the
#                       enclosing ERV's molecular architecture".
#
# This replaces the earlier triple-join cascade (gr_virus ⨝ ltr_valid_hits ⨝
# ltr_domain with brittle suffix juggling) with a clean two-step approach
# using findOverlaps + per-retrotransposon probe-set membership.

suppressMessages({
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})


# Subset `gr_hits` to those overlapping any retrotransposon.
find_candidate_hits <- function(gr_hits, retrotransposons) {
  if (length(gr_hits) == 0L || length(retrotransposons) == 0L) return(gr_hits[FALSE])
  # Strand-aware overlap: BLAST hits carry +/- strand; LTR retrotransposons
  # carry +/-/* (the * cases are LTRs whose strand LTRharvest could not infer).
  # GenomicRanges' findOverlaps with ignore.strand=FALSE treats * as wildcard,
  # so * retros still match either-strand BLAST hits while +/- mismatches
  # (a + BLAST hit vs a - retro) correctly do NOT overlap.
  ov <- GenomicRanges::findOverlaps(gr_hits, retrotransposons, ignore.strand = FALSE)
  gr_hits[unique(S4Vectors::queryHits(ov))]
}


# Build a list mapping each retrotransposon (by its ID attribute) to the set
# of probes assigned to its child Pfam domains. Domains carry a `Parent`
# attribute pointing to their enclosing retrotransposon's ID; we group by
# parent and union the probe assignments.
build_retrotransposon_probe_sets <- function(retrotransposons, domains_with_probes) {
  if (length(retrotransposons) == 0L) return(list())

  retro_ids <- as.character(retrotransposons$ID)
  # Initialise empty probe set for every retrotransposon — even ones with no
  # Pfam-domain children, so a candidate over a "domain-empty" ERV cleanly
  # gets the empty set and is dropped by the membership test.
  probe_sets <- setNames(rep(list(character(0)), length(retro_ids)), retro_ids)

  if (length(domains_with_probes) == 0L) return(probe_sets)

  parents <- as.character(domains_with_probes$Parent)
  probes  <- as.character(domains_with_probes$probe)
  for (i in seq_along(parents)) {
    pid <- parents[i]
    if (!is.null(probe_sets[[pid]])) {
      probe_sets[[pid]] <- unique(c(probe_sets[[pid]], probes[i]))
    }
  }
  probe_sets
}


# A candidate hit is "valid" if its probe is in the probe-set of at least one
# retrotransposon it overlaps. For ranges-with-multi-value probes (concatenate
# strategy emits "POL; GAG"), we split on the configured separator and pass
# if any of the candidate's probes is in the enclosing retrotransposon's set.
find_valid_hits <- function(gr_candidates, retrotransposons, domains_with_probes,
                            concat_separator = "; ") {
  if (length(gr_candidates) == 0L || length(retrotransposons) == 0L) {
    return(gr_candidates[FALSE])
  }

  probe_sets <- build_retrotransposon_probe_sets(retrotransposons, domains_with_probes)

  # For each candidate, find ALL overlapping retrotransposons and union their
  # probe-sets — a candidate may straddle two close retros, and matching either
  # one passes.
  ov <- GenomicRanges::findOverlaps(gr_candidates, retrotransposons, ignore.strand = FALSE)
  if (length(ov) == 0L) return(gr_candidates[FALSE])

  qhits <- S4Vectors::queryHits(ov)
  shits <- S4Vectors::subjectHits(ov)
  retro_ids_per_subj <- as.character(retrotransposons$ID)[shits]

  # Group: for each candidate index q, the union of probe-sets across all its
  # overlapping retrotransposons.
  per_candidate_set <- split(retro_ids_per_subj, qhits)
  per_candidate_set <- lapply(per_candidate_set, function(ids) {
    unlist(probe_sets[ids], use.names = FALSE)
  })

  # Probe of each candidate; concatenate strategy can pack multiple probes
  # into one string, so split before the membership check.
  probes_chr <- as.character(S4Vectors::mcols(gr_candidates)$probe)
  probes_split <- strsplit(probes_chr, concat_separator, fixed = TRUE)

  # Decide pass/fail for the candidates that actually overlap something.
  candidate_idx <- as.integer(names(per_candidate_set))
  valid_mask <- rep(FALSE, length(gr_candidates))
  for (k in seq_along(candidate_idx)) {
    q <- candidate_idx[k]
    valid_mask[q] <- length(intersect(probes_split[[q]], per_candidate_set[[k]])) > 0L
  }
  gr_candidates[valid_mask]
}
