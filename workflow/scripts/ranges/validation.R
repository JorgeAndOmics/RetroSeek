# -----------------------------------------------------------------------------
# ranges / validation.R
# -----------------------------------------------------------------------------
# Refinement steps applied to the reduced BLAST GRanges:
#
#   - candidate hits  = reduced BLAST hits that overlap an LTR retrotransposon
#                       (ERV); these are the LTR-flanked candidates for ERV identity.
#   - LTR-flanked hits   = the candidate set, LABELLED (not filtered): every
#                       candidate is kept and annotated with its enclosing
#                       element (`Parent`). Nothing is discarded - the
#                       "valid" tier is now the whole LTR-flanked set carrying the
#                       labels needed to judge domain support downstream (ADR-009).
#   - orphan hits     = the strand-aware complement: hits overlapping no element.
#
# This replaces the earlier find_element_hits filter (which dropped every LTR-flanked
# hit lacking a matching Pfam domain) with a findOverlaps + per-retrotransposon
# probe-set membership annotation.

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


# The strand-aware complement of find_candidate_hits: reduced BLAST hits that
# overlap NO retrotransposon. These are the non-LTR-associated orphans - solo
# ORFs, degraded proviruses, and candidate novel retroviruses whose LTRs are too
# diverged for LTRharvest to pair. They are recovered into the orphan tier and
# classified by their own sequence (taxonomy_classify_loci.py --source orphan).
# With no retrotransposons, every hit is an orphan.
find_orphan_hits <- function(gr_hits, retrotransposons) {
  if (length(gr_hits) == 0L) return(gr_hits[FALSE])
  if (length(retrotransposons) == 0L) return(gr_hits)
  ov <- GenomicRanges::findOverlaps(gr_hits, retrotransposons, ignore.strand = FALSE)
  ltr_flanked <- unique(S4Vectors::queryHits(ov))
  gr_hits[setdiff(seq_along(gr_hits), ltr_flanked)]
}


# Cluster orphan (non-LTR-associated) hits by physical OVERLAP only and stamp each
# with a synthetic `Parent`, so the taxonomic classifier's build_loci groups them
# into single, non-overlapping orphan loci - the grouping an LTR-flanked provirus gets
# from its LTR element, but keyed on co-location because an orphan has no element
# (ADR-010). Only hits whose ranges OVERLAP (or are book-ended) merge; gap-separated
# hits stay separate. Strand is ignored (a degraded provirus's genes may be on
# either strand); the locus strand is later taken as the member-strand mode.
#
# Overlap is EVIDENCE the hits are the same feature (vs. the retired proximity
# window, which INFERRED a provirus from nearness). Consequence: adjacent genes
# (gag/pol/env occupy distinct, non-overlapping positions) do NOT merge, so orphan
# loci are mostly single-gene - a deliberate, conservative deduplication rather
# than speculative multi-gene assembly. Orphan loci stay flagged `source=orphan`.
#
# `max_provirus_len` is the per-genome ground-truth cap: the widest LTRdigest
# LTR_retrotransposon. A cluster wider than any real provirus can't be one, so it
# is FLAGGED (`oversized`), not dropped, for downstream filtering. Inf = no cap.
cluster_orphan_hits <- function(orphan_hits, max_provirus_len = Inf) {
  if (length(orphan_hits) == 0L) {
    S4Vectors::mcols(orphan_hits)$Parent    <- character(0)
    S4Vectors::mcols(orphan_hits)$oversized <- character(0)
    return(orphan_hits)
  }
  # Overlap-only: reduce merges ranges separated by a gap < min.gapwidth, so
  # min.gapwidth = 1 keeps only overlapping / book-ended ranges together.
  clusters <- GenomicRanges::reduce(orphan_hits, min.gapwidth = 1L,
                                    ignore.strand = TRUE)
  cl_oversized <- BiocGenerics::width(clusters) > max_provirus_len
  ov <- GenomicRanges::findOverlaps(orphan_hits, clusters, ignore.strand = TRUE)
  cl <- rep(NA_integer_, length(orphan_hits))
  cl[S4Vectors::queryHits(ov)] <- S4Vectors::subjectHits(ov)
  S4Vectors::mcols(orphan_hits)$Parent <- sprintf(
    "orphan_%s_%d",
    as.character(GenomicRanges::seqnames(clusters))[cl],
    BiocGenerics::start(clusters)[cl]
  )
  # String bool to match the classifier's convention (parse_valid_full reads it).
  S4Vectors::mcols(orphan_hits)$oversized <- ifelse(cl_oversized[cl], "True", "False")
  orphan_hits
}


# Annotate every candidate (LTR-flanked) hit WITHOUT discarding any. Emits one
# mcols column:
#
#   Parent   greatest-overlap LTR_retrotransposon id - the anchor the taxonomic
#            classifier groups a locus's per-gene hits by.
#
# `domain_tier` no longer lives here (ADR-016). It was computed from LTRdigest's
# domain NAMES at ELEMENT grain, while the catalog's `domain_tier` is computed
# from the domain scan's ACCESSIONS at LOCUS grain. Two columns with one name and
# two answers is a trap, so the curated classification happens once, in the scan.
# What LTRdigest still contributes is `n_domains_total`, a count (stage_dataframe.R).
annotate_ltr_flanked_hits <- function(gr_candidates, retrotransposons) {
  if (length(gr_candidates) == 0L || length(retrotransposons) == 0L) {
    S4Vectors::mcols(gr_candidates)$Parent <- character(0)
    return(gr_candidates)
  }

  retro_ids <- as.character(retrotransposons$ID)
  ov    <- GenomicRanges::findOverlaps(gr_candidates, retrotransposons, ignore.strand = FALSE)
  qhits <- S4Vectors::queryHits(ov)
  shits <- S4Vectors::subjectHits(ov)
  retro_ids_per_subj <- retro_ids[shits]

  n         <- length(gr_candidates)
  parent_of <- rep(NA_character_, n)
  best_w    <- rep(-1L, n)

  ov_widths <- IRanges::width(GenomicRanges::pintersect(
    gr_candidates[qhits], retrotransposons[shits], ignore.strand = TRUE))
  for (i in seq_along(qhits)) {
    q <- qhits[i]
    if (ov_widths[i] > best_w[q]) {
      best_w[q] <- ov_widths[i]
      parent_of[q] <- retro_ids_per_subj[i]
    }
  }

  S4Vectors::mcols(gr_candidates)$Parent <- parent_of
  gr_candidates
}
