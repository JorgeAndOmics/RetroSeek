# -----------------------------------------------------------------------------
# ranges / validation.R
# -----------------------------------------------------------------------------
# Refinement steps applied to the reduced BLAST GRanges:
#
#   - candidate hits  = reduced BLAST hits that overlap an LTR retrotransposon
#                       (ERV); these are the anchored candidates for ERV identity.
#   - anchored hits   = the candidate set, LABELLED (not filtered): every
#                       candidate is kept and annotated with its enclosing
#                       element (`Parent`), a per-provirus `domain_tier`, and a
#                       per-hit `domain_hit_class`. Nothing is discarded — the
#                       "valid" tier is now the whole anchored set carrying the
#                       labels needed to judge domain support downstream (ADR-009).
#   - orphan hits     = the strand-aware complement: hits overlapping no element.
#
# This replaces the earlier find_valid_hits filter (which dropped every anchored
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
# overlap NO retrotransposon. These are the non-LTR-associated orphans — solo
# ORFs, degraded proviruses, and candidate novel retroviruses whose LTRs are too
# diverged for LTRharvest to pair. They are recovered into the orphan tier and
# classified by their own sequence (taxonomy_classify_loci.py --source orphan).
# With no retrotransposons, every hit is unanchored.
find_unanchored_hits <- function(gr_hits, retrotransposons) {
  if (length(gr_hits) == 0L) return(gr_hits[FALSE])
  if (length(retrotransposons) == 0L) return(gr_hits)
  ov <- GenomicRanges::findOverlaps(gr_hits, retrotransposons, ignore.strand = FALSE)
  anchored <- unique(S4Vectors::queryHits(ov))
  gr_hits[setdiff(seq_along(gr_hits), anchored)]
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


# Per-retrotransposon boolean: does the element carry at least one LTRdigest
# protein domain child (any Pfam, regardless of the config probe filter)?
# `all_domains` are the `protein_match` features (extract_all_domains); each
# points at its enclosing LTR_retrotransposon via `Parent`. This separates the
# `domain_unlisted` tier (has domains, none config-matched) from `non_domain`.
build_retrotransposon_domain_presence <- function(retrotransposons, all_domains) {
  retro_ids <- as.character(retrotransposons$ID)
  presence  <- setNames(rep(FALSE, length(retro_ids)), retro_ids)
  if (length(retro_ids) == 0L || length(all_domains) == 0L) return(presence)
  parents <- unique(as.character(all_domains$Parent))
  presence[intersect(parents, retro_ids)] <- TRUE
  presence
}


# Domain-tier precedence (strongest wins when a hit straddles several elements).
.DOMAIN_TIER_RANK <- c(non_domain = 0L, domain_unlisted = 1L, domain_selected = 2L)


# Per-hit positional domain class (hit_domain_mode == "positional"): a hit is
# `substring_match` when it physically overlaps a config-matched domain OF ITS
# OWN gene, `no_substring_match` when it overlaps some other protein domain, and
# `non_domain` when it overlaps none. `probes_split` is the per-candidate probe
# vector (concatenate strategy may pack several probes into one hit).
.positional_hit_class <- function(gr_candidates, domains_with_probes, all_domains,
                                  probes_split) {
  n   <- length(gr_candidates)
  cls <- rep("non_domain", n)
  if (length(all_domains) > 0L) {
    any_ov <- unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(
      gr_candidates, all_domains, ignore.strand = TRUE)))
    cls[any_ov] <- "no_substring_match"
  }
  if (length(domains_with_probes) > 0L) {
    ov <- GenomicRanges::findOverlaps(gr_candidates, domains_with_probes,
                                      ignore.strand = TRUE)
    qh <- S4Vectors::queryHits(ov)
    dom_probe <- as.character(S4Vectors::mcols(domains_with_probes)$probe)[
      S4Vectors::subjectHits(ov)]
    for (i in seq_along(qh)) {
      q <- qh[i]
      if (dom_probe[i] %in% probes_split[[q]]) cls[q] <- "substring_match"
    }
  }
  cls
}


# Annotate every candidate (LTR-anchored) hit WITHOUT discarding any. Emits three
# mcols columns (ADR-009):
#
#   Parent            greatest-overlap LTR_retrotransposon id — the anchor the
#                     taxonomic classifier groups a locus's per-gene hits by.
#   domain_tier       per-provirus, strongest across straddled elements:
#                       domain_selected   >=1 config-matched domain (any gene)
#                       domain_unlisted   has protein domains, none config-matched
#                       non_domain        no protein domain at all
#   domain_hit_class  per-hit, controlled by `hit_domain_mode`:
#                       membership  substring_match iff the hit's own gene has a
#                                   config-matched domain in an enclosing element
#                                   (co-occurrence; the old valid gate), else
#                                   no_substring_match
#                       positional  see .positional_hit_class (co-localization;
#                                   adds a non_domain level)
annotate_anchored_hits <- function(gr_candidates, retrotransposons,
                                   domains_with_probes, all_domains,
                                   hit_domain_mode = "membership",
                                   concat_separator = "; ") {
  if (length(gr_candidates) == 0L || length(retrotransposons) == 0L) {
    S4Vectors::mcols(gr_candidates)$Parent           <- character(0)
    S4Vectors::mcols(gr_candidates)$domain_tier      <- character(0)
    S4Vectors::mcols(gr_candidates)$domain_hit_class <- character(0)
    return(gr_candidates)
  }

  probe_sets <- build_retrotransposon_probe_sets(retrotransposons, domains_with_probes)
  has_domain <- build_retrotransposon_domain_presence(retrotransposons, all_domains)

  # Element-wise tier for every retrotransposon (config-matched > any-domain > none).
  retro_ids <- as.character(retrotransposons$ID)
  elem_tier <- vapply(retro_ids, function(id) {
    if (length(probe_sets[[id]]) > 0L) "domain_selected"
    else if (isTRUE(has_domain[[id]])) "domain_unlisted"
    else "non_domain"
  }, character(1))

  ov    <- GenomicRanges::findOverlaps(gr_candidates, retrotransposons, ignore.strand = FALSE)
  qhits <- S4Vectors::queryHits(ov)
  shits <- S4Vectors::subjectHits(ov)
  retro_ids_per_subj <- retro_ids[shits]

  # Concatenate strategy can pack multiple probes into one hit ("POL; GAG").
  probes_chr   <- as.character(S4Vectors::mcols(gr_candidates)$probe)
  probes_split <- strsplit(probes_chr, concat_separator, fixed = TRUE)

  n           <- length(gr_candidates)
  domain_tier <- rep("non_domain", n)   # every candidate overlaps >=1 element
  parent_of   <- rep(NA_character_, n)
  best_w      <- rep(-1L, n)
  membership_set <- vector("list", n)   # union of config probe-sets per candidate

  ov_widths <- IRanges::width(GenomicRanges::pintersect(
    gr_candidates[qhits], retrotransposons[shits], ignore.strand = TRUE))
  for (i in seq_along(qhits)) {
    q  <- qhits[i]
    id <- retro_ids_per_subj[i]
    t  <- elem_tier[[id]]
    if (.DOMAIN_TIER_RANK[[t]] > .DOMAIN_TIER_RANK[[domain_tier[q]]]) domain_tier[q] <- t
    if (ov_widths[i] > best_w[q]) { best_w[q] <- ov_widths[i]; parent_of[q] <- id }
    membership_set[[q]] <- c(membership_set[[q]], probe_sets[[id]])
  }

  if (identical(hit_domain_mode, "positional")) {
    domain_hit_class <- .positional_hit_class(gr_candidates, domains_with_probes,
                                              all_domains, probes_split)
  } else {
    domain_hit_class <- vapply(seq_len(n), function(q) {
      if (length(intersect(probes_split[[q]], membership_set[[q]])) > 0L)
        "substring_match" else "no_substring_match"
    }, character(1))
  }

  S4Vectors::mcols(gr_candidates)$Parent           <- parent_of
  S4Vectors::mcols(gr_candidates)$domain_tier      <- domain_tier
  S4Vectors::mcols(gr_candidates)$domain_hit_class <- domain_hit_class
  gr_candidates
}
