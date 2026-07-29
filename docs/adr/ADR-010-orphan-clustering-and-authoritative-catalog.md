# ADR-010: Orphan proximity-clustering, `fragment`->`orphan` rename, and the authoritative catalog

- **Status**: Accepted
- **Date**: 2026-07-10
- **Deciders**: Jorge González García
- **Refines**: [ADR-009](ADR-009-anchored-domain-tiering-and-structure-class.md) (anchored tiering + orphan tier)

## Context

After ADR-009, anchored proviruses were a clean per-locus catalog - one row per
LTR element (grouped by `Parent`), non-overlapping, each with an authoritative
`taxon_call` + `structure_class` + `is_mosaic`. **Orphans were not.** An orphan
(a non-LTR-associated hit) has no LTR element, so `build_loci` made **each
parentless gene-hit its own singleton locus** - meaning orphans were always
single-gene (`structure_class = gene`, never mosaic) and could **spatially
overlap** each other (a POL orphan and a GAG orphan on the same degraded provirus
became two overlapping rows: ~1,478 overlaps for *Antrozous*).

Two consequences: (1) there was no single, non-overlapping, authoritative locus
table across both tiers; (2) the vocabulary was split-brain - the concept was
"orphan" (the `source` value, loss metrics, plot labels since ADR-009) but the
plumbing still said "fragments" (dirs, files, rule names, CLI args, functions).

## Decision

1. **Orphan overlap-clustering, capped (synthetic `Parent`).** `cluster_orphan_hits`
   (ranges/validation.R) merges **only orphan hits whose ranges physically
   overlap** (`GenomicRanges::reduce`, `min.gapwidth = 1`, strand-ignored) and
   stamps each cluster with a synthetic `Parent` (`orphan_<seqname>_<clusterStart>`).
   The classifier's existing `build_loci`/`_assemble` then group orphans **verbatim
   the anchored path**. Overlap is *evidence* the hits are the same feature - a
   deliberate retreat from the earlier proximity window, which *inferred* a
   provirus from nearness. **Consequence (accepted):** adjacent genes (gag/pol/env
   occupy distinct, non-overlapping positions) do **not** merge, so orphan loci are
   mostly **single-gene** - conservative deduplication rather than speculative
   multi-gene assembly. A per-genome **ground-truth cap** = the widest LTRdigest
   `LTR_retrotransposon` (`max(width(retrotransposons))`); a cluster wider than any
   real provirus is **flagged `oversized` (kept, not dropped)** for downstream
   filtering. (The initial revision used a `parameters.orphan_merge_gap` proximity
   window, since removed.)

2. **`fragment` -> `orphan` rename (full).** Every tier-referring identifier:
   `TRACK_FRAGMENTS_DIR`->`TRACK_ORPHANS_DIR` (+ `tracks/orphans/`), rule
   `taxonomy_fragments`->`taxonomy_orphans`, `{genome}.fragments.*`->`.orphans.*`,
   `fragments_counts`->`orphans_counts`, `--fragments_ranges`->`--orphans_ranges`,
   `load_fragments`->`load_orphans`, `fragment_recovery_*`->`orphan_recovery_*`.
   Biological "ORF fragment" language (fixtures, placement alignment) is kept.

3. **Unified authoritative catalog (`catalog.csv`).** `taxonomy_plot_generator.R`
   writes the union of anchored  U  orphan loci as one non-overlapping record set -
   the single "this is what we found at this location, and here is everything
   about it" table. `reconcile_catalog` resolves the only remaining cross-tier
   overlap with **anchored precedence**: where an orphan cluster's *span* bridges
   over an anchored provirus, the LTR-confirmed anchored locus wins and the orphan
   is dropped **from the catalog** (it stays in the per-genome `.orphans` table).

4. **Loss funnel** gains an honest orphan grouping step: orphan **hits ->
   orphans_total (clustered loci) -> orphans_recovered (classified)**, mirroring the
   anchored `valid -> loci_total -> loci_classified`.

## Consequences

- **Orphans become deduplicated, non-overlapping loci.** 0 orphan-orphan overlaps
  (verified on 5 genomes). With overlap-only clustering they are **mostly
  single-gene** (co-located redundant hits collapse; distinct adjacent genes stay
  separate) - the conservative, evidence-based choice. Clusters exceeding the
  per-genome max-provirus cap are `oversized`-flagged for filtering.
- **The catalog is fully non-overlapping** (0 union overlaps; anchored precedence
  dropped ~1-1.4% bridging orphans), with `source` keeping the confidence
  gradient explicit (**anchored = LTR-confirmed, orphan = proximity-inferred**).
- **Honesty caveat (biological judgment, recorded here):** even overlap-clustering
  *infers* that co-located orphan hits are the same feature where the anchored path
  had structural *evidence* (a shared LTR element). Overlap is a much stronger cue
  than proximity, and the ground-truth cap bounds implausibly-long clusters, but
  `source=orphan` still marks the whole tier as lower-confidence than
  LTR-confirmed proviruses.
- **Out of scope (per the maintainer):** hotspot, pair-detection, circle plots,
  and solo-LTR were intentionally *not* re-integrated in this change - they are a
  dedicated follow-up session.
