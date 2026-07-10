# ADR-010: Orphan proximity-clustering, `fragment`→`orphan` rename, and the authoritative catalog

- **Status**: Accepted
- **Date**: 2026-07-10
- **Deciders**: Jorge González García
- **Refines**: [ADR-009](ADR-009-anchored-domain-tiering-and-structure-class.md) (anchored tiering + orphan tier)

## Context

After ADR-009, anchored proviruses were a clean per-locus catalog — one row per
LTR element (grouped by `Parent`), non-overlapping, each with an authoritative
`taxon_call` + `structure_class` + `is_mosaic`. **Orphans were not.** An orphan
(a non-LTR-associated hit) has no LTR element, so `build_loci` made **each
parentless gene-hit its own singleton locus** — meaning orphans were always
single-gene (`structure_class = gene`, never mosaic) and could **spatially
overlap** each other (a POL orphan and a GAG orphan on the same degraded provirus
became two overlapping rows: ~1,478 overlaps for *Antrozous*).

Two consequences: (1) there was no single, non-overlapping, authoritative locus
table across both tiers; (2) the vocabulary was split-brain — the concept was
"orphan" (the `source` value, loss metrics, plot labels since ADR-009) but the
plumbing still said "fragments" (dirs, files, rule names, CLI args, functions).

## Decision

1. **Orphan proximity-clustering (synthetic `Parent`).** `cluster_orphan_hits`
   (ranges/validation.R) clusters orphan hits by spatial proximity
   (`GenomicRanges::reduce`, `min.gapwidth = parameters.orphan_merge_gap + 1`,
   strand-ignored) and stamps each cluster with a synthetic `Parent`
   (`orphan_<seqname>_<clusterStart>`). The classifier's existing
   `build_loci`/`_assemble` then group orphans into **single, non-overlapping,
   multi-gene loci — verbatim the anchored path**, so orphans gain
   `structure_class` and mosaic detection for free. New config
   `parameters.orphan_merge_gap` (default **10000** bp ≈ a full provirus span,
   ~7–12 kb; also matches the solo-LTR `nearest_erv_max_distance` precedent).

2. **`fragment` → `orphan` rename (full).** Every tier-referring identifier:
   `TRACK_FRAGMENTS_DIR`→`TRACK_ORPHANS_DIR` (+ `tracks/orphans/`), rule
   `taxonomy_fragments`→`taxonomy_orphans`, `{genome}.fragments.*`→`.orphans.*`,
   `fragments_counts`→`orphans_counts`, `--fragments_ranges`→`--orphans_ranges`,
   `load_fragments`→`load_orphans`, `fragment_recovery_*`→`orphan_recovery_*`.
   Biological "ORF fragment" language (fixtures, placement alignment) is kept.

3. **Unified authoritative catalog (`catalog.csv`).** `taxonomy_plot_generator.R`
   writes the union of anchored ∪ orphan loci as one non-overlapping record set —
   the single "this is what we found at this location, and here is everything
   about it" table. `reconcile_catalog` resolves the only remaining cross-tier
   overlap with **anchored precedence**: where an orphan cluster's *span* bridges
   over an anchored provirus, the LTR-confirmed anchored locus wins and the orphan
   is dropped **from the catalog** (it stays in the per-genome `.orphans` table).

4. **Loss funnel** gains an honest orphan grouping step: orphan **hits →
   orphans_total (clustered loci) → orphans_recovered (classified)**, mirroring the
   anchored `valid → loci_total → loci_classified`.

## Consequences

- **Orphans become first-class loci.** 0 orphan-orphan overlaps (verified on 5
  genomes); many orphan loci are now multi-gene (e.g. Mus 2,176) and mosaic
  (Mus 581) — recombinant signal that singleton orphans could never express.
- **The catalog is fully non-overlapping** (0 union overlaps; anchored precedence
  dropped ~1–1.4% bridging orphans), with `source` keeping the confidence
  gradient explicit (**anchored = LTR-confirmed, orphan = proximity-inferred**).
- **Honesty caveat (biological judgment, recorded here):** orphan clustering
  *infers* a provirus unit from proximity where the anchored path had structural
  *evidence* (a shared LTR element). Nearby orphan hits may be one degraded
  provirus, or independent insertions. `source=orphan` + the tunable
  `orphan_merge_gap` keep this uncertainty visible and adjustable; the default is
  deliberately conservative (one provirus span).
- **Out of scope (per the maintainer):** hotspot, pair-detection, circle plots,
  and solo-LTR were intentionally *not* re-integrated in this change — they are a
  dedicated follow-up session.
