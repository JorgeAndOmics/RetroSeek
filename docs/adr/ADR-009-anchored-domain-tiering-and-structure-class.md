# ADR-009: Anchored domain-tiering, structural ERV class, and orphan relabeling

- **Status**: Accepted
- **Date**: 2026-07-09
- **Deciders**: Jorge González García
- **Refines**: [ADR-007](ADR-007-taxonomic-classification.md) (per-locus assembly), [ADR-008](ADR-008-rank-agnostic-classification.md)

> **Terminology note (2026-07-23).** What this ADR calls the **anchored** tier is
> now named **`ltr-flanked`** throughout the code, data values and current docs -
> see [ADR-011](ADR-011-phylogeny-aware-plots-and-rank-segmentation.md). This file
> keeps the original wording (and its filename) because an ADR records what was
> decided at the time; read "anchored" here as "ltr-flanked".

## Context

The candidate->valid step (`ranges/validation.R::find_valid_hits`) **discarded**
every LTR-anchored tBLASTn hit whose gene label did not match a Pfam domain in
its enclosing element. That conflated two independent kinds of evidence -
**structural anchoring** (the hit sits inside an LTRdigest element) and **domain
corroboration** (a curated Pfam HMM fired for that gene) - and silently lost
divergent/degraded ERV genes whose HMMs simply did not cross threshold: exactly
the old, diverged integrations the pipeline exists to find. A latent footgun
compounded it: with an empty `domains:` config the valid tier came out **empty**,
not "pass-through".

Separately, the per-locus assembly (ADR-007) carried only **continuous**
structural metrics (`completeness`, `canonical_order`, `n_main_genes`) - it never
committed a locus to a discrete full / partial / single-gene class, the natural
catalogue unit for an ERV survey.

## Decision

Stop discarding; start **labelling**. Separate the conflated evidence into two
orthogonal axes and add a discrete structural class.

1. **Valid = the whole anchored set.** `find_valid_hits` (filter) becomes
   `annotate_anchored_hits` (label): every candidate is kept, annotated, and
   exported. `is_valid` (boolean) is retired.

2. **Anchoring axis - `source`**: `anchored` (inside an element) vs `orphan`
   (outside every element; renamed from `fragment` - "orphan" reflects lost
   proviral context without implying measured degradation). Downstream `source`
   splits and the loss-funnel `orphan` branch follow.

3. **Domain evidence, two grains** (only meaningful within `anchored`, since
   LTRdigest domains are children of elements):
   - **`domain_tier`** (per provirus / locus, element-wise): `domain_selected`
     (>=1 config-matched domain of *any* gene) > `domain_unlisted` (has protein
     domains, none config-matched) > `non_domain` (no domain). Strongest tier
     wins across straddled elements.
   - **`domain_hit_class`** (per hit): `substring_match` / `no_substring_match`
     (+ `non_domain` in positional mode), controlled by
     `parameters.hit_domain_mode` (`membership` default = the old `valid_mask`
     co-occurrence test; `positional` = co-localization, reusing the existing
     `feature_class` overlap). The per-hit flag rides the valid GFF3 (IGV); the
     per-locus tier rides the loci table + plots.

4. **Structural class - `structure_class`** (per locus): `gene`
   (`n_main_genes <= 1`) > `full` (`completeness >= classification.structure_full_min`,
   default `1.0`) > `partial` (multi-gene, below the floor). Deliberately
   gene-content only - flanking-LTR structure stays in the anchoring axis and the
   solo-LTR module ([ADR-003](ADR-003-ltr-retriever-pre-filter.md),
   [ADR-005](ADR-005-ltr-retriever-runner.md)), so LTR evidence is not
   double-encoded.

Both new per-locus columns (`domain_tier`, `structure_class`) join
`LOCI_COLUMNS`, ride the IGV GFF3, and drive two new taxonomy plots
(`domain_tier_composition`, `structure_class_composition`) plus an
`erv_like_structure_class` panel. `HC`/`LC` confidence is unchanged.

## Consequences

- **Recall preserved with provenance.** Formerly-discarded anchored hits survive
  as `domain_unlisted` / `non_domain`, tagged so confidence is judged downstream
  rather than by silent deletion. `domain_selected` counts tick slightly above
  the old `valid` count (element-wise vs gene-wise gate).
- **Hotspot enrichment + solo-LTR prefilter now consume the full anchored set**
  (they read `valid_ranges.gff3`). This widens their input versus the old
  domain-filtered valid tier - an accepted behavior change; filter on
  `domain_tier == "domain_selected"` to recover the prior scope.
- **Config surface**: `parameters.hit_domain_mode` and
  `classification.structure_full_min` added (3-way config/schema/docs sync).
- **Empty `domains:`** now yields all-`non_domain` anchored loci (still exported)
  instead of an empty valid tier.
- On-disk `fragments/` track paths keep their historical names (label-only
  rename); a full path rename is a possible follow-up.
