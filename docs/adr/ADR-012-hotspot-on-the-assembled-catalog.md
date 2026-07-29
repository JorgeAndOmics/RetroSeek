# ADR-012 - Hotspot detection on the assembled ERV catalog

- **Status**: Accepted
- **Date**: 2026-07-29
- **Deciders**: Jorge González García
- **Builds on**: [ADR-008](ADR-008-rank-agnostic-classification.md), [ADR-010](ADR-010-orphan-clustering-and-authoritative-catalog.md), [ADR-011](ADR-011-phylogeny-aware-plots-and-rank-segmentation.md)

## Context

Hotspot detection was the last stage still speaking the pre-redesign vocabulary.
It was built when the pipeline's unit of analysis was the tBLASTn *hit*; since
ADR-007..011 the authoritative unit is the *locus*: one non-overlapping
integration event carrying `taxon_call`, `segment`, `structure_class`,
`domain_tier`, `confidence` and `source`.

The statistics were never the problem. N-masking, genome tiling, the
`MASS::glm.nb` fit with its `offset(log(effective_bp))`, the deliberate refusal
to fall back to Poisson, BH correction and region merging are all sound and are
**unchanged by this ADR**. The *input contract* was wrong in four ways:

1. **It counted hits, not integration events.** The reduced valid track carries
   4,106-32,800 features per genome against 905-7,270 loci. A GAG+POL+ENV
   provirus contributes 3+ hits and a single-gene fragment 1, so window counts
   were weighted by **gene content**. Two clusters of equal integration density
   scored differently purely because one held more complete proviruses - which
   is not what "integration hotspot" means.
2. **`group_split` split by `label`**, the legacy probe-derived tag that
   ADR-007/008 replaced with `taxon_call`. The redesign keeps `label` only as
   detection provenance, so the stage stratified by the thing it superseded.
3. **ADR-009 silently inflated the input**: `valid` became the whole LTR-flanked
   set (all domain tiers) where it had been domain-filtered. Accepted then,
   never resolved.
4. **`structure_class`, `segment`, `source` and `confidence` were unreachable**
   from a GFF3 track, and the ~30k orphan loci were invisible.

## Decision

### Count integration events, from the catalog
`hotspot.input` becomes `catalog` (default) or `original`:

- **`catalog`** reads `catalog.csv` - the non-overlapping per-locus assembly -
  and counts one row per integration event. This makes hotspot a **consumer of
  the classification stage**, the coupling ADR-009 deferred: a stale catalog now
  pulls `--classify` into the DAG.
- **`original`** (raw tBLASTn hits) is kept deliberately. Its density is what
  gives the NB model power on sparse assemblies, and it remains the right tool
  for genome-wide calling where taxonomic resolution is not the question.
- **`valid` is retired**: it is the same events as `catalog`, multi-counted.

The catalog keys on the config display name while the rule wildcard is the genome
stem, so `load_catalog_loci()` matches through the `species:` map,
separator- and case-insensitively, and **aborts when a genome matches zero rows** -
silence there is indistinguishable from "this genome has no ERVs".

### Annotate composition, do not fragment the model
`structure_class` describes an *event*, not a window, so it cannot be a window
covariate. Splitting the model by class was rejected: the count matrix is already
sparse and the NB already fails to converge on low-count genomes. Instead,
detection runs on **every** locus in the tier and each called region is annotated
with what it contains - `n_full` / `n_partial` / `n_gene`, `n_ltr_flanked` /
`n_orphan`, `dominant_taxon`, `mean_confidence` - which rides into the GFF3
attributes, a new per-region `{genome}.hotspots.csv`, and a
`{genome}_composition.pdf` panel. A hotspot is then readable as
intact-provirus-driven or fragment-driven without a second query, and no region
is added, dropped or re-scored by the annotation.

### Group by the calibrated call
`hotspot.group_by` = `segment` (default) | `taxon_call` | `none`, replacing the
`group_split` bool. Rank-agnostic: `segment` follows
`classification.segment_rank`, so switching the study to family-level grouping is
a config change. When the column is absent - notably under `input: original` -
the run logs a warning and pools, because a missing covariate must not abort a
detection run.

### Select the tier explicitly
`hotspot.source` = `ltr-flanked` (default) | `orphan` | `both`. The default
answers "where do LTR-confirmed proviruses cluster"; `both` is available when
density matters more than tier purity.

## Consequences

- **Positive**: counts are integration events; hotspots carry biological
  composition; grouping uses the authoritative taxonomy at any rank; the orphan
  tier is reachable; the retired `valid` tier removes a silent bias.
- **Negative, measured not hidden**: counting events is far sparser than counting
  hits. On the 5 model genomes at the default tier and 500 kb windows, only
  *Mus musculus* (7,270 events) yielded hotspots - 2 regions; the other four
  yielded none. That is the honest consequence of counting integrations instead
  of genes, not a regression to tune away. The levers are documented rather than
  applied by default: `source: both` for density, a larger `window_size`, or
  `input: original` for genome-wide calling.

  Measured head-to-head, same genomes and window size:

  | genome | tier | events | groups | regions called |
  |---|---|---|---|---|
  | *Mus musculus* | `catalog` | 7,270 | 7 taxa | **2** |
  | *Mus musculus* | `original` | 218,370 | 1 (pooled) | 1 |
  | *Antrozous pallidus* | `catalog` | 905 | 6 taxa | 0 |
  | *Antrozous pallidus* | `original` | 21,931 | 1 (pooled) | 0 |

  Sparser input did **not** mean less sensitive. On *Mus* the catalog tier found
  MORE regions despite 30x fewer events: 218k hits over 5,499 windows push the NB
  baseline so high that a genuine cluster struggles to exceed it, and raw hits
  carry no taxonomic call, so everything pools into one model where a
  lineage-specific cluster is diluted. On *Antrozous* both tiers agree on zero,
  which is the useful control - the empty result is a property of that genome at
  500 kb windows, not an artifact of counting events.
- **Coupling**: `--hotspot-detection` now depends on the classification stage
  under the default tier.

## Verification

`make check` (ruff, mypy strict, 239 pytest, 552 testthat). End-to-end on the 5
model genomes with cached tBLASTn/GenomeTools: event counts equal the catalog
LTR-flanked locus counts exactly (905 / 406 / 1,069 / 2,499 / 7,270), groups are
taxa (`Alpharetrovirus`, ..., `unassigned_at_genus`) rather than probe labels,
and the *Mus* regions carry composition (26 loci: 5 full, 11 partial, 10 gene).

**Also fixed here**: `hotspot_detector.R` still sourced
`scripts/range_analysis/exporters.R`, a directory renamed to `ranges/` in the
scripts reorg (`6b3ba8e`). The stage had been unrunnable since then; nothing
caught it because its integration test is skipped and every unit test sources the
modules directly.
