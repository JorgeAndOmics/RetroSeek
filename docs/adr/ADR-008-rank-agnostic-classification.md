# ADR-008: Rank-agnostic taxonomic classification (axis-taxon calls)

- **Status**: Accepted
- **Date**: 2026-07-01
- **Deciders**: Jorge González García
- **Refines**: [ADR-007](ADR-007-taxonomic-classification.md) (per-locus classification stage)

## Context

ADR-007 replaced best-bitscore probe-label transfer with calibrated, self-evidenced
per-locus calls (placement + weighted-LCA) - a real improvement. But it did so by
**hard-committing the classifier to the genus rank**, which quietly discarded the
rank-agnostic contract the probeset `Label` column was designed for:

- The reference axis was a **hard-coded list of retroviral genera** (`taxonomy_reference_builder.GENERA`).
- The hierarchy was **trimmed to the `Retroviridae` subtree** (`taxonomy_build_hierarchy.ROOT`), so any non-retroviral lineage was thrown away.
- "Resolved" was literally `rank == "genus"` in the assembler and the placement guard.
- Consequently any `Label` value that is not a retroviral genus (`Bornaviridae` = family, `Metaviridae` = LTR-retrotransposon family) was **silently dropped** (`taxonomy_lca.lca` ignores labels absent from the taxonomy) and could never be a first-class call.

The probeset `Label` was always deliberately **rank-flexible** - an opaque taxon tag at
whatever rank the probe represents. The classifier's genus-rooting broke that: it quantised
all identity to the genus grid and its ancestors, resolving only *upward* (genus -> subfamily
-> family) and never to a rank the retroviral genus-tree cannot express.

Two different "agnosticisms" were being conflated. ADR-007 achieved **probe-agnosticism**
(the classifier does not trust the probe's label; it re-derives identity from an independent
reference) but lost **rank-agnosticism** (identity may live at any rank).

## Decision

Generalise the classifier's atomic unit from *genus* to **axis taxon** - a seed taxon at
*whatever rank the reference is declared at*. The classification engine (`lca`, `ancestors`,
`rank_of`, `weighted_lca`, placement, mosaic) is already rank-general; only four
genus-specific pins are removed:

1. **Axis is declared, not hard-coded (hybrid).** The reference is built for a set of seed
   taxa at any rank, resolved as: `classification.reference_taxa` (config override) if set,
   else the **distinct probeset `Label` values**. This reconnects `Label` to the classifier
   as the *axis declaration*. The reference CSV column `genus` becomes `taxon`.
2. **The hierarchy root trim is removed.** `taxonomy_build_hierarchy` keeps each taxon's
   full NCBI lineage; `lca` finds the true common ancestor. Mixed-rank axes (e.g.
   `{Lentivirus, Bornaviridae}`) can now share one tree.
3. **"Resolved" becomes `node  in  axis`**, replacing `rank == "genus"` in the placement guard
   and the confident-call set. A locus resolves when its call lands on a seed taxon (at the
   axis's rank); otherwise it is an honest LCA-backoff to an interior ancestor, reported
   with its true `rank`.
4. **The output column `genus_call` is renamed `taxon_call`** (the existing `rank` column
   already carries genus/subfamily/family/...). Input `Label` (rank-agnostic axis seed) and
   output `taxon_call` + `rank` (resolved identity, or honest backoff) are now distinct,
   correctly-named objects.

### Behaviour is preserved by default

With the default axis (the current retroviral genera) the change is **byte-identical**: the
hierarchy's only genus-rank nodes *are* the axis genera, so `rank == "genus" <=> node  in  axis`.
`erv_class` remains an overlay keyed on the call (`ERV_CLASS.get(taxon_call, "")` -> `""`
off-Retroviridae). Root un-trimming adds nodes above `Retroviridae` that are never reached
when the axis is retroviral.

### Alternatives considered

- **Keep `genus_call`, generalise only the meaning** - rejected: the name lies whenever the
  axis is not genus-rank.
- **Add `taxon_call` as an alias of `genus_call`** - rejected: two columns for one value
  violates DRY.
- **Config-only axis (no Label link)** - rejected: leaves `Label` disconnected; the hybrid
  keeps the decoupling option (`reference_taxa`) while making `Label` the natural default.

## Consequences

**Positive**

- **True rank-agnosticism.** A `Bornaviridae`-seeded axis yields a first-class family-rank
  `taxon_call` instead of `unclassified`; the pipeline resolves to whatever rank you declare.
- **`Label` reconnected**, honestly: it seeds the axis, it is not the per-locus truth.
- **All ADR-007 work preserved** (placement, weighted-LCA, confidence, mosaic, fragments,
  loss, plots); defaults reproduce every count.

**Negative / costs**

- The reference CSV schema (`genus` -> `taxon`) and the output schema (`genus_call` ->
  `taxon_call`) change, so the reference must be **rebuilt** and loci **reclassified**;
  downstream R plots, tests, and committed demo figures are updated to the new column.
- Un-trimming the hierarchy means backoff can now report ranks above family (e.g. order)
  when a mixed-rank axis spans them - intended, but a new possibility to interpret.

## References

- Supersedes in part [ADR-007](ADR-007-taxonomic-classification.md) (genus-rooting).
- `docs/architecture.md`, `docs/configuration.md` (`classification.reference_taxa`),
  `docs/taxonomy_classification/`.
