# ADR-025: A placement weight floor for the gene call

- **Status**: Accepted
- **Date**: 2026-10-02
- **Deciders**: Jorge González García
- **Amends**: [ADR-007](ADR-007-taxonomic-classification.md) and
  [ADR-008](ADR-008-rank-agnostic-classification.md) (when a placement wins)

## Context

Each gene region of a locus is called two ways: phylogenetic placement on the
gene's reference tree, and weighted-LCA over its blastx hits. Until now any
placement that landed on an axis taxon won, however weak it was. The locus then
reports the placement's gene first (placement outranks `gene_priority`).

Measured on the model genomes, a placement's weight (aLWR, the share of
likelihood weight on the placed lineage) predicts whether it agrees with blastx
on the same DNA:

| POL placements | Loci | Placement and blastx differ at genus |
|---|---|---|
| aLWR 0.8 or below | 2,904 | 750 (19.6%) |
| aLWR above 0.8 | 3,824 | 26 (0.7%) |

GAG shows the same trend. The weak placements drive most of the rare-lineage
calls: Delta (68 loci over the five genomes) rests mostly on single GAG, ENV or
POL placements with median aLWR 0.13 to 0.62, and 205 of Molossus's 277 Epsilon
calls are orphan GAG placements with median aLWR 0.76. Gene priority itself
(POL first) is sound: only 1.5% of resolved loci are outvoted by their other
genes, and several probes cover the same DNA, so a vote would count it twice.

## Decision

A new setting, `classification.placement_min_weight` (0 to 1). A placement wins
over the gene's blastx call only when its aLWR is at least this value (inclusive);
below it the gene is called by weighted-LCA, exactly as when no placement exists.

- **Default 0**: every on-axis placement wins, as before. Existing configs and
  results are unchanged.
- **Measured choice 0.8**: on the model five, Delta 68 to 34, Epsilon 919 to
  782, Alpha 110 to 106; 173 loci lose their genus because their blastx call stops
  above genus rank.
- The setting acts per gene, in `_gene_call`, so everything downstream (the gene
  behind the call, `nearest_virus*` per ADR-024, the HC/LC tag) follows the call
  that survives.

## What it does not fix

Confident but wrong placements. 88 of the 156 Lentivirus calls are ENV regions
placed with aLWR near 1 on Lentivirus branches, because six Betaretrovirus ENVs
sit inside the Lentivirus clade of the ENV reference tree. A weight floor cannot
see that; dropping ENV from `placement_genes` or requiring blastx agreement for
ENV placements can. That is a separate decision.

## Consequences

- One more classification setting (`config.yaml`, `schema.yaml`,
  `docs/configuration.md`) and CLI flag (`--placement-min-weight`).
- Adding the parameter changes the classification rules' recorded params: the
  first run after this change reruns classification (hours) unless run with
  `--rerun-triggers mtime`.
- A lower floor keeps more genus calls but more of them disagree with blastx;
  `method` and `confidence` in the loci tables still show which evidence won.
