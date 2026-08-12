# ADR-014 - Publishing the placement evidence, and co-phylogeny from it

- **Status**: Accepted
- **Date**: 2026-08-12
- **Deciders**: Jorge Gonzalez Garcia
- **Builds on**: [ADR-007](ADR-007-taxonomic-classification.md),
  [ADR-008](ADR-008-rank-agnostic-classification.md),
  [ADR-011](ADR-011-phylogeny-aware-plots-and-rank-segmentation.md)

## Context

RetroSeek published exactly two tree figures and one reusable reference
phylogeny. Everything else phylogenetic was computed on every run and deleted.

`taxonomy_classify` and `taxonomy_orphans` each run EPA-ng, which writes
`epa_result.jplace` and `labelled_tree.newick` into
`data/tmp/{taxonomy,orphans}/{genome}/place_POL/` - a directory the data-layout
contract describes as "scratch; cleared between runs". Measured on the model 5:

| tier | files | queries placed |
|---|---|---|
| LTR-flanked | 5 | 9,492 |
| orphan | 5 | 14,212 |

All ten carry an **identical reference tree** (md5 `cc42886999de`, 64 leaves /
125 branches), which is the precondition for comparing samples to one another,
and it holds exactly.

That jplace is the evidence behind every `taxon_call` in the catalog: which
branch of the retroviral phylogeny each locus attached to, with likelihood
weights across candidate branches. The catalog kept the conclusion and discarded
the evidence. `.jplace` is also the interchange format iTOL and gappa read, so
"show me why this locus is Betaretrovirus" had no answer.

`gappa` was already installed and already invoked (`gappa examine assign`), so
its entire comparative suite was reachable at zero marginal cost.

## Decision

### Publish the evidence

`taxonomy_placement.export_placement` copies the artifacts to
`results/tracks/taxonomy/placements/{genome}.{tier}.{gene}.jplace` (+
`.labelled.newick`), declared as rule outputs.

Placement legitimately does not run in three cases - the gene has no tree
package, every query aligned to all-gaps, or no locus carried that gene - so a
**valid empty jplace carrying the real reference tree** is synthesised instead.
A rule declaring the output would otherwise break the DAG, and `tree_layout.py`
already sets this precedent by writing header-only CSVs so the DAG holds.

The reference tree is carried through verbatim rather than replaced by a
placeholder: downstream gappa commands read it to know what they are drawing.
Edge numbering is generated for the synthetic file; gappa reads it back as 125
branches / 64 leaves, matching what it derives from real EPA-ng output.

**The aggregate rules must demand these files.** Without that, Snakemake never
schedules the setup rule when its table outputs are current, and the jplace
stays in scratch. This is why `--classify` re-runs classification once.

### Heat-trees, not grafts

`gappa examine graft` draws one pendant edge per query and is unreadable past a
few hundred; *Mus musculus* places 5,690 in the LTR-flanked tier alone.
Heat-trees accumulate mass onto branches and stay legible at any count. SVG for
editing, Newick for reuse, Nexus because FigTree opens it directly.

**Every gappa call is gated on the placement count.** `heat-tree` does not
return quietly on a placement-free file: it aborts with `Invalid Color
Normalization with min >= max` and dumps core, because it builds a colour scale
from an empty mass range. An empty-state SVG naming the sample and the reason is
written instead, mirroring `empty_plot()` in the R generators.

### EDPL as a second uncertainty axis

Expected Distance between Placement Locations measures how far apart a query's
candidate placements are. It is independent of the `confidence` column: a locus
with confidence 1.000 and a high EDPL is one whose certainty is an artifact of
the LCA collapsing genuinely scattered placements.

Measured on *Antrozous* (LTR-flanked): 98% of placements are tight (EDPL < 0.1,
median 0.000), with three loci at 2.1-2.6. The tail is small and specific, which
is exactly what makes it worth having.

### Co-phylogeny

`gappa analyze squash` builds a tree of the **host genomes** from their
placement distributions; `krd` gives the distance matrix behind it. Congruence
against the host phylogeny is measured on **bipartitions (Robinson-Foulds),
never on branch lengths** - the supplied host tree is frequently a cladogram
with placeholder lengths, and comparing those to placement distances would be
meaningless. `tree_layout.py` now warns when branch lengths are uninformative
(all absent, or all whole numbers) so a topological comparison is never mistaken
for a dated one.

Inputs are staged under genome-derived names first. gappa labels samples by file
basename and EPA-ng writes every run to `epa_result.jplace`, so passed directly
all five tips come back named `epa_result` - a silently useless tree.

Tiers are analysed separately: whether the weaker orphan tier tells the same
story as the LTR-confirmed one is itself a check.

### Catalog-side composition panels

`tree_confidence_plot` was already generic over (tree, catalog column) but only
plotted confidence. `tree_composition_plot` is its counts sibling, driving
`species_composition_tree` (lineage composition per host, ordered by host
relatedness) and `taxon_tier_tree` (LTR-flanked vs orphan per lineage - a
lineage that is almost all orphan is one whose structural evidence has eroded).

No gene-discordance tree: `mosaic_gene_discordance_plot` already unpacks the
same per-gene calls.

## Consequences

- **Positive**: the evidence behind every taxon call is reviewable and
  iTOL-ready; placement uncertainty becomes a reportable quantity; the
  co-phylogeny question is answerable from data already produced; no new
  dependency.
- **Negative, measured not hidden**: demanding the placement outputs makes
  `--classify` re-run classification once for all genomes and tiers (blastx +
  mafft + epa-ng; hours, not days). This is a one-time cost - afterwards the
  files are cached like any other output.
- **Coupling**: `--placement-trees` depends on the classification stage.

## Measured result

Both tiers are **discordant with the host phylogeny** (RF = 4, 0 of 2 splits
shared), and running them separately is what makes the result interpretable:

```
ltr-flanked:  (Homo, (Molossus, (Mus, (Antrozous, Desmodus))))
orphan:       (Homo, ((Antrozous, Mus), (Desmodus, Molossus)))
host:         (((Antrozous, Molossus), Desmodus), (Homo, Mus))
```

**What replicates:** *Homo sapiens* is the outgroup in both ERV trees, on a long
branch (0.65 / 0.80). The host phylogeny pairs it with *Mus*. So the human
POL-placement profile is distinct from both the mouse and the bats, and that
survives changing the evidence tier.

**What does not:** the internal structure among the other four genomes differs
completely between tiers. Whatever groups *Antrozous* with *Desmodus* under
LTR-flanked evidence groups it with *Mus* under orphan evidence.

The KRD matrix argues the signal is compositional rather than a sample-size
artifact: on the LTR-flanked tier the closest pair by query count (*Homo* 2,039
/ *Molossus* 786) is the **most distant** by KRD (1.231), while the wildly
mismatched *Mus* 5,690 / *Antrozous* 672 is among the closest (0.626).

**Reported as a finding to investigate, not a conclusion.** Only the
Homo-outgroup result is stable across tiers; the rest should not be cited.
Placement is on POL alone against a 64-tip reference, which bounds how much
lineage structure any of this can resolve - and is the strongest argument for
widening `placement_genes`, since three genes would give the comparison
independent replicates rather than one draw.

## Alternatives considered

**Leave the artifacts in `data/tmp` and copy them ad hoc.** Rejected: the files
are the evidence for published calls, and reproducibility means the pipeline
produces them, not that someone remembers to rescue them before the next run.

**Avoid the re-classification by migrating the existing scratch files.** They
were present and valid at implementation time, so a one-off copy would have
worked. Rejected as unreproducible: it would leave the DAG asserting outputs
nothing in the pipeline had produced.

**Widen `placement_genes` to GAG and ENV.** Deferred, not rejected. Three
independent placements per locus would turn `is_mosaic` from an LCA
disagreement into genuine phylogenetic discordance, but unlike everything here
it costs two more IQ-TREE reference builds plus per-genome alignment and
placement.

## Revisit trigger

- Placement is widened beyond POL; the per-gene trees make the mosaic question
  answerable properly and the co-phylogeny gains independent replicates.
- A dated host timetree is supplied, making the co-phylogeny quantitative
  rather than topological.
- gappa changes its sample-naming or file-prefix behaviour.

## References

- `workflow/scripts/taxonomy/taxonomy_placement.py` - `export_placement`.
- `workflow/scripts/taxonomy/placement_figures.py` - heat-tree, EDPL, LWR.
- `workflow/scripts/taxonomy/placement_cophylogeny.py` - squash, KRD, congruence.
- `workflow/scripts/taxonomy/tree_layout.py` - branch-length warning.
- Matsen, Kodner & Armbrust 2010, pplacer. *BMC Bioinformatics* 11:538.
- Czech, Barbera & Stamatakis 2020, genesis and gappa. *Bioinformatics* 36:3263-3265.
