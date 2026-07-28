# ADR-011 — Phylogeny-aware plots, rank segmentation, and the `ltr-flanked` rename

- **Status**: Accepted
- **Date**: 2026-07-23
- **Supersedes (terminology only)**: [ADR-009](ADR-009-anchored-domain-tiering-and-structure-class.md)
- **Builds on**: [ADR-008](ADR-008-rank-agnostic-classification.md), [ADR-010](ADR-010-orphan-clustering-and-authoritative-catalog.md)

## Context

Running the pipeline on 102 bat genomes exposed limits that 5 model genomes never
reached:

1. **Plots were built for ~5 genomes.** `auto_dims()` existed and the older
   plot2sort/stage generators used it, but the taxonomy, structure and loss
   generators never called it — their `emit()` passed only the fixed config
   width/height, so every ADR-009/010 panel rendered on a 15x12in canvas
   regardless of genome count. Even where it did run, the canvas grew but the
   text did not, so labels still collided.
2. **No phylogenetic context.** At n=102 a per-species bar chart ordered
   alphabetically is not just cramped, it is the wrong question: the reader wants
   to know whether ERV burden tracks host relatedness, and which viral lineages
   dominate where.
3. **No way to segment results by genus** — a natural unit for reporting, but one
   that must not re-hard-code a rank the pipeline deliberately abstracted away.
4. **"Anchored" named the method, not the object.** The tier is defined by being
   flanked by LTRs; "anchored" said only that it attached to something.

## Decision

### 1. Two trees, from two deliberately different sources
- **Taxon tree = a cladogram derived from the reference's `taxonomy.tsv`.**
  It covers every axis taxon by construction, needs no representative-sequence
  choice, and is rank-agnostic. It is a **classification, not a phylogeny**, and
  therefore carries **no branch lengths** — stated plainly on the panel.
  *Rejected*: collapsing the real `POL.contree` ML phylogeny to one tip per taxon.
  It has real branch lengths, but requires picking an arbitrary representative
  accession (or an MRCA policy) and breaks on paraphyletic taxa — a defensible
  phylogeny bought with an indefensible choice.
- **Species tree = a user-supplied Newick pinned by `input.species_tree`.**
  Reproducible with no run-time network, and free to carry real divergence times
  if the user exports a dated tree. Unset is a normal state: the panels render an
  explanatory placeholder rather than failing.
  *Deferred*: an Open Tree of Life fetch helper. The API works and returns a real
  induced subtree, but auto-fetching at plot time makes topology drift between
  runs. The pinned-file contract is designed so a build-once fetch rule can write
  that same file later.

### 2. Zero new dependencies — Python lays out, R draws
`tree_layout.py` uses **Bio.Phylo** (biopython is already pinned for the tBLASTn
cache) to parse, prune and ladderize, then writes flat `x/y` segment and tip
CSVs. The R panels draw them with `geom_segment` and align tree-to-bars with
`patchwork`. No `ape`, `ggtree`, `treeio` or `aplot` — a prototype proved the
coordinate bridge is sufficient, and the environment stays lean.

Layout is **deterministic**: ladderized, with sibling sets ordered by a stable
name key, so identical inputs give byte-identical coordinates.

### 3. Rank segmentation that stays rank-agnostic
`segment_of(taxon_call, segment_rank)` walks the reference hierarchy (reusing the
existing `taxonomy_lca.ancestors()` / `rank_of()`) and returns the first ancestor
at the requested rank. `classification.segment_rank` is **any** NCBI rank, and no
taxon name is hard-coded — set it to `family` and the same code segments by
family.

A call **coarser** than the requested rank (e.g. `Retroviridae` when segmenting by
genus) becomes **`unassigned_at_<rank>`**. This is the biologically honest
outcome: ADR-008 exists precisely so a locus is not given precision its evidence
does not support, and segmentation must not smuggle that back in.

`segment` is derived **in the classifier**, so it rides the loci table -> GFF3 ->
catalog like `structure_class` and `oversized` before it. The `taxonomy_segments`
stage is then a pure split plus a **curated** 3-plot subset per segment — the full
20-plot panel per segment would be hundreds of PNGs for little gain.

### 4. `anchored` -> `ltr-flanked`
Data values, labels and `--source` use `ltr-flanked`; identifiers use
`ltr_flanked`; prose uses `LTR-flanked`. `find_unanchored_hits` becomes
`find_orphan_hits` — it returns the orphan tier, and a naive substitution would
have produced `unltr_flanked`.

**Forward-only.** Catalogs already written carry `source=anchored` until the run
is repeated. **ADR-009's filename and the ADR/progress-log prose keep the original
wording**: they are point-in-time records, and rewriting them would falsify the
history they exist to preserve.

### 5. Plot scaling
`scale_categorical_axis()` attaches the canvas `auto_dims()` computes **and**
applies matching text sizing/rotation (`categorical_text_size` shrinks past the
base with a 5pt legibility floor; x ticks rotate at >=20 categories). `auto_dims`
keeps its signature, so the 12 existing call sites are untouched. New
`plots.per_stratum` / `plots.max_dim` knobs. Studies at or below the base strata
render byte-identically to before.

## Consequences

- **Positive**: 102-genome panels are readable; results carry phylogenetic
  context; genus-level reporting exists without breaking rank-agnosticism; the
  vocabulary describes the biology; no new dependencies.
- **Negative / accepted**: the taxon tree shows relationships without distances
  (a classification); the species tree is the user's responsibility to supply and
  pin; the rename means old and new catalogs disagree on one column value until
  a re-run; per-segment plots are a curated subset, not the full panel.
- **Species-name coupling**: tree tips are matched to the config `species:`
  display names (case- and separator-insensitive), and unmatched names are
  reported in **both** directions rather than silently dropped. This is why
  `erv_like_plot_generator`'s missing `relabel_species()` call was fixed in the
  same arc — mismatched labels would have silently emptied the species panel.

## Verification

`make check` (ruff, mypy strict, 232 pytest, 525 testthat); Snakemake dry-run
resolves `taxonomy_tree_layout` and `taxonomy_segments`; `tree_layout.py` run
against the real 5-genome reference produces a 13-tip mixed-rank taxon cladogram
(family + subfamilies + genera as sibling tips) and a 5-tip host tree, and both
panels render. Scaling asserted synthetically at n=102, which the 5 model genomes
cannot exercise.
