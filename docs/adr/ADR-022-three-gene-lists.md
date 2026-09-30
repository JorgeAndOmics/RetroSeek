# ADR-022: Three gene lists, one job each

- **Status**: Accepted
- **Date**: 2026-09-30
- **Deciders**: Jorge González García
- **Supersedes in part**: [ADR-007](ADR-007-taxonomic-classification.md) (gene
  reliability from the `main_probes` order) and
  [ADR-009](ADR-009-anchored-domain-tiering-and-structure-class.md) (the
  `canonical_order` column)

## Context

`parameters.main_probes` was one ordered list doing four jobs: which genes are
main, the number `completeness` divides by, the reliability order of the locus
call, and the expected genomic gene order. The last two want different orders.
Reliability wants `POL` first; a provirus reads `5'-gag-pro-pol-env-3'`. The
config told the user to "order for whichever column you intend to read".

Measured on the five model genomes with the POL-first list:

- `canonical_order` was `False` for 1,730 of 1,740 full LTR-flanked proviruses,
  nearly all of them textbook gag, pol, env.
- The column said nothing for most other loci: `False` for 6,584 loci with no
  main gene, and `True` by construction for any locus with one or two (any two
  genes match a list or its reverse).
- The check ignored strand, accepting either direction on either strand.
- Adding a gene for one job changed the others: a fourth main gene lowers
  `completeness` from 1.0 to 0.75 and turns `full` into `partial`.
- R and Python read the list differently (upper-casing in one, de-duplication in
  the other).

## Decision

**Three settings.**

| Setting | Job |
|---|---|
| `parameters.main_probes` | Membership: main versus accessory, the mosaic gene set, the completeness count. Order ignored. |
| `classification.gene_priority` | Reliability order of the locus call, most reliable first. |
| `parameters.gene_order` | The genes 5' to 3', any probe name, for `canonical_order` only. |

**`canonical_order` is strand-aware and has three answers.** The listed genes
present in a locus must lie in `gene_order` along the locus's strand (reversed
on the minus strand): `True` or `False`. With fewer than two listed genes, or a
strand vote that is tied, the column is blank.

**Old configs keep working.** A missing new list falls back to `main_probes` in
its written order, and the launcher's preflight warns. The preflight also stops
on a repeated or lower-case name in any of the three lists, and warns about a
name that is not a probe.

## Consequences

- With `gene_priority` equal to the old `main_probes` order, every taxon call is
  unchanged.
- `canonical_order` changes. On the model genomes: 1,722 full proviruses go from
  `False` to `True`; 37 two-gene loci and 10 three-gene loci that run against
  their strand or are scrambled go from `True` to `False`; 6,216 LTR-flanked loci
  and 29,640 orphans with nothing to check become blank.
- The locus strand is the majority strand of its hits. On the model genomes the
  hits of an LTR-flanked locus never disagree, and 99% of full proviruses follow
  gag, pol, env in that direction, so the strand is a sound reference.
- Readers of the column must handle the blank: the Parquet tables carry `""`,
  the catalog CSV an empty field that reads back as missing. The gene-order page
  draws only the loci that could be checked and says how many.
- The ranges manifest now records `main_probes`. The classifier's job log already
  records its full command line, including the three lists.

## Alternatives considered

- **Keep one list and document the trade-off**: the status quo; it left one of
  two columns wrong whichever order was chosen.
- **Hard-code the genomic order**: gene names are the user's probe names, so the
  pipeline cannot know them.
- **A fourth list for the genes a full provirus must have**: would let `PRO` be
  main without lowering completeness. Not needed yet; `main_probes` keeps that
  job.
