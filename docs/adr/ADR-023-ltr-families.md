# ADR-023: LTR families from the solo-LTR bait arms

- **Status**: Accepted
- **Date**: 2026-09-30
- **Deciders**: Jorge González García
- **Builds on**: [ADR-017](ADR-017-native-solo-ltr-detection.md) (the bait arms and the seed of each solo)

## Context

A solo LTR is tied to one seed arm: the arm of an intact, ERV-bearing element it
matches at 95% identity or more. Nothing grouped arms into families, so there was
no solos-per-intact-copy figure per family, no way to ask which families are old,
and no basis for finding solos too diverged to match a single arm. The only
"families" were clades cut from each genome's sampled evidence tree, used by
nothing else.

## Decision

**Families are groups of bait arms at least `solo_ltr.families.identity`
identical over the whole shorter arm** (default 0.8, the field's 80% convention),
found by `cd-hit-est` with both strands compared (bait arms are not oriented) and
each arm placed in its most similar family. Every element joins the family of its
arms; every solo inherits the family of its seed arm, as it inherits a genus.

**Identifiers** are `<code>_F001`, numbered by size, ties broken by the
representative arm's position, so the same input gives the same names. The code is
one genus letter and three species letters (`Mmus`, `Mmol`, `Hsap`). Two genomes
can share a code (see Consequences). The column is `ltr_family`, to keep it apart from the
evidence tree's `family` clades.

**Measured on the model genomes before choosing the default** (16,430 arms, all
five genomes):

| Identity | Families | Single-arm families | Arm pairs together | Genus purity |
|---|---|---|---|---|
| 0.80 | 820 | 29 | 97.0% | 94.4% |
| 0.85 | 1,108 | 56 | 96.5% | 94.6% |
| 0.90 | 2,091 | 517 | 82.5% | 95.7% |
| 0.95 | 5,126 | 3,344 | 59.4% | 97.4% |

"Arm pairs together" is the method's own positive control: an element's two arms
were identical at insertion. Counting identity over 80% of the shorter arm instead
of all of it gave nearly the same families with slightly fewer pairs together
(95.0% at 0.80), so the simpler whole-arm rule with one setting was kept.

**Pooled families** (`Pool_F001`...) cluster every genome's arms at once, in
sorted genome order so the config's order cannot change them, so a family shared
by several genomes is one family.

**Optional Dfam labels** (`solo_ltr.families.dfam`, off by default). The families
have neutral identifiers, so a reader cannot tell that `Mmus_F001` is IAPLTR1_Mm.
When switched on, each family's representative arm is searched with every curated
Dfam model (`nhmmer`, E-value at most 1e-5 per model), and the best bit score
names the family in `ltr_family_dfam.csv`. Not Dfam's gathering thresholds: in
4.0 about 5,600 of the 30,646 curated models have no model-level GA line, which
stops `nhmmer --cut_ga`, and where there is one it equals the strict TC, while the
per-taxon thresholds sit in TH lines. Over the whole search 1e-5 expects well under
one chance hit. The models are fetched like Pfam:
a pinned release (`input.dfam_release`, 4.0), its own guarded stage
(`--download-dfam`), checked against Dfam's MD5 and kept beside `Pfam-A.hmm`.
Curated models only: the uncurated set is mostly raw repeat-finder output, the
kind of label the label is meant to check.

## Consequences

- On the five model genomes no family spans two genomes at 0.8: the 820 pooled
  families are exactly the per-genome ones. A direct `blastn` of each genome's
  arms against the others agrees. Only a few dozen arms have any significant
  cross-genome hit, and the best reach about 81% over the whole arm. The LTR
  families of these genomes are lineage-specific. The pooled tables will matter
  for closely related genomes, such as congeneric species.

- Family identifiers are unique within a genome, not across genomes: two genomes
  can share a code (Canis_lupus_familiaris and Canis_lupus_dingo are both Clup).
  Every table that holds several genomes (the pooled tables, the Dfam labels)
  carries the genome, so nothing stops or merges over it.
- New per-genome tables: `{genome}.ltr_families.csv`, `.ltr_family_summary.csv`,
  `.ltr_family_genus.csv`, and a manifest recording the identity used. Identifiers
  depend on it, so tables from runs with different settings must not be compared
  by identifier.
- The summary carries each family's genus mix and median LTRdigest arm-pair
  similarity (an age signal). Families sit inside genera (about 94% purity), and a
  mixed family is worth a look: a wrong genus call, a recombinant, or a
  contaminant.
- The Dfam label is evidence, never a filter. Dfam 4.0's curated models apply to
  about 1,400 families each for human and mouse but about 760 to 810 for the three
  bats, mostly through families shared across mammals, so bat families will more
  often stay blank. The download is 1.7 GB and unpacks to about 11 GB.

## Alternatives considered

- **All-against-all `blastn` plus our own grouping rule**: no new dependency, but a
  few dozen lines of our own clustering to test and maintain, and a less familiar
  method to a reader. `cd-hit-est` is the usual tool for exactly this 80% rule, runs
  in seconds and is deterministic; it is declared in `environment.yml`.
- **Cutting the evidence tree**: the tree is a sample of about 950 tips per genome,
  not every arm.
