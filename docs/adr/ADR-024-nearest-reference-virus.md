# ADR-024: Nearest reference virus, with its identity

- **Status**: Accepted
- **Date**: 2026-09-30
- **Deciders**: Jorge González García
- **Extends**: [ADR-007](ADR-007-taxonomic-classification.md) (the locus call)

## Context

The classifier blasts every gene region of a locus against the reference
proteins and keeps only each hit's taxon, from which it makes the genus call. The
reference proteins name their source virus (`gag protein [Mouse mammary tumor
virus]`), but that name was dropped at this step, so no table after it could say
which known virus a locus resembles. By `catalog.csv` the virus names were gone.

A bare name would mislead. Measured on the model genomes (blastx of the
classifier's own regions, best hit per region, headline gene by
`gene_priority`), the median amino-acid identity of a locus to its nearest
reference virus is 40 to 47% in every genome. Most mouse loci "are" Simian
retrovirus 8 or Mason-Pfizer monkey virus at about 45%: distant relatives, not
instances. Only two sets are near-identical (at least 90%): 34 Desmodus loci to
the Desmodus rotundus endogenous retrovirus (median 99%), and 51 mouse loci to
Moloney murine leukemia virus plus one to mouse mammary tumor virus, the known
endogenous copies of those viruses. No human locus reaches 90%.

## Decision

Four new columns on the loci and orphan tables, and in `catalog.csv`:

| Column | Meaning |
|---|---|
| `nearest_virus` | the virus of the headline gene's best blastx hit |
| `nearest_virus_identity` | that hit's amino-acid identity, in percent |
| `nearest_virus_gene` | the headline gene |
| `per_gene_nearest` | every gene's nearest virus and identity, e.g. `GAG:Mouse mammary tumor virus(55.0);POL:...`; the identity is always the last parenthesis, since a virus name can hold one (`Feline sarcoma virus (STRAIN HARDY-ZUCKERMAN 4)`) |

- **Best hit** = highest bit score among hits to the reference, the same hit that
  gives the placement frame. Its identity is that alignment's; the strongest
  alignment of a long region is its long one, since bit score grows with length.
- **Headline gene** is the gene behind the locus's `taxon_call`, so the virus
  name and the call rest on the same evidence (a Betaretrovirus call from a
  placed GAG must not sit beside POL's nearest gammaretrovirus). When that gene
  has no reference hit, or the locus has no call, `classification.gene_priority`
  decides (ADR-022), then the strongest hit.
- **The virus** is the last `[bracket]` of the reference defline (an earlier one
  can hold a strain); a reference without one is named by its accession (none of
  the 383 current references lacks one).
- The identity always travels with the name. The columns are evidence beside
  the genus call, never an input to it.

## Consequences

- blastx writes one more field (`pident`); hits and calls are unchanged, so every
  existing column keeps its value. Loci without a reference hit leave the four
  columns blank.
- Taking the headline from the call's gene matters: on Desmodus it moves the
  headline for 95 of 406 loci away from the first gene in `gene_priority`. One of
  the 34 DrERV loci then reports its call gene's lower identity; its POL match
  stays visible in `per_gene_nearest`.
- A reader can now separate "this locus is a copy of a known endogenous virus"
  (DrERV, MLV in the mouse) from "this locus is a distant relative of the nearest
  reference" (nearly everything else).

## Alternatives considered

- **Only the name, no identity**: shorter, but it invites reading 45% as "is
  that virus".
- **One column per gene**: fixed width for a variable gene set; the packed
  `per_gene_nearest` follows the existing `per_gene` column instead.
