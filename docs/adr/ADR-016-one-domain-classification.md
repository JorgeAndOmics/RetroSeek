# ADR-016: One domain classification, at locus grain

- **Status**: Accepted
- **Date**: 2026-09-18
- **Deciders**: Jorge González García
- **Refines**: [ADR-009](ADR-009-anchored-domain-tiering-and-structure-class.md), [ADR-015](ADR-015-domain-evidence-by-symmetric-scan.md)

## Context

ADR-015 introduced a domain scan that reads Pfam accessions and bit scores and applies
`data/config/pfam_domain_classes.tsv` at **locus** grain, for both tiers. It left the existing
`ranges` stage applying the **same curated table** to `gt ltrdigest`'s domain **names** at
**element** grain.

The result was two columns called `domain_tier`, in two tables, from two sources, answering two
questions:

| where | classifies | grain |
|---|---|---|
| catalog (`*.loci.csv`, `*.orphans.csv`) | the scan's accessions | locus |
| stage hits table + valid GFF3 track | LTRdigest's names | element |

They are not required to agree, and legitimately will not: an element (mean 6,919 bp) can carry a
domain in a region no catalogued locus (mean 2,996 bp) covers. A reader comparing the two would
reasonably assume otherwise. The name-keyed side is also the fragile one, since Pfam guarantees
accession stability but not names: 22 of the 6,175 names in the current tracks are already absent
from Pfam 38.2.

A second, related question was whether LTRdigest's domain search could simply be retired, which
would have removed the duplicate at the source. Measurement said no, and the reasons are recorded
below because they were not obvious.

## Decision

**The curated table is applied exactly once: in the domain scan, on accessions, at locus grain.**

1. The `ranges` stage stops classifying. `extract_all_domains` returns LTRdigest's `protein_match`
   features unclassified; `load_domain_classes` and `extract_selected_domains` are deleted.
2. `annotate_ltr_flanked_hits` emits **`Parent` only**. The element-grain `domain_tier` is gone from
   the valid track and from the stage hits table.
3. What LTRdigest still contributes at element grain is a **count**, not a judgement:
   `n_domains_total`, plus `element_domains`, a raw "; "-joined set of Pfam names. Deliberately not
   called `domain_names`, which is the catalog's locus-grain column.
4. `build_stage_probe_domain_df` reports `hit_probe` against `domain_name`, the raw Pfam name,
   rather than a class.
5. The two stage funnel plots lose their `domain_selected` tier and run `homology -> candidate`.
   The equivalent domain view lives in the catalog, which already emits
   `loci_domain_selected` / `loci_domain_unlisted` / `loci_non_domain` and the
   `domain_tier_composition` plot.
6. `defaults.PFAM_DOMAIN_CLASSES` becomes the single resolved path for the table, replacing three
   hardcoded copies in the Snakefile and a fourth read in `ranges/io.R`.

**`gt ltrdigest -hmms` is kept.** It is not redundant with the scan, for a reason that only showed
up under measurement.

## Why LTRdigest's domain search stays

Measured on the model genomes:

- **`-hmms` is essentially all of LTRdigest's 24-hour runtime.** PPT detection is a Viterbi decode
  over a 61 bp window; without the Pfam search the tool runs in minutes. Dropping it looked like a
  free saving.
- **But it is the pipeline's only accurate strand source.** LTRharvest assigns none at all (26,499
  of 26,499 Antrozous elements come out `?`). LTRdigest derives element strand from its domain
  hits, and scored against an independent physical signal (purine content of the polypurine tract,
  which owes nothing to either rule) it is right **93.7%** of the time. A bit-score vote over our
  own scan managed 87.7%, and two variants scored lower still. Against the majority strand of the
  tBLASTn hits inside an element, LTRdigest was right 403/403 and PPT orientation 144/161.
- Dropping `-hmms` would leave 26 to 47% of elements unoriented, and `findOverlaps(...,
  ignore.strand = FALSE)` treats `*` as a wildcard (`ranges/validation.R:35`), so those elements
  absorb hits on either strand: 1.1 to 2.5% of orphan hits migrate tiers. Removing LTRdigest
  entirely costs 1.5 to 3.4% and the polypurine tract as well.
- 12 R files read `protein_match`. Removing it degrades them silently rather than loudly.

A caution for anyone revisiting this: **PPT strand as written in the GFF3 is not independent
evidence.** LTRdigest forces the chosen strand onto every child feature, so 100% of `RR_tract`
strands equal their parent element's. Scoring a strand rule against it is circular. Use the tract's
purine content instead.

The route worth trying if the runtime ever needs reclaiming is **PBS**: `gt ltrdigest -trnas` takes
a tRNA library and is currently not passed one, so PBS emits zero features. It is LTRdigest's
designed strand evidence and would be independent of domains.

## The `candidate` tier goes with it

`find_candidate_hits` output was exported as its own track, counted as its own
funnel stage, and carried as its own variable. Once ADR-009 stopped `valid` from
filtering, `valid` became `candidate` plus a `Parent` attribute, and the two were
identical in every genome:

| genome | `candidate_ranges` | `valid_ranges` |
|---|---|---|
| Homo_sapiens | 36,348 | 36,348 |
| Mus_musculus | 153,248 | 153,248 |
| Antrozous_pallidus | 11,202 | 11,202 |

So `tracks/candidates/`, the `candidate_ranges` counters, the `candidate` column
of the overlap matrix and the duplicated funnel step are removed.
`build_stage_hits_df` takes one set instead of two, and `is_candidate` keeps its
name because it still says something true: the hit overlaps an LTR element.

The track is NOT renamed here. `valid` is a poor name for a set that no longer
excludes anything, and `element_hits` is the intended replacement, but
`valid_ranges` has 24 receivers and four of them belong to the solo-LTR
workstream, whose branch is unmerged and already edits the same lines. The
rename lands with that work, in one pass instead of two.

## Consequences

- One curated judgement, in one place, on the stable key. Editing the table changes the catalog and
  nothing else.
- `domain_tier` now means exactly one thing. The stage tables no longer carry a same-named column
  that can disagree with it.
- Two stage plots lose a tier and the structural panel loses one component. The information is not
  lost, it moved to the catalog where it is computed symmetrically across both tiers.
- LTRdigest keeps its 24 hours, which is paid once per genome added.
- Stop codons now translate to `X` rather than `*`, chosen by benchmark over `esl-translate`'s
  ORF-splitting and BATH: it covers 2,192 loci against 1,946 and 1,973 respectively, and 1,866
  retroviral-diagnostic loci against 1,557 and 1,653.
- **BATH v2.0.0-rc4** was benchmarked and not adopted: it has no `--cut_ga` so thresholds cannot be
  matched, and it is a pre-release with no conda recipe. It wins on depth (10,138 locus-family
  pairs against 9,316) and loses on breadth (2,069 loci against 2,192). Revisit at a stable release.
