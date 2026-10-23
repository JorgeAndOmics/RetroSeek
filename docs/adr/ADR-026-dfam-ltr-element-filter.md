# ADR-026: Drop LTR elements whose arms are LINE, SINE or DNA, using Dfam

- **Status**: Accepted
- **Date**: 2026-10-22
- **Deciders**: Jorge González García
- **Builds on**: [ADR-023](ADR-023-ltr-families.md) (LTR families, Dfam labels of families),
  [ADR-017](ADR-017-native-solo-ltr-detection.md) (the solo-LTR bait) and
  [ADR-010](ADR-010-orphan-clustering-and-authoritative-catalog.md) (the orphan tier that receives the dropped elements' hits)

## Context

LTRharvest pairs any two similar stretches into an "LTR element". Two copies of an L1
3' end or two Alus around an ERV probe hit make an element that looks ERV-bearing, and
its arms then become the solo-LTR bait, so the solo search counts every L1 or Alu copy
in the genome as a solo LTR. Measured on the model genomes (backlog 16, 2026-09-30):
83 to 98% of all solos came from families whose arms Dfam calls non-LTR; after R7 the
human genome had 16,042 LINE and 5,305 SINE solos against 699 ERV-LTR ones.

The first fix is not in LTRharvest: `-motif tgca -motifmis 1 -vic 10` keeps only 43 to
63% of real ERV elements and still 11 to 29% of LINE/SINE-armed ones, and a TSD
requirement changes nothing (every default element already has a chance-level 4 bp
TSD). A solos-per-intact ratio filter is wrong too, since real MaLR families are mostly
solo. The fix belongs at the bait: know what the arms are.

Dfam labels of arms are reliable (measurements, 2026-10-06, 2,500 sampled arms):
Dfam's non-LTR calls match UCSC RepeatMasker for 75 of 75 arms, no ERV LTR in 906 was
called non-LTR, and the two arms of an element reach the same LINE/SINE-or-LTR verdict
in 1,249 of 1,250 elements. Elements with a non-LTR arm are 7.2% (human) to 17.2%
(Molossus) of the ERV-bearing ones, about 1,100 of 11,709 across the model genomes.

## Decision

**After LTRdigest, label both arms of every LTR element that overlaps a probe hit with
Dfam, and drop the element from the LTR set when either arm is a non-LTR repeat.**
Arm classes `LTR/*` (ERV or not: ERV1, ERVK, ERVL, Gypsy, Copia, MaLR) stay; arms with
no Dfam label stay; anything else (LINE, SINE, DNA, satellite, ...) drops the element.
The element's ERV hits then overlap no element and become orphans, classified by their
own sequence like every orphan.

- **Where.** A new rule, `ltr_arm_filter_setup`, between LTRdigest and
  `ranges_analysis`; `ranges_analysis.R` removes the dropped elements (with their arms,
  protein matches, TSDs and container) right after loading, before hits are sorted into
  LTR-flanked and orphan. No other stage changes; they simply see fewer elements.
- **Which elements.** Those overlapping any raw probe hit, whatever the strand: a
  superset of the ERV-bearing ones (10 to 14% more: 449 against 406 in Desmodus, 1,208
  against 1,069 in Molossus, 8,253 against 7,270 in mouse). A superset is safe, since a
  verdict for an element that ends with no hit is never read. The catalog does not exist
  yet at this point, so it cannot name the ERV-bearing set.
- **How the arms are labelled.** By family, as in ADR-023: the arms are clustered with
  cd-hit-est at `solo_ltr.families.identity`, one representative per family is searched
  with every curated Dfam model (nhmmer, E 1e-5, the same code as the solo labels), and
  each arm takes its family's class. That is a few hundred nhmmer queries per genome;
  labelling every arm would be about 23,400 for the model genomes (about 22 hours,
  extrapolated from the measurement). The family label holds for its members in the
  measurement: 1,428 of 1,432 arms Dfam calls ERV LTR belong to ERV LTR families, and
  every LINE/SINE arm to a LINE/SINE family.
- **One switch.** `solo_ltr.families.dfam`, the setting that already turns on the Dfam
  labels of the solo stage. With it off, the rule is not in the DAG, `ranges_analysis`
  gets no flag and nothing is filtered; the launcher warns once when a run reaches
  ranges analysis with Dfam off. No new CLI flag.
- **Evidence kept.** `results/tables/ltr_arm_filter/{genome}.ltr_arm_filter.csv`: per
  element, each arm's family and class, the verdict and the reason.

## Result on the model 5 (2026-10-22, rerun from ranges analysis, 4 h 52 min, 0 failed)

| Genome | Probe-hit elements | Dropped | LTR-flanked loci | Orphans | Solos |
|---|---|---|---|---|---|
| Antrozous pallidus | 1,154 | 202 | 902 to 779 | 2,757 to 2,979 | 27,246 to 753 |
| Desmodus rotundus | 449 | 64 | 406 to 367 | 1,605 to 1,679 | 11,197 to 1,888 |
| Homo sapiens | 2,287 | 255 | 2,062 to 1,909 | 5,773 to 6,088 | 22,047 to 700 |
| Molossus molossus | 1,208 | 246 | 1,069 to 899 | 3,552 to 3,925 | 22,992 to 1,320 |
| Mus musculus | 8,253 | 647 | 7,270 to 6,799 | 16,074 to 16,977 | 67,539 to 9,583 |

Solos fell by 83 to 97%, the share the 2026-09-30 measurement attributed to families with
non-LTR arms (97.2, 83.1, 98.0, 94.3 and 84.9% predicted). Desmodus: 64 of 449 elements
dropped, 62 with an L1 arm and 2 with a DNA/hAT arm; family-level labels agreed with direct
per-arm labels for all 250 elements checked both ways. Three genomes warn that their
evidence tree is less tight than expected (86 to 89% of solos within 0.1 of their seed arm,
90% expected) with the much smaller solo sets; the tables are unaffected.

## Consequences

- Positive:
  - Solo-LTR counts stop including L1 and Alu copies; the solo stage needs no filter of
    its own, and the per-class reporting of backlog 16 becomes a check, not a patch.
  - Fewer fake LTR-flanked loci; their ERV hits are still reported, as orphans.
  - Dfam off leaves the pipeline exactly as it was.
- Negative:
  - Dfam must be downloaded (11 GB) for the filter to run; the stage fires only with
    Dfam on, as the labels already did.
  - A new compute step of minutes to about an hour per genome (nhmmer over 30,646 models).
  - Labelling by family can mislabel a mixed family: 4 of 1,432 arms in the measurement.
  - An element whose arms match no Dfam model stays, even if fake: the filter removes
    what Dfam can name, not everything LTRharvest got wrong.
  - Catalog counts change (fewer LTR-flanked loci, more orphans), so everything built on
    the old catalog is rerun from `ranges_analysis`: the notebooks of the validation
    folder, the solo tables and plots, the hotspot and placement stages.
- Neutral:
  - `solo_ltr.families.dfam` now means "Dfam is in use", for labels and filter alike.
  - The hotspot and ranges counts record `ltr_elements_dropped_dfam`.

## Alternatives considered

- **Label every arm with nhmmer** (the measurement's procedure). Exact, but about 22
  hours for the model genomes and weeks for 102 genomes.
- **Label one arm per element.** Half the cost, about one element in 1,250 different
  from the either-arm rule, but a departure from the written design for a saving the
  family route already gives.
- **Report per Dfam class and filter nothing** (backlog 16, 2026-10-01). Kept as the
  fallback when Dfam is off; it leaves the contamination in every count.
- **LTRharvest settings** (motif, TSD). Measured and ruled out, above.
- **A ratio filter on solos per intact copy.** Removes real MaLR families.
- **A curated LINE/SINE library inside the pipeline.** A second source of truth beside
  Dfam, which RepeatMasker also uses.

## Revisit trigger

If a new Dfam release changes class names or the arm agreement with RepeatMasker falls,
rerun the reliability measurement. If the mixed-family rate rises above a few percent,
label the arms of the families that mix individually.

## References

- Backlog 16 (solo contamination), `measurements/2026-10-06/dfam_reliability/REPORT.md`
  and `ltrharvest_motif/RESULTS.md` (local), gotcha 94.
- `workflow/scripts/ltr_arm_filter.py`, `workflow/scripts/ranges/granges_build.R`
  (`remove_ltr_elements`), `workflow/scripts/ranges/io.R` (`load_dropped_elements`).
