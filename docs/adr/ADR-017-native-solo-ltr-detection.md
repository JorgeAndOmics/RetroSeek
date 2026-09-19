# ADR-017: Native solo-LTR detection, using detected LTRs as bait

- **Status**: Accepted (thresholds provisional, pending the LTRharvest sweep)
- **Date**: 2026-09-19
- **Deciders**: Jorge González García
- **Supersedes**: [ADR-013](ADR-013-solo-ltrs-on-the-assembled-catalog.md) (archived at `origin/solo-ltr-v1-archive`)
- **Builds on**: [ADR-015](ADR-015-domain-evidence-by-symmetric-scan.md), [ADR-016](ADR-016-one-domain-classification.md)

## Context

Solo LTRs are the single-LTR remnants left when a provirus's two LTRs recombine
homologously and excise the internal region. In mammalian genomes they outnumber
intact proviruses by one to two orders of magnitude and each marks one ancestral
integration, so omitting them omits most of the history. LTRharvest cannot find
them at any parameter setting: its search primitive requires a **pair**.

ADR-013 obtained them from LTR_retriever's undocumented `solo_finder.pl`. That ran
and produced a plausible result for Desmodus rotundus: 9,504 solos against 406
intact loci, a 23:1 ratio, with the coordinate taxonomy join firing on 9,504 of
9,504. The detection worked. **The classification attached to it did not.**

Measured across every genome LTR_retriever ran on:

- **62 to 81% of its library is labelled `LTR/Gypsy`**, against 15 to 26%
  `LTR/Retrovirus`. Ty3/Gypsy elements are LTR retrotransposons but not
  retroviruses, so on its face that is a contamination finding.
- But the label does not track *env*. Pooled over all genomes, `LTR/Gypsy` entries
  sit on elements carrying *env* **58%** of the time against `LTR/Retrovirus`'s
  **33%**, and `LTR/Copia` entries carry it **62%** of the time even though Copia
  elements have no *env* at all.
- Our own join was ruled out first: each library entry was mapped to the single
  LTRharvest element whose span it best matches, and that mapping checked. On
  Desmodus, 32 entries matched an element exactly and 45 more at Jaccard >= 0.9,
  median Jaccard 1.00.

The mechanism is documented. LTR_retriever's README lists **TEsorter** as a
dependency, so its superfamily call is TEsorter's. TEsorter classifies on the
`GAG-PROT-RT-RH-INT` backbone and has **no ENV model**; neither does DANTE. Both
were built for plants. Mammalian ERVs are phylogenetically *inside* the Ty3/Gypsy
lineage, having evolved from gypsy-like elements by acquiring *env*, so a
classifier asking "Copia or Gypsy?" from RT/RH/INT is asking a question whose
correct answer for a mammalian ERV is largely "Gypsy". RepeatMasker's own
documentation concedes the same limit: it is "especially hard to predict if an LTR
is derived from an endogenous retrovirus or a non-autonomous LTR element."

**No available tool distinguishes retroviral from Gypsy LTR elements in mammals.**
Swapping detectors would not fix this. What RetroSeek has that they do not is *env*
evidence (ADR-015) and a map of surviving retroviral coding sequence in both tiers.

## Decision

**Detect solo LTRs natively, using the LTR arms RetroSeek already extracts as
bait, and subtract everything that is not a solo. No new dependency.**

### The idea

An intact element's LTR arm is a sequence exemplar of that family's LTR. RetroSeek
already writes every arm to `results/tracks/flanking_ltr/`. Use them as blastn bait
to find other copies of the same LTR elsewhere in the genome. This is the move the
pipeline already makes one level up, where the probe CSV holds retroviral proteins
used as tBLASTn bait; here the bait is nucleotide and one level down.

### The three fates that must be separated

Every LTR began as one of a pair flanking a provirus.

| class | second LTR | internal region | solo? |
|---|---|---|---|
| flanking | present, still pairable | present | no, the element is intact |
| monoLTR at an orphan | lost or too diverged to pair | **present** | no, the coding sequence survives |
| solo | recombined with | **excised** | **yes** |

The middle class is the one **LTR_retriever structurally cannot exclude**, because
it subtracts only its own intact elements and has no concept of an orphan.
RetroSeek's orphan tier is precisely a map of surviving retroviral coding sequence
that LTRharvest missed.

### Criteria

1. matches an LTR arm of an element hosting a catalogued ERV locus (retroviral by
   construction, so no classifier is involved)
2. does not overlap any LTRharvest element (excludes class 1)
3. is not within a proviral distance of an orphan locus (excludes class 2)
4. alignment covers 0.8 to 1.2 of the bait arm and is at least 80 bp
5. **bait and hit are both at least 300 bp, at >= 95% identity**

Criteria 1, 2 and 4 are Ou & Jiang's published rules applied as written. Criterion 3
is the one only RetroSeek can apply. Criterion 5 is ours, and is required: with the
published rules alone the method returns 199,816 candidates, a 492:1 ratio that is
an order of magnitude too high. The reason is that Ou & Jiang's coverage rule is
against a family **consensus**, which is full length by construction, while our bait
is individual arms running 102 to 999 bp, where 80% of a short arm is meaningless.

### Validation

At the chosen operating point the method returns **11,197 solos, 27.6:1**, against
LTR_retriever's **9,504 and 23:1**. The two share no code, no library and no
external tool and land within about 18% of each other.

A phylogeny of 1,212 LTR sequences across all three classes (MAFFT, IQ-TREE
`GTR+G`) carries two results:

- **Control**: an element's two arms were identical at insertion, so they must be
  sister tips. 284 of 406 elements (70.0%) recover them as sisters.
- **Solos form their own clades**: 87.0% of tips have a same-class sister against a
  permutation null of 61.0% (sd 1.4). Solos are enriched ~2.2x beside other solos
  and **depleted** (~0.55x) beside flanking arms. There are LTR families in this
  genome existing predominantly or entirely as solos, with no intact representative
  for LTRharvest to find, which is evidence the method reaches genuinely new
  material rather than re-finding what we already had.

Full working, with every calibration curve and both negative results, is in
`notebooks/solo_ltr_native_method.Rmd`.

## What this gives up

**Sensitivity to ancient solos.** Criterion 5's identity requirement is an age
filter: it keeps solos whose family still has a close modern relative. That is the
same bias as LTRharvest's `-similar 85`, which is itself an age ceiling at roughly
34 My (15% divergence at ~2.2e-9 substitutions/site/year). A curated Dfam family
consensus sits closer to the ancestral sequence than any surviving copy does and
would reach further back. That is the one capability no native implementation
matches, and testing it is blocked on the installed Dfam being partition 0 only:
213 curated families for the whole mammalian ancestry.

**Target-site duplication does not rescue it.** Recombination happens between the
LTRs and leaves the TSD intact, so a flanking direct repeat should be positive
evidence. Measured against a matched random null it is enriched 2.23x, confirming
the candidates are not noise, but the absolute rate is 1.7% against 0.8%. A 4 to
6 bp exact repeat does not survive tens of millions of years. DANTE_LTR reaches the
same conclusion, grading solos `SL` versus `SL_noTSD` rather than requiring one.

**Independent family names.** No `MLT1A` or `LTR12C`; families are our own.

## Consequences

- **`ltr_retriever` leaves `environment.yml`**, and with it `repeatmasker` and
  `rmblast`. That retires **gotcha 43**, the standing pin conflict where rmblast
  installs its own `blastn`/`tblastn`/`makeblastdb` over blast's copies and forces
  both to be pinned in lockstep at 2.14.1. It also removes whole-genome
  RepeatMasker from the critical path, which ADR-013 listed as its own worst
  negative consequence.
- **Evidence is recorded, not gated** (the ADR-015 principle). Distance to the
  nearest orphan is a column, not a hard filter, because in repeat-dense regions a
  genuine solo can sit near an unrelated orphan by chance.
- **Thresholds are provisional.** An LTRharvest `-similar` sweep is running to test
  whether older bait lifts the identity ceiling. `-similar 85` has completed and
  reproduced the production element count exactly (9,893 for Desmodus), which is
  the control that sweep needed to pass.
- **Every number here is Desmodus only.** A second genome is required before this
  becomes pipeline code.

## Alternatives considered

**Keep LTR_retriever and ignore its classification.** Rejected: the classification
is what the dependency was bought for, `solo_finder.pl` is undocumented, and
RepeatMasker stays on the critical path.

**"One code to find them all"** (Bailly-Bechet, Haudry & Lerat 2014, *Mobile DNA*
5:13). The only solo detector benchmarked on human, and it sidesteps classification
entirely by assembling RepeatMasker fragments and reporting what will not assemble.
Rejected for now: it needs a real Dfam download, a per-genome RepeatMasker pass, and
a hand-curated LTR-to-internal dictionary that the paper is explicit cannot be
automated for human. Their curated human dictionary is from the 2012 library.
Revisit if ancient solos become the binding constraint.

**`dante_ltr_solo`** (Novák et al. 2024, *NAR Genomics and Bioinformatics* 6:113).
Purpose-built, TSD-validated, emits the solo/complete ratio directly. Rejected:
plant-benchmarked, self-labelled work in progress, and installs from a personal
conda channel rather than bioconda.

## Revisit trigger

- Ancient solos become the binding constraint, making the Dfam route worth its cost.
- The sweep shows relaxed-similarity bait materially lifts the identity ceiling, in
  which case criterion 5 should be re-derived rather than inherited.
- A second genome fails to reproduce the convergence with LTR_retriever.

## References

- `notebooks/solo_ltr_native_method.Rmd` - the full calibration and validation.
- Ou & Jiang 2018, *Plant Physiology* 176:1410. doi:10.1104/pp.17.01310.
- Bailly-Bechet, Haudry & Lerat 2014, *Mobile DNA* 5:13. doi:10.1186/1759-8753-5-13.
- Novák et al. 2024, *NAR Genomics and Bioinformatics* 6:113.
- `origin/solo-ltr-v1-archive` - the superseded LTR_retriever implementation.
- Gotcha 60 - LTR_retriever output traps, preserved from that branch.
