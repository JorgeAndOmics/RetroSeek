# What was built, what it found, and why it's trustworthy

Three explanations of the same work at increasing depth, then a section on why the
findings reliably reflect real biology. This describes the **implemented and tested**
classifier and its **empirical results** on the 5 model genomes — companion to the design
theory in `../taxonomy_classification_redesign.md` and the raw numbers in `04_FINDINGS.md`.

> Scope note. Genome-wide runs for mouse/human/Molossus characterise the *composition* of
> retroviral-homologous regions (a relative readout), not an exhaustive ERV census. The
> two bats with assembly-matched tracks (Desmodus, Antrozous) are classified at the
> domain-validated locus level. Both modes are described below.

---
---

# PART 1 — Plain language

Earlier we diagnosed why RetroSeek could say *where* a viral fossil (ERV) sits but not
reliably *what* it is: it picked the single best-scoring reference probe, and because the
key viral gene is similar across virus groups, that "winner" was close to a coin toss.

We built the fix and ran it on five mammal genomes (three bats, mouse, human). Instead of
copying a probe's label, the new method reads each fossil's **own** sequence, compares it
to a panel of known retroviruses, and — crucially — when several virus groups match, it
reports the **smallest group they all belong to** (e.g. "this is clearly a gamma-type" or,
honestly, "this is a retrovirus but I can't pin the genus") rather than forcing a guess.

What we found:
- On the bats, the method confidently named the genus for **~92%** of fossils, and the
  picture matches what biologists have published: **gamma-type (Class I) fossils dominate,
  beta-type (Class II) are a strong second, everything else is rare.**
- It automatically **corrected a known mislabel**: 23 vampire-bat fossils the old probes
  tagged "gamma" are really "beta" — exactly matching a published discovery of a
  *betaretrovirus* in vampire bats.
- The **mouse** genome came out **beta-dominant**, the opposite lean from bats — which is
  correct, because the mouse's most abundant ERVs (IAP, MusD) are beta-type. The method
  reproduced a real species difference, not a canned answer.
- **Human** came out strongly **gamma-type (Class I)** — the opposite lean from mouse,
  and exactly right: the human genome's most numerous ERVs (HERV-H, HERV-W) are Class I.
  The method reproduced the human↔mouse difference without being told the species.
- **Molossus** (third bat) matched the other two bats: **gamma-dominant, beta second.**

In short: the method reads the fossil itself, says how sure it is, admits when it can't
tell, and its answers line up with the published biology — including correcting a mistake
the old approach made.

---
---

# PART 2 — Intermediate

**The engine.** For each ERV locus we take its own genomic sequence (not the probe that
found it), search it with DIAMOND against an **independent**, genus-balanced reference of
retroviral proteins (built fresh from RefSeq — *not* the study probes, so reclassification
is non-circular), and assign the **bitscore-weighted lowest common ancestor (LCA)** of the
strong hits over the ICTV retroviral taxonomy. Only hits within 10% of the best score count
toward the LCA, so a locus with a clear best genus resolves to that genus, while a locus
with genuinely balanced cross-genus support backs off to subfamily/family or `unclassified`.
This is the MEGAN/detectEVE paradigm and is probe-agnostic by construction.

**Two run modes.** (A) *Locus reclassification* on the two bats whose valid tracks match
the testing-genome assembly (Desmodus via a RefSeq→GenBank seqname remap; Antrozous
natively). (B) *Genome-wide composition* via windowed DIAMOND for mouse/human/Molossus,
where no assembly-matched track exists and whole-genome tblastn was infeasible.

**Findings.**
- *Engine accuracy (leave-one-out):* 145/147 = **98.6%** correct genus on held-out
  reference proteins, high across **all** markers (GAG/PRO 100%, ENV/POL 98%) — quantitative
  evidence the method is correct and marker-agnostic.
- *Bats (Mode A):* Desmodus 91.7% genus-resolved, **Gamma 422 / Beta 274**; Antrozous 93.7%,
  **Gamma 767 / Beta 486**. Class I dominant, Class II strong second, rest trace.
- *DrERV correction:* 23 loci the probes called Gammaretrovirus are independently
  Betaretrovirus — matching the literature betaretrovirus in *Desmodus*. Robust to the
  score threshold (Gamma:Beta 440:292 / 422:274 / 377:261 at 5/10/20%).
- *Mouse (Mode B):* **Beta 442 / Gamma 325** — beta-dominant, matching IAP/ETn/MusD
  (Class II) abundance. Opposite lean from bats = species-specific signal.
- *Human (Mode B):* 96.8% genus, **Gamma 689 / Beta 186** (≈3.7:1) — strongly Class I,
  matching HERV-H/HERV-W numerical dominance; the inverse of mouse.
- *Molossus (Mode B):* 98.2% genus, **Gamma 378 / Beta 144** — bat-like, like the other two.

**Limits.** Mode B is composition, not abundance; the reference is exogenous-heavy and
currently lacks spuma (Class III undercounted, e.g. MuERV-L/HERV-L); exact proportions
depend on reference breadth. The *class-level shape* and *direction* are the robust claims.

---
---

# PART 3 — Technical

## Implementation (all in `workflow/scripts/`, unit-tested in `tests/unit/test_taxonomy_lca.py`)
- `taxonomy_lca.py` — curated ICTV retroviral taxonomy (parent map + ranks + ERV-class
  map); pure primitives `ancestors`, `lca`, `rank_of`, `weighted_lca`; GFF3 evidence parser.
  `weighted_lca(hits, top_percent=0.10)` keeps hits within `top_percent` of the best
  bitscore, sums score-mass per genus, returns `(lca_node, dominant_genus_fraction)`.
- `taxonomy_reference_builder.py` — Entrez RefSeq fetch per genus, gene-binned by defline,
  capped per (genus,gene) for balance → pinned `data/taxonomy_reference/{retro_reference.faa,csv}`
  (151 proteins, 6 genera; spuma absent — ICTV genus rename, known gap).
- `taxonomy_classify_loci.py` — Mode A: `bedtools getfasta` (strand-aware) → `diamond blastx`
  vs independent reference → per-locus `weighted_lca`.
- `taxonomy_classify_genome.py` — Mode B: 5 Mb windows → `diamond blastx` → coordinate
  remap → sweep-merge hits into loci → `weighted_lca`. (Windowing required: DIAMOND blastx
  returns nothing on chromosome-length queries.)
- `tblastn_to_loci.py` — helper to merge tblastn HSPs into candidate loci (for assemblies
  with no valid track, when tblastn is affordable).

## Validation suite
1. **13 unit tests** pin the taxonomy/LCA contract (single genus→genus; cross-genus→correct
   backoff; all-7→family; weighted near-tie→subfamily; abstention).
2. **Leave-one-out** on the reference: 98.6% genus accuracy, all markers (ground-truthed on
   RefSeq labels; self-hits excluded).
3. **Mode-A bats**: literature-consistent composition + locus-level probe-mislabel
   correction (DrERV) + threshold robustness.
4. **Mode-B mouse**: beta-dominant, matching IAP/MusD biology; species contrast vs bats.
5. **Literature cross-check** against named primary sources (see `04_FINDINGS.md`).

## Key numerical results
| Genome | mode | genus-resolved | dominant genera | lit. expectation | consistent? |
|---|---|---|---|---|---|
| Desmodus | A | 91.7% | Gamma 422, Beta 274 | Class I + II; beta present (DrERV) | yes |
| Antrozous | A | 93.7% | Gamma 767, Beta 486 | Class I + II | yes |
| Mus musculus | B | 99.9% | Beta 442, Gamma 325 | Class II (IAP/MusD) dominant | yes |
| Homo sapiens | B | 96.8% | Gamma 689, Beta 186 | Class I (HERV-H/W) most numerous | yes |
| Molossus | B | 98.2% | Gamma 378, Beta 144 | Class I + II (bat) | yes |

---
---

# Why these findings reliably represent true biological information

The reliability argument rests on **convergent independent evidence**, not a single run:

1. **Self-evidencing, not label-transfer.** Every call is made from the locus's *own*
   translated sequence searched against references, so it cannot merely echo the probe that
   detected it. The 23 DrERV loci where the independent call (Beta) *overrides* the probe
   label (Gamma) prove the classifier is reading the sequence, not the tag.

2. **Non-circular reference.** The classification reference is built from RefSeq independently
   of the study probe set, so a locus can be reassigned away from its detecting probe's genus
   (and is).

3. **Ground-truthed accuracy.** Leave-one-out gives a hard number — 98.6% correct genus —
   on sequences with known labels, and it holds across GAG/POL/ENV/PRO, supporting the
   probe-agnostic claim directly.

4. **Two independent methods agree.** Mode A (probe-detected valid loci, reclassified) and
   Mode B (reference-driven genome scan) both yield the same class-level shape on bats; they
   share only the taxonomy and reference, not the detection path.

5. **Species-specific signal (the strongest single line of evidence).** The method returns
   *different* answers for different species in the direction biology predicts: all three
   bats and human are gamma-lean (Class I), while mouse is beta-lean (Class II, IAP/MusD).
   The human↔mouse Class I/II **inversion** is reproduced from sequence alone, with no
   species label given. An artifact-producing method could not track a known species
   difference this specific.

6. **Literature concordance, verified.** The class-level composition (Class I + II dominant,
   Class III trace, alpha/delta/epsilon/lenti ≈absent as mammalian ERVs) matches named
   primary sources, and the DrERV reassignment matches a specific published discovery.

7. **Calibrated honesty.** Confidence + LCA backoff + abstention mean the method declines to
   over-call: ~6–8% of bat loci stay at subfamily, and it reports `unclassified` rather than
   forcing a genus — so the confident calls are the ones it can defend.

8. **Robustness.** Composition is stable across the weighting threshold; results are
   deterministic given the pinned reference.

### What is NOT claimed (honest boundaries)
- **Absolute abundance** — Mode B is a composition readout; the small exogenous-heavy
  reference undercounts highly-diverged and Class III (spuma-like: MuERV-L, HERV-L) elements
  (spuma is currently absent from the reference).
- **Exact proportions / minor genera** — the precise Gamma:Beta ratio and trace-genus counts
  depend on reference breadth and (Mode B) DIAMOND sensitivity; treat them as approximate.
- **Per-locus certainty for backed-off loci** — subfamily/family/unclassified calls are
  deliberately non-committal.

**Net:** the *direction and class-level composition* of each genome's ERV taxonomy — and the
specific probe-mislabel corrections — are reliable and literature-concordant. Finer
quantitative claims require a broader, spuma-inclusive reference (and, for the gold standard,
phylogenetic placement — the documented Phase-2 upgrade).
