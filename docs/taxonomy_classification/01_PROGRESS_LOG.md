# Progress Log - Taxonomy Classification

Append-only, newest at the bottom. Each entry dated.

---

## 2026-06-15 - P0 setup & environment recon

**Done**
- Confirmed branch `docs/taxonomy-classification-redesign`; all work to stay here.
- Inventoried tools/envs/data -> `03_ENVIRONMENT.md`. Key facts: DIAMOND + HMMER + BLAST present;
  placement tools absent; NCBI reachable; 5 model genomes staged in `testing-genomes`; 3 model bats
  already have local `valid_ranges`.
- Pivoted engine to **DIAMOND + weighted-LCA** (no new deps) with placement as a Phase-2 upgrade ->
  `02_DECISIONS.md` D1.
- Wrote `00_PLAN.md` (phases P0-P6) and `02_DECISIONS.md` (D1-D4).

**Open questions / blockers**
- Human + mouse have no local `valid_ranges`; deciding between running detection vs a known-reference
  validation shortcut (D3). Will resolve at P4.
- Curated reference accession list not yet assembled (P1, next).

**Next**
- P1: assemble a genus-labelled retroviral reference accession list (per ICTV genera + ERV classes),
  fetch via Entrez reusing `probe_extractor.py` patterns, build taxonomy map + DIAMOND DB, pin it.

---

## 2026-06-15 - P2/P3 iteration 1 (unweighted LCA) + literature verification

**Done**
- Wrote `workflow/scripts/taxonomy_lca.py` (curated ICTV retroviral taxonomy + pure
  `lca`/`rank_of` primitives + GFF3 evidence parser + summary CLI) and
  `tests/unit/test_taxonomy_lca.py` (9 tests, all pass).
- Ran on real valid tracks of all 3 model bats (reduced/merged) -> consistent
  **Gammaretrovirus(ClassI)-dominant, Betaretrovirus(ClassII)-second, rest trace**.
  Full (raw per-probe) Desmodus track -> uniform noise (diagnostic: unweighted LCA over
  singletons echoes each probe; degenerates without upstream aggregation).
- Verified literature (web): class-level shape CONFIRMED (gamma+beta dominant, spuma
  trace, alpha/delta/epsilon/lenti ~absent as mammalian ERVs). Found DrERV probe-label
  discrepancy (CSV=Gammaretrovirus vs literature Betaretrovirus) - a live mislabeling
  example; flags that exact gamma:beta ratio is probe-biased. See `04_FINDINGS.md`.

**Conclusion of iter 1.** LCA principle validated and class-level result is
literature-consistent, but unweighted LCA is reference/aggregation-bound. Need weighted
LCA vs an independent reference for trustworthy per-locus genus + proportions.

**Next (iter 2)**
- Build independent genus-balanced retroviral RefSeq reference (POL/GAG/ENV x 7 genera).
- Extract locus sequences (bedtools getfasta), diamond blastx, bitscore-weighted LCA.
- Re-run 3 bats; test whether gamma:beta balances and DrERV resolves to beta.

---

## 2026-06-15 - P3 iteration 2 (weighted LCA, independent reference) on bats

**Done**
- Built independent genus-balanced RefSeq reference (`taxonomy_reference_builder.py`):
  151 proteins, 6 genera (spuma=0: ICTV renamed spuma genera - known gap).
- Added `weighted_lca` (bitscore top-percent) to `taxonomy_lca.py` (+4 tests, 13 total pass).
- Wrote `taxonomy_classify_loci.py` (getfasta -> diamond blastx -> weighted LCA).
- **Desmodus** (RefSeq->GenBank seqname remap, same assembly): 91.7% genus, Gamma 422 /
  Beta 274; 23 loci probe-called Gamma are independently Beta (DrERV correction);
  threshold-robust. **Antrozous** (already GenBank-native): 93.7% genus, Gamma 767 /
  Beta 486. Two bats agree: Class I gamma dominant + Class II beta strong, rest trace.
- See `04_FINDINGS.md` iteration 2.

**Assembly-provenance blockers (the documented gotcha)**
- Molossus local track = older scaffold assembly (JACASF01); testing-genomes = newer
  chromosome assembly (CM138266). 0 seqname overlap -> cannot reuse; needs re-detection.
- Human + mouse: no local valid tracks at all -> need detection.

**Next (P4)**
- For Molossus/human/mouse: detect loci on the testing-genomes assemblies (tblastn of
  study probes vs prebuilt BLAST DB -> candidate loci), then classify. If tblastn on
  3 GB genomes is too slow this session, validate the engine on known human/mouse
  reference ERVs (HERV-K=Beta, MLV=Gamma, IAP=Beta, ERV-L=ClassIII) as a correctness
  ground-truth, and note genome-wide composition as pending compute.

---

## 2026-06-16 - P4 genome-wide (human/mouse/Molossus) + engine validation

**Done**
- tblastn whole-genome re-detection infeasible (>300 s timeout). Pivoted to DIAMOND
  windowed genome classification (`taxonomy_classify_genome.py`): tile->blastx->merge->LCA.
- **Mouse** (Mode B, more-sensitive, 126 min): 768 loci, Beta 442 > Gamma 325 -
  beta-dominant, matches IAP/ETn/MusD (Class II) dominance in mouse. Web-verified.
- **Leave-one-out engine accuracy: 98.6% correct genus** (145/147), high across ALL
  markers (GAG/PRO 100%, ENV/POL 98%) -> quantitative proof of probe-agnostic correctness.
- Human + Molossus genome runs launched (background, --sensitive). Human predicted
  Class I (gamma)-leaning (opposite of mouse) per HERV-H/W numerical dominance.

**Caveats recorded** (04_FINDINGS.md): Mode B is composition not abundance; reference is
exogenous-heavy + lacks spuma (Class III undercounted); mouse vs human use different
sensitivity (composition shape robust). Two validation modes kept distinct.

**Next**
- Collect human + Molossus results; finalise 04_FINDINGS; write three-register explainer
  + reliability justification (the requested deliverable).

---

## 2026-06-16 - P4-P6 complete: all 5 genomes classified, deliverables written

**Done**
- **Human** (Mode B): 938 loci, 96.8% genus, Gamma 689 >> Beta 186 (Class I-leaning) -
  matches HERV-H/W numerical dominance; inverse of mouse.
- **Molossus** (Mode B): 552 loci, 98.2% genus, Gamma 378 > Beta 144 - bat-like.
- All 5 genomes literature-consistent; human<->mouse Class I/II inversion reproduced from
  sequence alone (see 04_FINDINGS cross-species table).
- Code lint-clean (ruff) + formatted; 13 unit tests pass; LOO 98.6%.
- Result CSVs copied to `results/`; wrote `EXPLAINER_three_registers.md` (3 registers +
  reliability justification).

**Status: GOAL MET.** Engine developed (probe-agnostic, weighted-LCA), validated, and run on
all 5 model genomes with literature-concordant, self-verified findings; tracking + explainer
files in place.

**Known follow-ups (documented, not blocking):** add spuma (Class III) to reference;
broaden reference for absolute abundance; Phase-2 phylogenetic placement (needs new deps);
optional Snakemake integration.

---

## 2026-06-18 - EVEN FIELD: Mode A (valid-tier) for all 5; canonical data home established

**Key discovery:** pre-computed RetroSeek valid tracks for all 5 testing-genomes already exist at
`/mnt/v/workshop/testing-genomes/results/tracks/valid/` (GenBank namespace) - no ltrdigest rerun.

**Done**
- Established canonical dev-data home `data/taxonomy_dev/` (reference/ scans/ loci/ valid_tracks/
  + README); repointed all script defaults there; gitignored; memorized location + valid-track source.
- Ran **Mode A on all 5** (even field): Human 88.5% genus (Gamma 2765/Beta 1103); Mouse 94.0%
  (Beta 7841/Gamma 4171 - beta-dominant); Molossus 91.3% (Gamma 1627/Beta 792); Desmodus 91.9%
  (Gamma 460/Beta 302); Antrozous 94.3% (Gamma 1033/Beta 567). 3 bats + human gamma-dominant;
  mouse beta-dominant. Matches literature; reproduces human<->mouse Class I/II inversion.
- Results pinned in `data/taxonomy_dev/loci/<Genome>.modeA.csv`. Mode B deprecated (moved to
  `_deprecated_modeB/`). 04_FINDINGS updated with the definitive even-field table.

**Net:** all 5 model genomes on a uniform valid-tier Mode-A footing, one findable canonical location.

---

## 2026-06-18 - data-derived taxonomy + blastx default (dependency-clean)

- Taxonomy hierarchy now derived from NCBI (`taxonomy_build_hierarchy.py` -> `reference/taxonomy.tsv`);
  `taxonomy_lca.load_taxonomy()` loads it; runners call it. Built-in dict = fallback only.
  Behaviour-preserving (Desmodus identical). ERV class stays curated (not an NCBI rank). [D7]
- Switched Mode-A search default DIAMOND -> **blastx** (already in environment.yml; no new dep).
  Benchmark Desmodus: blastx 6.4s vs diamond 3.6s; blastx slightly more sensitive. `--search`
  flag keeps diamond optional. Regenerated all 5 canonical loci with blastx (composition
  unchanged: bats+human gamma-dominant, mouse beta-dominant; mouse genus-resolution up to 97.9%).
- 14 unit tests pass; ruff clean. Canonical data in data/taxonomy_dev/ all consistent.

---

## 2026-06-18 - TRIAL of redesign (pre-integration; gated on user report)

Env: created `retroseek-trial` conda env (mafft 7.526, iqtree 3.1.2, epa-ng 0.3.8, gappa 0.9.0).

**P1 Reference v2 (genus-comprehensive, data-derived):** `taxonomy_reference_builder.py` now keeps
ALL proteins per genus (gene best-effort incl OTHER; REX/TAX patterns) and uses the modern spuma
genera -> 383 proteins, 11 genera incl Class III. Data-derived `taxonomy.tsv` (Ortho+Spuma subfamilies)
+ `erv_class.tsv` data file with `load_erv_class()` loader. ERV class is the only curated piece (not
an NCBI rank). All lint clean, 14 LCA tests pass.

**P2 trees:** `taxonomy_build_tree.py` (mafft -> iqtree -> hmmbuild -> gappa taxon map + alignment-quality
flag). MFP ModelFinder was far too slow (~hr on POL); switched to fixed `LG+F+G4`. POL/GAG building.

**P3/P4 classifier rewrite** (`taxonomy_classify_loci.py`): LTR-element-anchored loci + gene-partition
+ per-gene dispatcher (placement POL/GAG else LCA; REX/TAX by presence) + mosaic + provenance.
`taxonomy_placement.py` (mafft --add -> epa-ng -> gappa) written.

**Two real bugs the trial caught & fixed (value of trialing):**
1. *ENV-bias*: locus summary picked highest-confidence gene, but LCA confidence is competition-
   dependent so ENV (few divergent refs) won spuriously (e.g. ENV->Alpharetrovirus@1.000). Fixed:
   choose locus call by marker reliability (POL>GAG>...>ENV), placement preferred; mosaic over main genes.
2. *LTR-anchoring*: workshop valid track is gene-merged with NO Parent attr; Parent-grouping fell back
   to 1-gene-per-locus (0 mosaics). Fixed: anchor loci by OVERLAP with `ltrdigest` LTR_retrotransposon
   elements (workshop, GenBank namespace, all 5 genomes).

**Desmodus LCA-only after fixes:** 355 LTR-anchored proviral loci, multi-gene (ENV,GAG,POL=92; GAG,POL=76...),
2 mosaics, Gammaretrovirus 194 > Beta 136 (Class I dominant - literature-consistent for vampire bat).

---

## 2026-06-18 - TRIAL COMPLETE (report delivered; awaiting user review)

All 5 genomes classified end-to-end (placement POL,GAG + LCA, LTR-anchored, mosaic-aware):
- bats + human Class I (gamma) dominant; mouse Class II (beta) dominant - human<->mouse inversion
  reproduced. Genus resolution 89-98%. Mosaics 3-296.
Verification: LOO 96.7% (all genes incl REX/TAX); determinism byte-identical; plot<->table concordance
PASS; literature-consistent. 5 design bugs caught & fixed (ENV-bias, LTR-anchoring, blastx qframe,
placement all-gap crash, placement over-backoff).
Outputs pinned in data/taxonomy_dev/ (loci/ summary/ plots/ tracks_igv/). Report:
docs/taxonomy_classification/TRIAL_REPORT.md. Compute-intensive detection outputs preserved.
**STOP - integration not started; awaiting user say-so.**

---

## 2026-06-18 - Refinement round R1 (user-approved: gene-agnostic + raxml-ng + revert to Parent)

- **R1.1 Parent-based loci (revert from overlap).** New `taxonomy_annotate_parent.py` writes
  `Parent=<LTR_retrotransposon>` into the valid tracks (overlap with ltrdigest elements, computed
  once); `build_loci` reverted to group by `Parent`. Readable tracks. Note: the Parent relationship
  is natively present in the pipeline's ltrdigest track (protein_match->element); only the gene-merged
  workshop *valid* track dropped it - production `ranges_analysis.R` should emit it natively.
  *Off-by-one fix:* overlap was missing `+1` for 1-based-inclusive coords, so features ABUTTING an
  element (gap=0) were dropped as orphans (5-12/genome); fixed -> **0 orphans** all 5 genomes.
- **R1.2 gene-agnostic.** Removed hard-coded GENE_PRIORITY/MAIN_GENES/DIAGNOSTIC_GENE_GENUS. Priority
  + mosaic set now from ordered `--main-probes`; diagnostic genes auto-detected from the reference
  (`auto_diagnostic`: a gene whose members are all one genus). No gene privileged in code.
- **R1.3 raxml-ng.** Installed raxml-ng 2.0.2 in retroseek-trial. `taxonomy_build_tree.py` runs
  `raxml-ng --evaluate` -> `<gene>.raxml.bestModel`+`.bestTree`; `taxonomy_placement.py` uses them
  (calibrated model -> less over-backoff). Trees rebuilding.
- Next: re-run all 5 on parented tracks; re-verify (incl. gene-agnosticism permutation check); update report.

## 2026-06-18 - R1 COMPLETE (refined trial; awaiting user review)

All 5 re-run with Parent-based loci + gene-agnostic ordering + raxml-optimized placement:
Desmodus 196:135, Antrozous 462:269, Molossus 597:296, Human 1376:431 (all Class I);
Mouse 1205:4479 (Class II). Inversion holds; composition stable vs pre-R1.
Verify: determinism byte-identical (PASS); gene-agnosticism config-driven (PASS, POL-first vs
ENV-first changes calls); Parent 0 orphans; concordance PASS; LOO 96.7%; lint clean; 14 tests pass.
raxml-ng increased placement contribution (human 172->269, mouse 3900->4558). GAG tree flagged
low-confidence (mean_pident 19.7%). TRIAL_REPORT.md updated. STOP - integration awaits user say-so.
