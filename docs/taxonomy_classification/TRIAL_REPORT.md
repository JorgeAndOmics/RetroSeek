# Trial Report - ERV taxonomic-classification redesign on the 5 model genomes

**Status: trial complete (incl. refinement round R1). Awaiting your review before any integration.**

**Refinement round R1 applied (user-approved):** (1) **gene-agnostic** - removed the hard-coded
gene priority / mosaic-set / diagnostic-gene literals; gene ordering now from `--main-probes`,
diagnostic genes auto-detected from the reference (verified config-driven). (2) **raxml-ng-optimized
placement** - `raxml-ng --evaluate` fits the model per gene; EPA-ng uses the `.bestModel`/`.bestTree`
(placement contribution rose, e.g. human 172->269, mouse 3900->4558). (3) **Parent-based loci** -
valid tracks annotated with `Parent=<LTR_retrotransposon>` (readable; 0 orphans after a 1-based-
inclusive off-by-one fix); classifier groups by `Parent`. Composition unchanged (robust).

This reports a pre-integration trial of the redesigned ERV taxonomic classifier (per-gene,
mosaic-aware, placement + LCA, data-derived, no hard-coding), run end-to-end on the 5 model
genomes, with verification on every axis. All artifacts are pinned under `data/taxonomy_dev/`.

## What was run

Pipeline (per genome): preserved `valid` tracks -> **LTR-element-anchored loci** (overlap with
`ltrdigest` `LTR_retrotransposon` elements) -> **gene-partitioned** evidence -> per-gene
**dispatcher**: POL/GAG **phylogenetic placement** (MAFFT->EPA-ng->gappa) refining to genus when it
resolves, **weighted-LCA** otherwise, REX/TAX by **presence** -> mosaic-aware locus summary.

- Reference v2: **genus-comprehensive**, 383 RefSeq proteins, 11 genera incl. the modern spuma
  genera (Class III) and accessory genes (REX/TAX); data-derived `taxonomy.tsv`; `erv_class.tsv` data
  file. Only the ERV-class map is curated (not an NCBI rank); nothing else is hard-coded.
- POL + GAG reference trees (IQ-TREE `LG+F+G4`, 1000 UFBoot), pinned.
- Tools: `retroseek-trial` conda env (mafft, iqtree, epa-ng, gappa); `blastx` (default search, no new
  dep); the committed `environment.yml` was NOT modified.
- Compute-intensive detection outputs (suffixerator/ltrharvest/ltrdigest/valid) were **preserved,
  never deleted** (`/mnt/v/workshop/testing-genomes/results/`).

## Results - all 5 genomes

| Genome | loci | genus % | mosaic | dominant class | Gamma : Beta | method place/lca |
|---|---|---|---|---|---|---|
| Desmodus_rotundus (bat) | 355 | 94% | 3 | **Class I** | 196 : 135 | 210 / 123 |
| Antrozous_pallidus (bat) | 757 | 98% | 20 | **Class I** | 462 : 269 | 466 / 275 |
| Molossus_molossus (bat) | 928 | 97% | 26 | **Class I** | 597 : 296 | 726 / 170 |
| Homo_sapiens | 2035 | 90% | 34 | **Class I** | 1376 : 431 | 269 / 1572 |
| Mus_musculus | 6436 | 89% | 171 | **Class II** | 1205 : **4479** | 4558 / 1149 |

(R1 numbers - Parent-based loci, gene-agnostic ordering, raxml-optimized placement.)

Per-genome detail in `data/taxonomy_dev/loci/<genome>.modeA.csv`; summary in
`data/taxonomy_dev/summary/per_genome_summary.csv`.

## Verification (all axes)

| Axis | Method | Result |
|---|---|---|
| **Engine accuracy** | leave-one-out genus on Reference v2 (self-hit excluded) | **96.7%** (350/362); per gene ENV 95%, GAG 97%, POL 98%, PRO/REX/TAX 100%, OTHER 96% |
| **Literature consistency** | composition vs published ERV biology | bats+human **Class I (gamma) dominant**; mouse **Class II (beta) dominant** - the human<->mouse inversion reproduced; spuma/Class III trace; alpha/delta/epsilon/lenti trace - all as expected |
| **Determinism / end-to-end** | regenerate classification layer from preserved detection; re-run + diff | **byte-identical** (PASS) |
| **Plot<->table concordance** | plotted genus sums vs summary counts, all 5 | **PASS** (exact) |
| **Probe-agnosticism** | REX/TAX (non-POL genes) classify | REX/TAX 100% in LOO; present in reference -> REX-only loci now resolvable |
| **Gene-agnosticism (R1)** | POL-first vs ENV-first `--main-probes` (placement off) | locus calls change with config (PASS) - ordering data-driven, no gene hard-coded |
| **Parent integrity (R1)** | annotate valid features -> ltrdigest elements | 0 orphans across all 5 (valid  subset-or-equal  elements, as expected) |

Literature anchors (web-verified earlier, see `04_FINDINGS.md`): mouse high-copy ERVs are IAP/ETn/MusD
(Class II, betaretrovirus-like) - Maksakova 2006 PLoS Genet, Zhang 2012 Genome Biol; human's most
numerous HERVs are Class I (HERV-H/W) - Belshaw 2005 MBE; bat ERVs Class I gamma-dominant - Hayward
2013 PNAS / 2015 Sci Rep; DrERV vampire-bat betaretrovirus - J Virol 2015.

## Outputs produced (and checked)

- **Tables:** `loci/<g>.modeA.csv` (per-locus: genus_call, rank, confidence, method, per_gene,
  is_mosaic, mosaic_composition, erv_class, detection provenance, ref_version); `summary/`.
- **IGV tracks:** `tracks_igv/<g>.classified.gff3` (genome-coordinate, genus/rank/method/mosaic in
  attributes) - loadable in IGV over the genome.
- **Plots:** `plots/genus_composition.png`, `rank_resolution.png`, `method_mix.png` (built from the
  same `plot_long.csv` the summary derives from -> concordant by construction, and verified).

## Issues caught and fixed during iteration (value of trialing)

1. **ENV-confidence bias** - LCA confidence is competition-dependent, so ENV (few refs) won locus
   calls spuriously. Fixed: locus call chosen by marker reliability (POL>GAG>...>ENV), placement
   preferred; mosaic over main genes only.
2. **LTR-anchoring** - the namespace-matching valid track lacks `Parent`; fixed by anchoring on
   `ltrdigest` element overlap (enabled multi-gene loci + mosaics).
3. **blastx frame** - used `sframe` (always 0 for blastx) -> garbage translations -> empty placement;
   fixed to `qframe`.
4. **Placement crash** - all-gap query after `mafft --add` aborted EPA-ng (killed Antrozous/Mouse);
   fixed with a gap-guard (such queries fall back to LCA).
5. **Placement over-backoff on divergent HERVs** - placement-as-override sent 89% of human loci to
   subfamily; fixed so placement refines **only when it resolves a genus**, else LCA - recovering
   human to 90% genus while keeping placement precision on bats/mouse.

## Honest caveats / open issues (for the integration discussion)

- ~~Placement model unoptimized~~ - **RESOLVED in R1** via `raxml-ng --evaluate` (per-gene
  `.bestModel`); placement contribution rose (human 172->269, mouse 3900->4558). raxml-ng is now a
  required tool for the placement path.
- **GAG tree is low-confidence** - the builder flags GAG alignment quality LOW (mean pairwise
  identity 19.7%, 72% gaps); GAG cross-genus placement should be treated cautiously (POL is the
  reliable tree). ENV is intentionally not a placement gene.
- **Placement contribution varies by divergence** - strong on bats/mouse, weaker on human (LCA still
  carries most human calls). Expected and handled by the hybrid; placement's value is genome-dependent.
- **ENV remains the unreliable marker** (kept in LCA; not a placement gene by default).
- **Reference breadth** - 383 proteins is modest; broader curation would help rare genera and
  sub-genus resolution. POL/GAG trees use a fixed model and modest taxon sampling.
- **Mode B (genome-wide DIAMOND scan) was abandoned** - everything is Mode A (valid-tier).

## Verdict

The redesign runs end-to-end on all 5 genomes, is **reproducible (byte-identical), literature-
consistent (incl. the human<->mouse inversion), 96.7% accurate (leave-one-out), probe-agnostic, and
hard-coding-free**, with mosaic detection and IGV-ready outputs whose plots match the tables. The
trial also surfaced and fixed 5 real design bugs. Recommended next step *if you approve*: production
integration (Snakemake rules, `environment.yml` incl. mafft/iqtree/epa-ng/gappa/raxml-ng, config
surface, docs/schema sync, tests) - **not started, awaiting your say-so.**

## How to reproduce / inspect (production path)

The trial has been **integrated into the pipeline** (see [ADR-007](../adr/ADR-007-taxonomic-classification.md)).
Reproduce via the production rules in the `RetroSeek` env (the standalone trial scripts -
`taxonomy_trial_report.py`, `taxonomy_trial_plots.R`, `taxonomy_classify_genome.py`,
`taxonomy_annotate_parent.py` - were retired; their function is now the Snakemake stage + native
`Parent` in `validation.R`):

```bash
conda activate RetroSeek
# 1. Build the reference once (Entrez fetch + POL placement tree):
make reference                       # or: ./RetroSeek --build-reference
# 2. Classify every valid locus + render the taxonomy plot panel:
./RetroSeek --classify --cores all --skip-validation
```

Outputs: `results/tables/taxonomy_classification/<G>.loci.csv` (per-locus genus calls - the
genus-founded ERV assembly), `results/tracks/taxonomy/<G>.gff3` + `.bed` (IGV), and
`results/plots/taxonomy/` (5 PNGs). To classify a single genome ad-hoc, call
`taxonomy_classify_loci.py <valid.gff3> <genome.fa> --ref-dir data/taxonomy_reference
--placement-genes POL --out <out.csv>` directly. The original trial scratch lived under
`data/taxonomy_dev/` (gitignored); the production reference is `data/taxonomy_reference/`.
