# Taxonomy Classification — Working Plan

**Branch:** `docs/taxonomy-classification-redesign` (all work stays here).
**Goal:** Replace best-bitscore probe-label transfer with a principled, self-evidencing,
probe-agnostic taxonomic classifier for ERV loci; iterate on the 5 model genomes until findings
are literature-consistent and self-verified.

Companion files:
- `01_PROGRESS_LOG.md` — dated running log (append-only).
- `02_DECISIONS.md` — key technical decisions + rationale.
- `03_ENVIRONMENT.md` — environment / tool / data inventory (ground truth of what's runnable).
- `04_FINDINGS.md` — empirical results vs literature (created when first results land).
- design source: `../taxonomy_classification_redesign.md` (the three-register theory write-up).

---

## The 5 model genomes (in `/mnt/v/databases/testing-genomes/`)

| Genome | Group | Role | Local `valid_ranges`? |
|---|---|---|---|
| `Homo_sapiens` | mammal (primate) | literature gold standard (HERV-K/W/L, ERV1/2/3) | no — needs detection |
| `Mus_musculus` | mammal (rodent) | literature gold standard (IAP, MuERV-L, MLV/ETn) | no — needs detection |
| `Desmodus_rotundus` | bat (Phyllostomidae) | bat truth anchor (DrERV = Gammaretrovirus) | **yes** |
| `Antrozous_pallidus` | bat (Vespertilionidae) | bat generalisation | **yes** |
| `Molossus_molossus` | bat (Molossidae) | bat generalisation | **yes** |

Data boundary: `testing-genomes` is the sanctioned dataset. The real 102/103-bat study data and the
`*.local.yaml` real catalog remain off-limits — not read.

---

## Engine (see 02_DECISIONS.md D1)

Translated-marker homology search + **weighted-LCA** taxonomic assignment (MEGAN/detectEVE paradigm),
using tools already installed (DIAMOND, HMMER, BLAST). Phylogenetic placement is the documented
Phase-2 upgrade, gated on approval to add `mafft`/`epa-ng`/`gappa`/`iqtree`.

Core flow (per valid hit = per-gene marker):
1. Extract the locus's own genomic interval (from `valid_ranges.gff3` coords + genome FASTA).
2. `diamond blastx` it against a genus-labelled retroviral protein reference DB.
3. Weighted-LCA over hits above a score threshold → `genus | higher-rank | unclassified` + confidence.
4. Per-gene calls → locus summary + mosaic flag (recombination made explicit).

This satisfies the design invariants: probe-agnostic (G1 — keys off the locus sequence, not the
probe), self-evidencing (G2), calibrated + abstaining (G3 — LCA backoff), mosaic-aware (G4),
deterministic (G5 — pinned reference, fixed thresholds), additive (G6 — consumes existing loci).

---

## Phases  (ALL COMPLETE — see 01_PROGRESS_LOG.md / 04_FINDINGS.md / EXPLAINER_three_registers.md)

- **P0 — Setup & tracking.** ✅ Files in place; environment + dataset inventoried.
- **P1 — Reference package.** ✅ Independent genus-balanced RefSeq reference (151 proteins,
  `data/taxonomy_reference/`); spuma gap noted.
- **P2 — Classifier prototype.** ✅ `taxonomy_lca.py` (+`weighted_lca`), `taxonomy_classify_loci.py`,
  `taxonomy_classify_genome.py`, `taxonomy_reference_builder.py`, `tblastn_to_loci.py`; 13 unit tests.
- **P3 — 3 bats.** ✅ Desmodus + Antrozous (Mode A, valid loci) + Molossus (Mode B). Gamma-dominant.
- **P4 — Human/mouse.** ✅ Mode-B genome composition. Human Class I-leaning, mouse Class II-leaning.
- **P5 — Outputs.** ✅ Per-locus/per-genome CSVs in `results/`; cross-species summary in 04_FINDINGS.
- **P6 — Verification & write-up.** ✅ Literature verified (named sources); LOO 98.6%; reliability
  argument + three-register explainer written.

Outcome: all 5 model genomes classified with literature-concordant, self-verified findings.
Follow-ups (non-blocking): add spuma/Class III to reference; broaden for abundance; Phase-2
phylogenetic placement (needs `mafft`/`epa-ng`/`gappa`/`iqtree` — a new-dependency decision);
optional Snakemake integration.

## Definition of done

Findings on the 5 genomes that match published mammalian ERV genus composition (verified against
named sources), with a written justification of why they reliably represent true biology, plus the
three-register explainer — all captured in these text files.
