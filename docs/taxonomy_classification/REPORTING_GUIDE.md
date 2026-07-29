# What the new ERV taxonomic-classification step does - a plain-language guide

*A reporting/presentation companion. It explains, in simple terms, the change made to
RetroSeek, where it sits in the pipeline, which tools do what, and which figures to show.*

---

## 1. The problem we fixed (one sentence)

RetroSeek was **excellent at finding** endogenous retrovirus (ERV) loci in a genome, but
**poor at saying what each one is** - i.e. which retroviral **genus** it belongs to.

**Why it was poor.** The old method just copied the label of the single best-scoring *probe*
hit onto each locus. But the key retroviral gene (`pol`, the reverse transcriptase) is so
conserved that one ERV looks similar to *many* reference viruses at once. Picking the single
best hit from that cloud is an uncalibrated guess, and the probe set is a *detection* tool, not
a classification system - it can't name a genus it has no probe for.

## 2. What the new step does (one sentence)

For every ERV locus we already found, we now **re-identify it from its own sequence** against
an **independent, genus-comprehensive reference**, gene by gene, and assign a **calibrated genus
call** - with a confidence, a method, and a flag for *mosaic* loci (where different genes point
to different genera, i.e. recombinants).

Think of the old way as *"which of my fishing lures caught it?"* and the new way as
*"let's actually run the catch through a reference library and identify the species."*

---

## 3. Where it sits in the pipeline

The new step is **purely additive** - nothing upstream changed. It runs *after* the existing
detection, on the validated ERV loci:

```
genome -> LTR discovery (LTRharvest/LTRdigest) -> probe BLAST -> valid ERV loci
                                                                   |
                                              (existing pipeline --+)
                                                                   v
                                            * NEW: taxonomic classification *
                                                                   v
                            per-locus genus calls  ->  tables + IGV tracks + plots
```

Two new things had to exist for it to work:

1. **A reference** (built once, then reused) - the "library" we identify loci against.
2. **The classifier** (runs per genome) - does the identification.

---

## 4. Step by step - what happens and which tool does what

### Part A - Build the reference (run once: `make reference`)

| Step | Tool | What it's for (plain terms) |
|------|------|------------------------------|
| Download reference proteins | **NCBI Entrez** (Biopython) | Pull a balanced set of known retroviral proteins, one batch per genus, straight from NCBI - so we identify loci against viruses that are *independent* of our probes (no circular reasoning). |
| Build the family tree of relationships | **NCBI Taxonomy** | Work out the genus -> subfamily -> family hierarchy automatically, so the method can "back off" to a higher rank honestly when it can't pin a genus. |
| Add the ERV class labels | curated `erv_class.tsv` | The one hand-made file: maps each genus to ERV **Class I / II / III** (the standard biological grouping). |
| Align the reference proteins | **MAFFT** | Line up the reference proteins of a gene (e.g. POL) so they're comparable position-by-position. |
| Build a reference tree | **IQ-TREE** then **raxml-ng** | Build the evolutionary tree of those reference proteins - the "map" onto which we'll place each new locus. raxml-ng fine-tunes it for accurate placement. |
| Make a search profile | **HMMER** | A statistical profile of the gene, used to slot query sequences into the reference alignment. |

**Reproducibility:** the tree-building is seeded from the global `seed` in the config, and a
**provenance manifest** records exactly which proteins/genera went in (counts + a content hash),
so any result can be traced to the exact reference snapshot.

### Part B - Classify each genome's loci (`RetroSeek --classify`)

| Step | Tool | What it's for (plain terms) |
|------|------|------------------------------|
| Group hits into loci | (RetroSeek) | Collect the per-gene hits belonging to the same ERV element (its LTR "`Parent`"), so each locus is one provirus with its genes. |
| Cut out each gene's sequence | **bedtools** | Extract the actual DNA of each gene region from the genome. |
| Match against the reference | **BLAST (blastx)** | Compare each locus gene to all reference proteins - the raw evidence of "what it looks like". |
| Decide the genus - main method | **weighted-LCA** (RetroSeek) | Instead of one best hit, take the *cluster* of top hits and find their lowest common ancestor on the taxonomy. A clean single-genus cluster -> confident genus; a mixed cluster -> an honest higher rank. |
| Decide the genus - precision method | **EPA-ng + gappa** | For the reliable gene (POL), *place* the locus onto the reference tree and read off the genus from where it lands. Used when it resolves cleanly; otherwise we fall back to weighted-LCA. |
| Combine + flag mosaics | (RetroSeek) | Merge the per-gene calls into one locus call; if genes disagree, flag it as a **mosaic** (recombinant) and record the composition. |

The output of Part B **is** the genus-founded "ERV assembly": each locus is one provirus with
its consensus genus, confidence, mosaic status, ERV class, and structural completeness.

---

## 5. Where to find the outputs (for figures + tables)

All paths are under the run's results folder
(`results/` locally, or `/mnt/v/workshop/testing-genomes/results/` for the 5-genome run):

| Output | Location | Use it for |
|--------|----------|------------|
| **Plots** (5 PNGs) | `results/plots/taxonomy/` | The figures for your talk - see Section 6. |
| **Per-locus table** (the headline data) | `results/tables/taxonomy_classification/<genome>.loci.csv` | Each locus: `genus_call`, `rank`, `confidence`, `method`, `is_mosaic`, `mosaic_composition`, `erv_class`, per-gene calls. Open in Excel/R. |
| **Genome browser track** | `results/tracks/taxonomy/<genome>.gff3` (+ `.bed`) | Load in **IGV** to show genus calls in genome context, colour-by-genus. |
| **Reference provenance** | `data/taxonomy_reference/manifest.yaml` | Cite which reference snapshot was used (genera + counts). |

---

## 6. The five figures and what each one shows

In `results/plots/taxonomy/`:

| File | What it shows | The point to make |
|------|---------------|-------------------|
| `erv_class_composition.png` | Class I/II/III mix per species | **The headline.** Bats + human are Class I-dominant; mouse is Class II-dominant - the known human<->mouse inversion, recovered from sequence. |
| `genus_composition.png` | Genus mix per species (Gamma, Beta, ...) | The finer-grained version: bats Gammaretrovirus-rich, mouse Betaretrovirus-rich (IAP/MusD biology). |
| `rank_resolution.png` | How deep each call goes (genus / subfamily / family) | Shows the method is *honest* - it backs off rather than forcing a wrong genus. |
| `method_mix.png` | How calls were made (placement vs weighted-LCA) | Shows both engines contributing; placement strongest on bats/mouse. |
| `mosaic_alluvial.png` | Gene->genus flows within mosaic loci | Shows we can detect **recombinant** ERVs (different genes, different ancestry). |

> **Suggested single slide:** `erv_class_composition.png` next to a one-line "before vs after"
> (old = best-probe label; new = calibrated genus call from an independent reference).

---

## 7. Headline results (numbers you can quote)

Verified on the 5 model genomes (3 bats, human, mouse), using the freshly built reference:

- **Bat (Desmodus rotundus):** 355 loci, **Class I dominant** (Gammaretrovirus 196 vs Betaretrovirus 135).
- **Mouse (Mus musculus):** 6436 loci, **Class II dominant** (Betaretrovirus 4480 vs Gammaretrovirus 1217).
- -> The **human<->mouse Class I/II inversion** is reproduced - a known, published biological fact,
  recovered here purely from sequence (a strong sanity check that the method is right).
- **Accuracy:** ~**96.7%** genus accuracy in a leave-one-out test on the reference.
- **Reproducible:** re-running gives **byte-identical** results (seeded trees + pinned reference).
- **Mosaic ERVs detected:** e.g. ~150 in mouse - recombinants the old method couldn't surface.

---

## 8. Before vs after (the slide-ready summary)

| | **Before** | **After** |
|---|-----------|-----------|
| Basis of the call | Best-scoring *probe* label | The locus's *own* sequence vs an independent reference |
| Calibrated? | No (single argmax over a cross-genus cloud) | Yes (cluster-based LCA + phylogenetic placement) |
| Can it name genera with no probe? | No | Yes (reference covers all retroviral genera) |
| Honest about uncertainty? | No | Yes (backs off to subfamily/family) |
| Detects recombinant (mosaic) ERVs? | No | Yes |
| Reproducible / traceable? | - | Seeded + provenance manifest |

---

## 9. How to regenerate everything (commands)

```bash
conda activate RetroSeek
make reference                                   # build the reference once (network)
./RetroSeek --classify --configfile data/config/config.local.yaml --cores all
# -> tables in results/tables/taxonomy_classification/
# -> tracks in results/tracks/taxonomy/
# -> figures in results/plots/taxonomy/
```

*Design & rationale: [ADR-007](../adr/ADR-007-taxonomic-classification.md). Full method detail:
the other files in this folder and [`docs/configuration.md` -> `classification`](../configuration.md).*
