# Architecture

RetroSeek is a Snakemake-orchestrated bioinformatics pipeline combining Python, R, and Bash to detect endogenous retroviral (ERV) integrations in eukaryotic genomes.

## Stack

- **Orchestrator**: Snakemake ≥ 8.
- **Python 3.10**: probe extraction, BLAST execution, object serialisation, validation, CLI.
- **R 4.3**: range analysis, segmentation, hotspot detection, pair detection, plotting.
- **External tools** (provided by the conda env):
  - BLAST+ (`tblastn`, `blastn`, `makeblastdb`).
  - GenomeTools (`gt suffixerator`, `gt ltrharvest`, `gt ltrdigest`, `gt gff3`).
  - NCBI Datasets CLI (`datasets`).
  - Pfam HMMs (downloaded on demand).
- **Reproducibility**: single conda/mamba env at [`data/config/environment.yml`](../data/config/environment.yml).

## Data flow

```
probe CSV ──► probe_extractor (Entrez) ──► probe_dict.pkl
                                              │
genome list ──► genome_downloader (datasets) ──► {genome}.fa
                       │
                       ├─► blast_db_generator (makeblastdb)  ─► per-genome BLAST DB
                       ├─► ltr_index_generator (suffixerator)
                       │     └─► ltr_harvester                ─► {genome}.gff3 + fa
                       │           └─► ltr_digester (Pfam)    ─► {genome}.gff3 (annotated)
                       │
                       └─► full_genome_blaster (tblastn)      ─► {genome}.pkl
                                                                  │
                             obj2dict ◄─── all {genome}.pkl ──────┘
                                   │
                                   ├─► species_segmenter.R
                                   └─► ranges_analysis.R
                                             │
                                             ├─► GFF3 tracks (original, candidate, valid [+valid _reduced])
                                             ├─► overlap matrix (CSV)
                                             └─► plot dataframes (Parquet)
                                                   │
                                                   ├─► plot2sort.R          (global plots)
                                                   ├─► circle_plot_generator.R
                                                   ├─► hotspot_detector.R
                                                   └─► pair_detector.R
```

## Rule overview

All rules follow `<name>_setup` (per-wildcard) + `<name>` (aggregate via `expand`). Groups:

**Acquisition** — `genome_downloader`, `pfam_hmm_downloader`.
**Indexing** — `blast_db_generator`, `ltr_index_generator`.
**LTR discovery** — `ltr_harvester`, `ltr_digester` (depends on both LTR_harvest outputs and the Pfam HMM download).
**BLAST search** — `probe_extractor`, `full_genome_blaster`, `blast_pkl2parquet`.
**Integration & segmentation** — `species_segmenter`, `ranges_analysis` (phase-modular: `workflow/scripts/ranges_analysis.R` is a thin orchestrator over sibling modules in `workflow/scripts/range_analysis/{io,granges_build,filtering,reductions,validation,erv_assembly,plot_dataframe,exporters}.R`).
**Solo-LTR detection** — `ltr_retriever_prefilter`, `ltr_retriever`, `solo_ltr_integrator`, `solo_ltr_detector` (aggregate). Runs LTR_retriever over LTRharvest output pre-filtered by `valid_ranges.gff3`; propagates RetroSeek probe labels onto discovered solo LTRs. See [`docs/solo_ltr.md`](solo_ltr.md) for the full mechanism and [ADR-003](adr/ADR-003-ltr-retriever-pre-filter.md) for the pre-filter rationale.
**Taxonomic classification** — `taxonomy_reference` + `taxonomy_reference_trees` (build-once reference: taxon-comprehensive proteins via Entrez, NCBI taxonomy, curated `erv_class.tsv`, per-gene placement trees; `make reference` / `--build-reference`), then `taxonomy_classify` (per-locus taxon calls), `taxonomy_orphans` (recovered, proximity-clustered orphans tier), `taxonomy_plot_generator`, and `loss_analysis` (`--classify`). The classification is **rank-agnostic** (ADR-008): the axis of taxa a locus can resolve to is declared by `classification.reference_taxa` (any rank; defaults to the distinct probeset `Label` values). Reclassifies each valid LTR-element locus from its own marker sequence — POL/GAG by phylogenetic placement, weighted-LCA otherwise — into a calibrated `taxon_call` + `rank` (resolved to an axis taxon, or an honest backoff) with confidence, `confidence_tag` HC/LC vs `classification.confidence_min`, `n_blastx_hits`, method, mosaic, and ERV class. The per-locus loci table is the taxon-founded ERV assembly, carrying a discrete `structure_class` (full / partial / gene) and per-provirus `domain_tier` (domain_selected / domain_unlisted / non_domain) alongside the taxon call (ADR-009); legacy probe `virus`/`label` are kept as detection provenance (`probe_label_set`). Note the **valid tier is now the whole LTR-flanked set** — every LTR-overlapping hit is labelled (with `domain_tier` + a per-hit `domain_hit_class`) rather than filtered, so no LTR-flanked hit is discarded (ADR-009). A parallel **orphan tier** recovers non-LTR-associated hits — reduced BLAST hits overlapping no retrotransposon — by running the same classifier on them and keeping only those that earn a taxonomic call (the novel-retrovirus path; `source=orphan`). Orphans have no LTR element to group by, so their hits are clustered by **physical overlap** (co-location = same feature) into single non-overlapping loci — assembled by the same machinery, flagged `source=orphan` (proximity-*inferred* vs LTR-*confirmed*), and capped at the widest LTRdigest element per genome (clusters wider than any real provirus get an `oversized` flag, kept for filtering). Overlap-only means orphan loci are mostly single-gene (adjacent genes don't overlap) — a conservative dedup, not speculative multi-gene assembly (ADR-010). The LTR-flanked ∪ orphan union is written as one fully non-overlapping authoritative `catalog.csv` (LTR-flanked-precedence resolves the cross-tier edge). `loss_analysis` then unions the ranges-analysis and blastx-stage counts into one per-stage attrition funnel and exports per-genome novel candidates (loci with `n_blastx_hits == 0`). See [configuration.md `## classification`](configuration.md), [`docs/taxonomy_classification/`](taxonomy_classification/), [ADR-007](adr/ADR-007-taxonomic-classification.md), and [ADR-008](adr/ADR-008-rank-agnostic-classification.md).
**Downstream analyses** — `plot_generator` (phase-modular: `workflow/scripts/plot2sort.R` is a slim orchestrator over sibling modules in `workflow/scripts/plot2sort/{helpers,io,plots_distribution,plots_categorical,plots_sankey}.R`; emits 22 PNGs with auto-scaled canvases driven by `auto_dims()`, opt-in long-tail collapse via `sankey_top_n`, and waffle auto-scaling, including a detected-virus stacked bar per species), `stage_plot_generator` (18 middle-stage PNGs: concordance/structure/funnel/multiplicity + pre-reduction overlap + LTR-interaction), `erv_like_plot_generator` (7 *structural* PNGs read from the taxon-founded taxonomy loci table: completeness, canonical order, gene combinations, length, main-gene count, taxon×gene composition, structure_class), `circle_plot_generator`, `hotspot_detector`, `pair_detector`. `--generate-global-plots` drives the ranges panel (`plot_generator` + `stage_plot_generator`) plus the structure panel (`erv_like_plot_generator`).

## Outputs

- `results/tracks/` — GFF3 per stage (`original`, `candidate`, `orphans`, `valid`, `flanking_ltr`, `ltrdigest`, `ltrharvest`, `solo_ltr`, `ltr_retriever/`, `taxonomy`). `orphans/{genome}.gff3` is the orphan (non-LTR) hit set, carrying a synthetic `Parent=` from proximity clustering; `orphans/{genome}.classified.{gff3,bed}` is the recovered, taxonomically-called subset. Only `valid` carries a `_reduced` variant (merged overlapping ranges); `original`/`candidate` are unreduced. The composite ERV "assembly" is no longer a separate `erv_like` track — it is now the taxon-founded `taxonomy/{genome}.gff3` (one feature per LTR-element locus, labelled by `taxon_call`). `solo_ltr/{genome}.gff3` carries probe_labels propagated from valid ERVs (see `docs/solo_ltr.md`).
- `results/tables/<name>/` — user-facing **CSV** copy of every table group, each in its own subdirectory (no loose files): `ranges_analysis` (`{genome}.{final_loci,homology_loci,ltr_structure,reduction_multiplicity,counts,provirus_overlap,ltr_interaction,probe_domain_overlap,reduction_coverage}.csv`), `overlap_matrix`, `segmented_species`, `probe_pairs`, `solo_intact_ratio`, `hotspots`, `probe_dict`, `full_genome_blast`, `taxonomy_classification`.
- `data/tables/<name>/` — the pipeline-internal **Parquet** copy of each of those table groups, mirroring the `results/tables/` layout. `data/tables/_input/` holds the user-provided probe CSV.
- `results/manifest/` — per-genome run manifest YAMLs: provenance only (generator build, timestamp, input md5s, resolved parameters, seed).
- `results/plots/` is laid out to mirror the pipeline stages (so the folder tree *is* documentation):
  - `plots/ranges/` — the **ranges stage** (pre-classification): `homology/` (22 PNGs from `plot2sort` — density, raincloud, query-coverage, bar, **virus-by-species stacked bar**, balloon, heatmap-probe×species, virus waffle, three Sankey variants, split into `main`/`accessory`/`full` where applicable) and `integration/` (18 middle-stage PNGs from `stage_plot_generator` — homology↔LTR concordance, structure, refinement funnel, multiplicity, **pre-reduction overlap/redundancy** — self-overlap degree, reciprocal-overlap fraction, reduction fold, total-bp before/after — and **LTRharvest/LTRdigest interaction** — distance-to-nearest-retrotransposon, position-within-element metagene, probe×Pfam-domain overlap, per-feature breakdown, strand concordance, element length vs hits). Every plot stamps its range tier in the subtitle.
  - `plots/classification/` — the **taxonomy stage**: `taxonomy/` (18 PNGs, below), `structure/` (7 structural PNGs from `erv_like_plot_generator` — completeness, canonical_order, gene_combinations, length_distribution, n_main_genes, composition_heatmap, structure_class; renamed from the `erv-like/` panel), and `loss/` (5 attrition-funnel PNGs, below).
  - `plots/circle/` (per-genome Circos-style PNG + PDF) and `plots/hotspot/` (enrichment PDFs) sit alongside.
- **Taxonomic classification** — `data/taxonomy_reference/` holds the build-once reference (`retro_reference.{faa,csv}`, `taxonomy.tsv`, `erv_class.tsv`, `trees/`, `manifest.yaml`); gitignored, rebuilt via `make reference`. Per-genome outputs: `results/tables/taxonomy_classification/{genome}.loci.csv` (+ Parquet under `data/tables/`) — the per-locus taxon-founded ERV assembly (`taxon_call`, `rank`, `confidence`, `confidence_tag`, `n_blastx_hits`, `method`, `is_mosaic`, `mosaic_composition`, `erv_class`, `completeness`, `canonical_order`, `structure_class`, `domain_tier`, `oversized`, per-gene calls, `probe_label_set` provenance, `source`); `{genome}.orphans.csv` — the recovered, proximity-clustered orphan tier (same schema, `source=orphan`, gated to classified-only); `{genome}.classification_counts.csv` + `{genome}.orphans_counts.csv` — blastx-stage loss counters (incl. `structure_*` / `loci_domain_*` tallies + `orphans_total`/`orphans_recovered`); `results/tables/taxonomy_classification/classification_report.csv` — tidy counts by taxon / confidence / method + mosaic + integration totals, split by tier; **`catalog.csv`** — the unified authoritative catalog: LTR-flanked proviruses ∪ clustered orphan loci as one fully non-overlapping record set (LTR-flanked-precedence), `source` keeping the confidence gradient (ADR-010); `results/tracks/taxonomy/{genome}.gff3` + `.bed` (IGV, colour-by-taxon; GFF3 attrs carry `structure_class` + `domain_tier`); `results/plots/classification/taxonomy/` — 20 PNGs: the composition panel (taxon composition, rank resolution, method mix, ERV-class composition, confidence + confidence_count + confidence_gradient [stacked count bars, HC/LC and a viridis confidence gradient]) + a **mosaic/recombination sub-panel** (mosaic_alluvial, mosaic_burden, mosaic_taxon_pairs recombination-partner heatmap, mosaic_gene_discordance, mosaic_composition_by_species) + the evidence/structure panel (evidence_depth, confidence_density, confidence_vs_evidence, structure_by_tier, source_yield, taxon_by_source, domain_tier_composition, structure_class_composition), derived from the loci + orphan tables so they stay concordant. **Loss analysis** — `results/tables/loss_analysis/loss_analysis.csv` (per-genome per-stage attrition funnel) + `novel_candidates/{genome}.novel_candidates.csv` (loci with `n_blastx_hits == 0`) + `results/plots/classification/loss/` — 5 PNGs (loss_funnel, step_retention heatmap, orphan_recovery yield, novel_burden, loss_waterfall).
- `results/tables/hotspots/` + `results/tracks/hotspots/` — per-window CSV/Parquet + provenance manifest, merged-region GFF3/BED, and `results/plots/hotspot/` (Manhattan / karyotype / Q-Q / summary). Detection is a deterministic Negative-Binomial GLM (`hotspot_detector.R` + `hotspot_analysis/*`); see [configuration.md `## hotspot`](configuration.md).
- `data/ltr_scn/` — LTRharvest screen-format intermediates: `{genome}.scn` (raw LTRharvest output) plus `{genome}_retroviral.scn` (Coupling-A filtered) and `{genome}_full.scn` (byte-equal passthrough). The prefilter rule emits both filtered files from one read pass; one of them feeds LTR_retriever.

## Key configuration

All user-tunable parameters live in [`data/config/config.yaml`](../data/config/config.yaml); validation rules are in [`data/config/schema.yaml`](../data/config/schema.yaml). See [`docs/usage.md`](usage.md) for field-by-field reference.

Decisions captured under [`docs/adr/`](adr/):
- ADR-001: single conda/mamba env vs per-rule `--use-conda`.
- ADR-002: configurable metadata aggregation strategies across merged ranges.
- ADR-003: retroviral-only pre-filter for LTR_retriever (Coupling A).
- ADR-004: `SPECIES_POST` → Snakemake checkpoint + runtime `species_with_hits(wildcards)` resolver.
- ADR-005: wrap LTR_retriever invocation in a Python runner script (`workflow/scripts/run_ltr_retriever.py`) for testability + log capture + fail-loud behaviour.
- ADR-006: canonicalise genome FASTA filenames to `.fa` via symlink (handles `.fna`/`.fasta`/`.ffn` inputs) — owned by `genome_fasta_normalizer_setup`.

## Design principles

- **Reproducibility**: single pinned env; deterministic rules where possible; stochastic steps (hotspot permutations, Entrez fetches) documented.
- **Resilience**: Snakemake checkpointing; retry logic on Entrez; validator checks tool versions before run.
- **Modularity**: one script per responsibility; strict separation of Python (I/O, orchestration) and R (analysis, plotting).
- **Elegance**: SOLID + DRY; GenomicRanges-first interval ops in R; pathlib-first path handling in Python.
