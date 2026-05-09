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
                                             ├─► GFF3 tracks (original / candidate / valid + _reduced)
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
**Integration & segmentation** — `species_segmenter`, `ranges_analysis` (phase-modular: `workflow/scripts/ranges_analysis.R` is a thin orchestrator over sibling modules in `workflow/scripts/range_analysis/{io,granges_build,filtering,reductions,validation,plot_dataframe,exporters}.R`).
**Solo-LTR detection** — `ltr_retriever_prefilter`, `ltr_retriever`, `solo_ltr_integrator`, `solo_ltr_detector` (aggregate). Runs LTR_retriever over LTRharvest output pre-filtered by `valid_ranges.gff3`; propagates RetroSeek probe labels onto discovered solo LTRs. See [`docs/solo_ltr.md`](solo_ltr.md) for the full mechanism and [ADR-003](adr/ADR-003-ltr-retriever-pre-filter.md) for the pre-filter rationale.
**Downstream analyses** — `plot_generator` (phase-modular: `workflow/scripts/plot2sort.R` is a slim orchestrator over sibling modules in `workflow/scripts/plot2sort/{helpers,io,plots_distribution,plots_categorical,plots_sankey}.R`; emits 20 PNGs with auto-scaled canvases driven by `auto_dims()`, opt-in long-tail collapse via `sankey_top_n`, and waffle auto-scaling), `circle_plot_generator`, `hotspot_detector` (Negative-Binomial GLM on per-window counts, mask-aware via `log(effective_bp)` offset and chromosome-stratified; phase-modular under `workflow/scripts/hotspot_analysis/{io,masking,windowing,models,postprocess,plots}.R`; `regioneR` permutation pass available as opt-in validation via `parameters.hotspot_validate_permutation`), `pair_detector`.

## Outputs

- `results/tracks/` — GFF3 per stage (`original`, `candidate`, `valid`, `flanking_ltr`, `ltrdigest`, `ltrharvest`, `solo_ltr`, `ltr_retriever/`, `hotspots/`), each with a `_reduced` variant (merged overlapping ranges) where applicable. `solo_ltr/{genome}.gff3` carries probe_labels propagated from valid ERVs (see `docs/solo_ltr.md`). `hotspots/{genome}.gff3` + `hotspots/{genome}.bed` carry merged hotspot regions with NB GLM q-values.
- `results/tables/` — overlap matrices, plot dataframes (Parquet), probe-pair tables, segmented species data, solo/intact ratios (`solo_intact_ratio/{genome}.csv` + `all_species.csv`), hotspot per-window tables (`hotspots/{genome}.csv` + `.parquet` + `.manifest.yaml`).
- `results/plots/` — 20 global PNGs from `plot2sort` (density, raincloud, query-coverage, bar, balloon, heatmap-probe×species, virus waffle, three Sankey variants, each split into `main` / `accessory` / `full` where applicable), plus per-genome Circos-style PNG + PDF, and hotspot Manhattan / karyotype / Q-Q / summary PDFs (+ histogram / density when permutation validation is on).
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

- **Reproducibility**: single pinned env; deterministic rules where possible; stochastic steps (the opt-in hotspot permutation pass and Entrez fetches) seeded via the global `parameters.seed` (default `67`) where seedable.
- **Resilience**: Snakemake checkpointing; retry logic on Entrez; validator checks tool versions before run.
- **Modularity**: one script per responsibility; strict separation of Python (I/O, orchestration) and R (analysis, plotting).
- **Elegance**: SOLID + DRY; GenomicRanges-first interval ops in R; pathlib-first path handling in Python.
