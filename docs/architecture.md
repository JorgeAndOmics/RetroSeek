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
**Integration & segmentation** — `species_segmenter`, `ranges_analysis`.
**Solo-LTR detection** — `ltr_retriever_prefilter`, `ltr_retriever`, `solo_ltr_integrator`, `solo_ltr_detector` (aggregate). Runs LTR_retriever over LTRharvest output pre-filtered by `valid_ranges.gff3`; propagates RetroSeek probe labels onto discovered solo LTRs. See [`docs/solo_ltr.md`](solo_ltr.md) for the full mechanism and [ADR-003](adr/ADR-003-ltr-retriever-pre-filter.md) for the pre-filter rationale.
**Downstream analyses** — `plot_generator`, `circle_plot_generator`, `hotspot_detector`, `pair_detector`.

## Outputs

- `results/tracks/` — GFF3 per stage (`original`, `candidate`, `valid`, `flanking_ltr`, `ltrdigest`, `ltrharvest`, `solo_ltr`, `ltr_retriever/`), each with a `_reduced` variant (merged overlapping ranges) where applicable. `solo_ltr/{genome}.gff3` carries probe_labels propagated from valid ERVs (see `docs/solo_ltr.md`).
- `results/tables/` — overlap matrices, plot dataframes (Parquet), probe-pair tables, segmented species data, solo/intact ratios (`solo_intact_ratio/{genome}.csv` + `all_species.csv`).
- `results/plots/` — density, raincloud, bar, Sankey, balloon, Circos-style PNG + PDF.
- `results/hotspots/` — CSV + GFF3 + PDFs (histogram / density / heatmap).
- `data/ltr_scn/` — LTRharvest screen-format intermediates (`{genome}.scn` + `{genome}_retroviral.scn`); consumed by LTR_retriever.

## Key configuration

All user-tunable parameters live in [`data/config/config.yaml`](../data/config/config.yaml); validation rules are in [`data/config/schema.yaml`](../data/config/schema.yaml). See [`docs/usage.md`](usage.md) for field-by-field reference.

Decisions captured under [`docs/adr/`](adr/):
- ADR-001: single conda/mamba env vs per-rule `--use-conda`.
- ADR-002: configurable metadata aggregation strategies across merged ranges.
- ADR-003: retroviral-only pre-filter for LTR_retriever (Coupling A).
- ADR-004: `SPECIES_POST` → Snakemake checkpoint + runtime `species_with_hits(wildcards)` resolver.

## Design principles

- **Reproducibility**: single pinned env; deterministic rules where possible; stochastic steps (hotspot permutations, Entrez fetches) documented.
- **Resilience**: Snakemake checkpointing; retry logic on Entrez; validator checks tool versions before run.
- **Modularity**: one script per responsibility; strict separation of Python (I/O, orchestration) and R (analysis, plotting).
- **Elegance**: SOLID + DRY; GenomicRanges-first interval ops in R; pathlib-first path handling in Python.
