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
**Downstream analyses** — `plot_generator`, `circle_plot_generator`, `hotspot_detector`, `pair_detector`.

## Outputs

- `results/tracks/` — GFF3 per stage (`original`, `candidate`, `valid`, `flanking_ltr`, `ltrdigest`, `ltrharvest`, `solo_ltr`), each with a `_reduced` variant (merged overlapping ranges).
- `results/tables/` — overlap matrices, plot dataframes (Parquet), probe-pair tables, segmented species data.
- `results/plots/` — density, raincloud, bar, Sankey, balloon, Circos-style PNG + PDF.
- `results/hotspots/` — CSV + GFF3 + PDFs (histogram / density / heatmap).

## Key configuration

All user-tunable parameters live in [`data/config/config.yaml`](../data/config/config.yaml); validation rules are in [`data/config/schema.yaml`](../data/config/schema.yaml). See [`docs/usage.md`](usage.md) for field-by-field reference.

Decisions captured under [`docs/adr/`](adr/):
- ADR-001: single conda/mamba env vs per-rule `--use-conda`.

## Design principles

- **Reproducibility**: single pinned env; deterministic rules where possible; stochastic steps (hotspot permutations, Entrez fetches) documented.
- **Resilience**: Snakemake checkpointing; retry logic on Entrez; validator checks tool versions before run.
- **Modularity**: one script per responsibility; strict separation of Python (I/O, orchestration) and R (analysis, plotting).
- **Elegance**: SOLID + DRY; GenomicRanges-first interval ops in R; pathlib-first path handling in Python.
