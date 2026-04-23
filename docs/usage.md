# Usage

End-user reference for configuring and running RetroSeek.

## Setup

```bash
make env                    # Create the conda/mamba env
conda activate retroseek
```

The env installs BLAST+, GenomeTools, NCBI Datasets CLI, Python 3.10, R 4.3, Bioconductor, and every library the pipeline needs. No system-level tools are required beyond mamba/conda itself.

## CLI

```
./RetroSeek [STAGE_FLAG] [SNAKEMAKE_FLAGS]
```

Exactly one **stage flag** selects which pipeline section runs. Any additional arguments are passed through to Snakemake.

### Stage flags

| Flag                        | Snakemake target              | Inputs (short)                          | Outputs (short)                              |
|-----------------------------|-------------------------------|-----------------------------------------|----------------------------------------------|
| `--download-genomes`        | `genome_downloader`           | `config.species` accessions             | `{genome}.fa` in `SPECIES_DB`                |
| `--download-hmm`            | `pfam_hmm_downloader`         | —                                       | `Pfam-A.hmm`                                 |
| `--blast-dbs`               | `blast_db_generator`          | `{genome}.fa`                           | BLAST nucleotide DB files                    |
| `--suffix-arrays`           | `ltr_index_generator`         | `{genome}.fa`                           | Suffix-array index files                     |
| `--ltr-candidates`          | `ltr_harvester`               | Suffix arrays                           | LTR candidate GFF3 + FASTA                   |
| `--ltr-domains`             | `ltr_digester`                | LTR GFF3 + Pfam HMMs + FASTA            | Domain-annotated LTR GFF3                    |
| `--probe-extractor`         | `probe_extractor`             | `config.input.probe_csv`                | `probe_dict.pkl` (+ CSV / Parquet)           |
| `--blast`                   | `full_genome_blaster`         | BLAST DBs + probe dict                  | `{genome}.pkl` (tBLASTn results)             |
| `--ranges-analysis`         | `ranges_analysis`             | BLAST Parquet + LTRdigest GFF3          | GFF3 tracks + overlap matrices + dataframes  |
| `--generate-global-plots`   | `plot_generator`              | Plot dataframes                         | PNG plots (density, raincloud, bar, Sankey)  |
| `--generate-circle-plots`   | `circle_plot_generator`       | Valid GFF3 + LTRdigest GFF3 + FASTA     | Per-genome Circos-style PNG + PDF            |
| `--hotspot-detection`       | `hotspot_detector`            | Original GFF3 tracks + FASTA            | Hotspot CSV + GFF3 + histogram/density PDFs  |
| `--pair-detection`          | `pair_detector`               | Valid ranges GFF3                       | Per-species pair tables (CSV + Parquet)      |
| `-skp`, `--skip-validation` | —                             | —                                       | Bypass pre-run validation (debug only)       |

### Snakemake pass-through

Any argument not recognised as a stage flag is forwarded. Common examples:

- `--cores N` (or `--cores all`) — parallelism.
- `--profile <name>` — HPC/cluster profile.
- `--keep-going` — continue on rule failure.
- `--latency-wait N` — filesystem latency tolerance.
- `--dry-run` / `-n` — DAG-only, no execution.
- `--configfile <path>` — override default config.

## Configuration

### Pipeline config — [`data/config/config.yaml`](../data/config/config.yaml)

Top-level sections (fields inside each section — see file for full list):

- **`blast`** — `e_value`, `optional_parameters`.
- **`genome_tools`** — `suffix_array_parts`, per-subcommand optional parameters.
- **`parameters`** — core thresholds and filters:
  - `identity_threshold`, `bitscore_threshold` — BLAST hit filters.
  - `probe_min_length` — per-probe minimum alignment length.
  - `main_probes` — probes subject to Pfam-domain validation.
  - `merge_option` — how overlapping ranges collapse (`virus` or `label`, strict enum).
  - Hotspot settings: `hotspot_window_size`, `hotspot_permutations`, `hotspot_pvalue_threshold`, etc.
  - Pair settings: `probe_to_pair`, `pair_max_gap`.
- **`logging`** — colour styles for console logging.
- **`plots`** — DPI, dimensions, Sankey omission threshold, circle-plot bitscore cutoff.
- **`execution`** — parallelism and API politeness:
  - `num_cores`, `max_threadpool_workers`.
  - `retrieval_time_lag` (Entrez delay), `max_retrieval_attempts` (retries).
  - `entrez_email` — **required** (NCBI ToS).
- **`input`** — `probe_csv` (absolute path to your probe metadata CSV).
- **`display`** — verbosity toggles.
- **`root`** — base directories for DB, data, results, logs.
- **`domains`** — per-probe Pfam domain/regex lists used in validation.
- **`species`** — map of genome ID → scientific name.

### Validation — [`data/config/schema.yaml`](../data/config/schema.yaml)

`schema.yaml` defines types, ranges, and enum constraints (e.g., `merge_option` must match `^(virus|label)$`). `validator.py::validation_run()` checks the config against this schema before any stage runs (unless `--skip-validation` is passed).

### Probe CSV

The path specified by `config.input.probe_csv` points to a CSV describing probes. Expected columns are parsed by `workflow/scripts/probe_extractor.py::table_parser()`. A template lives under `data/tables/` (not tracked — user-provided).

Probe name strings are **uppercased** on load; downstream comparisons (including config matching) are case-sensitive. Use uppercase in `config.parameters.main_probes`, `config.parameters.probe_min_length`, `config.domains`, and `config.parameters.probe_to_pair`.

## Resuming after interruption

Snakemake tracks rule completion. After a crash, re-running the same stage flag skips completed outputs and continues from the last checkpoint. No extra action needed.

## Troubleshooting

- **Entrez soft-bans** — ensure `execution.retrieval_time_lag ≥ 0.3` and `execution.entrez_email` is a real address. Optional: register an NCBI API key.
- **OOM in `ltr_index_generator`** — raise `genome_tools.suffix_array_parts`.
- **`merge_option` validation failure** — check for typos; must be exactly `virus` or `label`.
- **Silent missing species** — if `use_species_dict: false`, ensure `{genome}.fa` exists under `SPECIES_DB`. If `true`, ensure genome IDs in `config.species` match the expected file stems.
- **`git-lfs` not installed** — repo uses LFS for `*.pkl`. Operations on code files still work; for LFS-backed files install `git-lfs` or bypass hooks (`git -c core.hooksPath=/dev/null ...`).
