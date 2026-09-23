# Usage

End-user reference for configuring and running RetroSeek.

## Setup

```bash
make env                    # Create the conda/mamba env
conda activate RetroSeek
```

The env installs BLAST+, GenomeTools, NCBI Datasets CLI, Python 3.11, R 4.3, Bioconductor, and every library the pipeline needs. No system-level tools are required beyond mamba/conda itself.

## CLI

```
./RetroSeek [STAGE FLAGS] [RUN OPTIONS] [SNAKEMAKE OPTIONS]
```

Stage flags choose what to build. All requested stages run as **one** Snakemake
workflow, so their order on the command line does not matter, and a stage whose
inputs are stale pulls its upstream stages in. `./RetroSeek -h` lists them by
phase (the stage table lives in `workflow/scripts/stages.py`, ADR-020).

### Stage flags

**Setup**: fetch and prepare inputs, once per study.

| Flag | Snakemake target | Makes |
|---|---|---|
| `--download-genomes` | `genome_downloader` | `{genome}.fa` in `SPECIES_DB`, from `config.species` accessions |
| `--download-hmm` | `pfam_hmm_downloader` | `Pfam-A.hmm`, the pinned release `input.pfam_release` (ADR-019) |
| `--build-reference` | `taxonomy_reference_trees` | `data/taxonomy_reference/`: reference proteins, NCBI taxonomy, per-gene placement trees (network, built once) |
| `--probe-extractor` | `probe_extractor` | `probe_dict.pkl` (+ CSV / Parquet) from `input.probe_csv` |

**Indexing**: per-genome search indexes.

| Flag | Snakemake target | Makes |
|---|---|---|
| `--blast-dbs` | `blast_db_generator` | BLAST nucleotide database per genome |
| `--suffix-arrays` | `ltr_index_generator` | GenomeTools suffix array per genome |

**Discovery**: the heavy searches (LTRdigest takes about a day per large genome).

| Flag | Snakemake target | Makes |
|---|---|---|
| `--ltr-candidates` | `ltr_harvester` | LTR element candidates (LTRharvest GFF3 + FASTA) |
| `--ltr-domains` | `ltr_digester` | LTR element annotation with LTRdigest: protein domains and polypurine tracts inside each element |
| `--blast` | `blast_pkl2parquet` | tBLASTn of every probe against every genome, `full_genome_blast.parquet` |

**Analysis**: everything that reads the discovery results.

| Flag | Snakemake target(s) | Makes |
|---|---|---|
| `--ranges-analysis` | `ranges_analysis` | Element-hit, orphan and flanking-LTR tracks, overlap matrices, stage tables |
| `--domain-scan` | `domain_scanner` | Curated Pfam domains on element AND orphan loci (`hmmsearch --cut_ga`, ADR-015), plus `Pfam.version` |
| `--classify` | `taxonomy_classify`, `taxonomy_orphans`, `taxonomy_plot_generator`, `loss_analysis` | Per-locus genus calls for both tiers, `catalog.csv`, `tracks/taxonomy/`, `taxonomy.pdf`, `loss.pdf`. Runs the domain scan first. |
| `--segment` | `taxonomy_segments` | The catalog split by taxon at `classification.segment_rank`, one PDF per segment + `overview.pdf` |
| `--solo-ltr-detector` | `solo_ltr_detector` | Solo LTRs: tracks, tables, the evidence tree and one PDF per genome (`docs/solo_ltr.md`) |
| `--hotspot-detection` | `hotspot_detector` | Hotspot CSV + GFF3 + one PDF per genome |
| `--pair-detection` | `pair_detector` | Per-species probe-pair tables |
| `--placement-trees` | `placement_trees` | Heat-trees, EDPL/LWR tables, co-phylogeny against the host tree |

**Figures**

| Flag | Snakemake target(s) | Makes |
|---|---|---|
| `--generate-global-plots` | `plot_generator`, `stage_plot_generator`, `erv_like_plot_generator` | `homology.pdf`, `integration.pdf`, `structure.pdf` |
| `--generate-circle-plots` | `circle_plot_generator` | Per-genome circle plots (currently broken) |

> **Note (ADR-012):** hotspot detection counts **integration events** from the
> per-locus catalog (`hotspot.input: catalog`), not tBLASTn hits, so a stale
> `catalog.csv` pulls `--classify` into the run. `hotspot.source: both`, a larger
> `hotspot.window_size` or `hotspot.input: original` recover density.

**Taxonomic classification** needs the reference: `./RetroSeek --build-reference`
(or `make reference`) fetches it once into `data/taxonomy_reference/`; `--classify`
builds it first if it is missing. `--segment` rolls each call up to
`classification.segment_rank`; a call coarser than that rank is reported as
`unassigned_at_<rank>`. Set `input.species_tree` to a Newick of the host phylogeny
to put species in tree order.

### Presets and run options

| Option | Meaning |
|---|---|
| `--downstream` | Every Analysis and Figures stage except circle plots: the usual run once the discovery searches exist. |
| `-skp`, `--skip-validation` | Skip the slow checks: NCBI lookups of every probe accession and the prompts. The fast checks always run. |
| `--allow-heavy` | Let a heavy rule run although its own stage was not requested (see below). |
| `--stop-on-error` | Stop at the first failed job. By default Snakemake's `--keep-going` lets independent jobs (other genomes) finish. |
| `--config-help [KEY]` | Print the documentation of one config field, or list them all, and exit. |

### What happens before anything runs

1. **Fast checks, always**: the config against `schema.yaml`, the tools the chosen
   stages call, and, for `--domain-scan` and `--classify`, that `Pfam-A.hmm` holds
   every family in the curated class table.
2. **Slow checks, unless `-skp`**: NCBI lookups of the probe accessions and the API
   key prompt (only for stages that talk to NCBI), then a confirmation prompt. The
   prompts fall back to their defaults when no terminal is attached.
3. **The heavy-rule guard**: a dry run of the requested stages. If it would run a
   heavy rule (the downloads, suffix arrays, LTRharvest, LTRdigest, tBLASTn or the
   probe fetch) that no requested stage owns, nothing runs; the message names each
   rule, Snakemake's reason and the usual fix. A newly downloaded `Pfam-A.hmm`, for
   example, makes every LTRdigest output look stale (ADR-019).

The command exits non-zero whenever a check, the guard or a job failed.

### Reading a run

The screen opens with a banner (commit, config, genomes, cores, stages, the run
log's path), then the checks, then the run: one line per finished job (`OK`), every
warning and error, and a progress bar. It closes with a summary:

```
status    failed (exit 1)
time      7h41m
warnings  12   (all of them: logs/runs/2026-09-23_220114.warnings.txt)
  solo_finder: N baits under N bp (9)
  hotspot_detector: no windows passed (3)
failed    1
  taxonomy_classify Myotis_lucifugus: blastx failed with exit code 2. Fix: its error output is in the job log
    log: logs/taxonomy_classify/Myotis_lucifugus.log
run log   logs/runs/2026-09-23_220114.log
Fix the cause, then rerun the same command: finished work is kept.
```

How much the screen shows is `display.verbosity` (`quiet`, `normal`, `verbose`),
or `--verbosity` for one run. Nothing is lost at any setting: every job keeps its
own log, `LOG_DIR/<step>/<genome>.log`, and the run log keeps everything the run
printed. The line format, the levels and the colours are described in
[console_style.md](console_style.md).

### Snakemake options

Any other option goes to Snakemake unchanged, after the launcher's own:

- `-n` / `--dry-run`: show what would run. The guard's findings are printed too.
- `--configfile <path>`: the study config (see below).
- `--forcerun RULE...`: put it last; it takes every name after it.
- `--rerun-triggers mtime`: decide reruns on file dates only.
- `--unlock`, `--cleanup-metadata FILE...`: maintenance; the guard stays out of these.
- `--profile <name>`, `--latency-wait N`: cluster use.

### Inspecting config fields

`./RetroSeek --config-help` prints the field reference in the terminal and exits
without running anything (no validation, no Snakemake, no directories created).
Pass a key for one field, or omit it to list every field:

```bash
./RetroSeek --config-help                 # list every field, grouped by section
./RetroSeek --config-help merge_option    # one field: type, default, meaning
./RetroSeek --config-help classification.placement_genes
```

The text is sourced directly from [`docs/configuration.md`](configuration.md), so
the terminal help and the written reference cannot diverge.

## Configuration

### Local overrides for production paths

`data/config/config.yaml` ships with **repo-relative defaults** (`data/`, `results/`, `logs/`) so a fresh clone runs portably and the toy-genome smoke tests work out of the box. For production runs pointing at external storage:

```bash
cp data/config/config.example.yaml data/config/config.local.yaml
# edit data/config/config.local.yaml - set the four `root` paths and input.probe_csv
./RetroSeek --probe-extractor --configfile data/config/config.local.yaml
```

`config.local.yaml` is in `.gitignore`. It must be a **complete** config, not a list of changes: the pipeline reads the file given to `--configfile` on its own and never falls back to `config.yaml` for missing fields. Start from a copy and edit it. Absolute paths in the local config are honoured as-is; relative paths resolve against the repo root.

### Pipeline config - [`data/config/config.yaml`](../data/config/config.yaml)

`config.yaml` is **values-only** - every field's type, default, and meaning lives in [`docs/configuration.md`](configuration.md) (the canonical reference), also reachable via `./RetroSeek --config-help [KEY]`. Top-level sections at a glance:

- **`blast`** - `e_value`, `optional_parameters`.
- **`genome_tools`** - `suffix_array_parts`, per-subcommand optional parameters.
- **`parameters`** - core thresholds and filters:
  - `identity_threshold`, `bitscore_threshold` - BLAST hit filters.
  - `probe_min_length` - per-probe minimum alignment length.
  - `main_probes` - probes subject to Pfam-domain validation.
  - `merge_option` - how overlapping ranges collapse (`virus` or `label`, strict enum).
  - `aggregation` - per-field strategy (`list` / `concatenate` / `best` / `majority` / `first` / `strict`) applied when merged ranges collapse. See [`docs/configuration.md`](configuration.md#aggregation-strategies) for the vocabulary and [ADR-002](adr/ADR-002-aggregation-strategies.md) for the rationale.
  - Pair settings: `probe_to_pair`, `pair_max_gap`.
  - (The composite ERV assembly is no longer a `parameters.erv_like` tier - it is now the genus-founded loci table from the `classification` stage; the erv-like plot panel reads that table.)
- **`hotspot`** - deterministic NB-GLM hotspot detection (its own top-level config section): `input` (`catalog` | `original`), `group_by`, `source`, `window_size`, `mask_size` / `mask_mismatch`, `pvalue_threshold`, `min_hits`, `merge_gap`, `strata_by_chromosome`, `unplaced_min_factor`. See [`docs/configuration.md`](configuration.md#hotspot).
- **`solo_ltr`** - native solo-LTR detection (ADR-017): the acceptance thresholds (`min_bait_length`, `min_hit_length`, `min_identity`, the coverage window), the monoLTR-at-orphan distance (`orphan_pad`), the blastn settings, and a `tree:` subblock for the evidence phylogeny. Every value the method depends on lives here; the scripts carry no defaults of their own. See [`docs/configuration.md`](configuration.md#solo_ltr) and [`docs/solo_ltr.md`](solo_ltr.md).
- **`placement`** - colour scale for the published heat-trees: `mass_norm` (`absolute` | `relative`).
- **`classification`** - per-locus ERV taxon calls, rank-agnostic since ADR-008 (`reference_taxa` sets the axis, `segment_rank` the roll-up): `enable`, `placement_genes` (default `[POL, GAG, ENV]`), `search` (`blastx`), `evalue`, `top_percent` (weighted-LCA band), `min_orf`, `confidence_min`, `structure_full_min`, `segment_rank`, `reference_taxa`. Reuses `parameters.seed` / `parameters.main_probes` / `execution.entrez_email`. See [`docs/configuration.md`](configuration.md#classification) and [ADR-007](adr/ADR-007-taxonomic-classification.md).
- **`plots`** - segment page selection, axis scales, Sankey and waffle settings, page growth per genome (`per_stratum`).
- **`execution`** - parallelism and API politeness:
  - `num_cores`, `max_threadpool_workers`.
  - `retrieval_time_lag` (Entrez delay), `max_retrieval_attempts` (retries).
  - `entrez_email` - **required** (NCBI ToS).
- **`input`** - `probe_csv` (path to your probe metadata CSV; relative paths resolve against the repo root), `species_tree`, `pfam_domain_classes` and `pfam_release`.
- **`display`** - `verbosity` (`quiet` | `normal` | `verbose`): how much the terminal shows.
- **`root`** - base directories for DB, data, results, logs.
- **`species`** - map of genome ID -> scientific name.

### Validation - [`data/config/schema.yaml`](../data/config/schema.yaml)

`schema.yaml` defines types, ranges, and enum constraints (e.g., `merge_option` must match `^(virus|label)$`). The launcher checks the config against it before every run, `-skp` or not (`validator.py::preflight`); an unknown key or a wrong type stops the run with the offending field named.

### Probe CSV

The path specified by `config.input.probe_csv` points to a CSV describing probes. Expected columns are parsed by `workflow/scripts/probe_extractor.py::table_parser()`. A template lives under `data/tables/_input/` (not tracked - user-provided).

Probe name strings are **uppercased** on load; downstream comparisons (including config matching) are case-sensitive. Use uppercase in `config.parameters.main_probes`, `config.parameters.probe_min_length`, `config.domains`, and `config.parameters.probe_to_pair`.

## Resuming after interruption

Snakemake tracks rule completion. After a crash, re-running the same stage flag skips completed outputs and continues from the last checkpoint. No extra action needed.

The pipeline also uses a true Snakemake **checkpoint** (`blast_pkl2parquet`) to defer downstream DAG evaluation until the BLAST stage has produced its parquet; rules that used to depend on a parse-time `SPECIES_POST` list now resolve their inputs via a runtime `species_with_hits(wildcards)` call against the checkpoint's output. See [ADR-004](adr/ADR-004-species-post-checkpoint.md).

## Observability on long-running rules

`gt suffixerator` and `gt ltrharvest` write their primary outputs to stdout (redirected to files) and emit nothing to stderr, so their Snakemake logs are silent for 30-90 min on mammalian genomes. RetroSeek wraps both rules in a trap-backed background heartbeat that emits `[heartbeat:<rule>:<genome>] still running at Nm elapsed` to stderr every 60 s. Useful for distinguishing a live suffixerator run from a wedged one during multi-hour executions. No configuration needed - the heartbeat is unconditional and zero-cost when rules are fast.

## Troubleshooting

- **Entrez soft-bans** - ensure `execution.retrieval_time_lag >= 0.3` and `execution.entrez_email` is a real address. Optional: register an NCBI API key.
- **OOM in `ltr_index_generator`** - raise `genome_tools.suffix_array_parts`.
- **`merge_option` validation failure** - check for typos; must be exactly `virus` or `label`.
- **Silent missing species** - if `use_species_dict: false`, ensure `{genome}.fa` exists under `SPECIES_DB`. If `true`, ensure genome IDs in `config.species` match the expected file stems.
