# Configuration reference

Every field in [`data/config/config.yaml`](../data/config/config.yaml) documented. Validation rules live in [`data/config/schema.yaml`](../data/config/schema.yaml) and are enforced by `validator.py::validation_run()` at pipeline start.

## `blast`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `e_value` | number ≥ 0 | `0.01` | E-value threshold for `tblastn`. Lower = stricter. |
| `optional_parameters` | string | `""` | Extra flags appended to every `tblastn` invocation (e.g., `-word_size 3`). |

## `genome_tools`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `suffix_array_parts` | int | `30` | `gt suffixerator -parts` value. Higher = lower memory, more temp files. Tune to 20–40 for mammalian genomes. |
| `suffixerator_optional_parameters` | string | `""` | Extra flags for `gt suffixerator`. |
| `ltrharvest_optional_parameters` | string | `""` | Extra flags for `gt ltrharvest` (e.g., `-minlenltr 100 -maxlenltr 1000`). |
| `ltrdigest_optional_parameters` | string | `""` | Extra flags for `gt ltrdigest`. |

## `parameters`

### Thresholds

| Key | Type | Default | Meaning |
|---|---|---|---|
| `identity_threshold` | int ≥ 0 | `0` | Minimum % identity for BLAST hits. `0` disables the filter. |
| `bitscore_threshold` | number ≥ 0 | `0` | Minimum bit score. `0` disables the filter. |
| `ltr_resize` | int ≥ 0 | `0` | Padding (bp) added to each LTR retrotransposon on both sides before overlap detection. |
| `ltr_flank_margin` | int ≥ 0 | `0` | Tolerance (bp) used when classifying flanking LTRs as left vs right. |
| `merge_option` | `virus` \| `label` | `virus` | How overlapping ranges group before `plyranges::reduce_ranges_directed`. **Strict enum** — typos fail validation. |
| `main_probes` | list of strings | `[POL, GAG, ENV, PRO]` | Probe names treated as *main* (as opposed to *accessory*). Semantically a set — duplicates ignored. Drives the `probe_type` column on plot dataframes and the `probe_category` attribute on GFF3 tracks. |
| `probe_min_length` | map (string → int) | `{ GAG: 200, POL: 400, ... }` | Per-probe minimum alignment length in residues. Ranges shorter than the probe-specific threshold are filtered out. |

### Hotspot detection

| Key | Type | Default | Meaning |
|---|---|---|---|
| `hotspot_window_size` | int ≥ 1 | `10000` | Sliding window width (bp). |
| `hotspot_mask_size` | int ≥ 0 | `20` | Length of the all-N mask used to exclude unsequenced regions. |
| `hotspot_mask_mismatch` | int ≥ 0 | `3` | Mismatches allowed within the mask pattern. |
| `hotspot_permutations` | int ≥ 1 | `1` | Number of permutations for `regioneR::permTest`. `1` is testing-only; use `1000+` for real results. |
| `hotspot_pvalue_threshold` | number | `0.05` | p-value cutoff for calling a window a hotspot. |
| `hotspot_group_split` | bool | `false` | Whether to split hotspot groups during permutation testing. |

### Pair detection

| Key | Type | Default | Meaning |
|---|---|---|---|
| `probe_to_pair` | string | `"ENV"` | Probe name used as the anchor for `pair_detector.R`. |
| `pair_max_gap` | number ≥ 0 | `300000` | Maximum bp distance between paired probes. |

### Aggregation strategies

When overlapping ranges are collapsed via `plyranges::reduce_ranges_directed`, their metadata (virus, label, probe, species) must be reduced to a single value *per merged range*. The strategy is configurable per field.

#### Vocabulary

| Strategy | Output shape | Semantics |
|---|---|---|
| `list` | native multi-value (`CharacterList` in R; `list<string>` in parquet; comma-separated in GFF3 attributes) | Every unique contributor preserved. **Default** — no information lost. |
| `concatenate` | single string joined by `concat_separator` | Current pre-refactor behaviour. Human-readable in a text editor. |
| `best` | single value | Row with the highest `best_tiebreaker` wins. |
| `majority` | single value | Mode (most frequent value). Ties broken by first alphabetical. |
| `first` | single value | Alphabetical first unique value. Deterministic but arbitrary. |
| `strict` | single value or `strict_marker` | Returns the value if all contributors agree, else the configured marker. |

#### `parameters.aggregation`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `virus` | strategy | `list` | Strategy for `virus` column. |
| `label` | strategy | `list` | Strategy for `label` column. |
| `probe` | strategy | `list` | Strategy for `probe` column. Usually single-valued in practice; `list` is safe. |
| `species` | strategy | `first` | Per-genome pipeline — all rows share one species. |
| `best_tiebreaker` | `bitscore` \| `identity` \| `align_length` | `bitscore` | Column used to rank contributors when strategy is `best`. |
| `concat_separator` | string | `"; "` | Separator used by `concatenate`. |
| `strict_marker` | string | `"ambiguous"` | Value emitted by `strict` when contributors disagree. |

#### `parameters.solo_ltr_aggregation`

Separate block for solo-LTR probe-label propagation (produced by the LTR_retriever workstream).

| Key | Type | Default | Meaning |
|---|---|---|---|
| `probe` | strategy | `list` | Strategy for propagating probe labels from contributing ERVs to solo LTRs. |
| `best_tiebreaker` | `consensus_members` \| `bitscore` \| `identity` | `consensus_members` | Column used when strategy is `best`. `consensus_members` = number of ERVs that seeded the consensus family. |

#### Choosing a strategy

- **Default (`list`)** — recommended. Downstream code always has the full contributor set and can compute any summary it wants (mode, first, best, etc.) without needing to parse strings.
- **`concatenate`** — pick this if your downstream tooling expects semicolon-delimited strings (e.g., grep-based inspection, legacy scripts).
- **`best`** — pick this if you want a single deterministic "dominant" label per merged range based on alignment strength.
- **`majority`** — pick this if you trust hit counts more than hit strength.
- **`strict`** — pick this when you want to flag ambiguity explicitly; useful for high-confidence result tables.
- **`first`** — mostly for testing/comparison; rarely the right production choice.

## `logging`

`level_styles` and `field_styles` are passed through to `coloredlogs`. See `coloredlogs.install()` documentation for accepted style dicts. Keys: `color`, `bold`, `background`.

## `plots`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `dpi` | int ≥ 1 | `300` | Output resolution. |
| `width` | int ≥ 1 | `15` | Plot width in inches. |
| `height` | int ≥ 1 | `12` | Plot height in inches. |
| `omit_lower_percent` | number ≥ 0 | `0.05` | Sankey paths representing less than this fraction are dropped. |
| `circle_plot_bitscore_threshold` | number ≥ 0 | `0` | Bit-score cutoff for circle-plot display. |

## `execution`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `num_cores` | int ≥ 1 | `64` | Snakemake parallelism cap. |
| `use_species_dict` | bool | `true` | If `true`, `SPECIES` comes from the `species:` map below. If `false`, derive by scanning `SPECIES_DB` for `*.fa` files. |
| `retrieval_time_lag` | number ≥ 0 | `0.3` | Seconds between Entrez API calls. Required > 0 for NCBI politeness. |
| `max_retrieval_attempts` | int ≥ 1 | `9` | Entrez retry count. |
| `max_threadpool_workers` | int ≥ 1 or `null` | `null` | Python `ThreadPoolExecutor` size. `null` = CPU default. |
| `entrez_email` | string | — | **Required.** Contact email for NCBI Entrez API. NCBI ToS requirement. |

## `input`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `probe_csv` | string | — | **Required.** Absolute path to the probe metadata CSV. Columns expected: `Label, Name, Abbreviation, Probe, Accession`. |

## `display`

Toggle verbosity flags used by logging and UI code.

| Key | Type | Default |
|---|---|---|
| `display_snakemake_info` | bool | `false` |
| `display_requests_warning` | bool | `false` |
| `display_operation_info` | bool | `true` |

## `root`

Base directories. Must be absolute paths. Everything under `data/`, `results/`, `logs/` is resolved relative to these.

| Key | Meaning |
|---|---|
| `db_root_folder` | Where genome FASTAs and BLAST / suffix-array indices live. |
| `data_root_folder` | Input data, pickles, tmp, intermediate tables. |
| `results_root_folder` | Output tracks, tables, plots. |
| `logs_root_folder` | Per-rule log files. |

## `domains`

Map of probe name → list of domain patterns. Patterns are treated as case-insensitive regexes by `ranges_analysis.R` when assigning probes to LTRdigest Pfam hits. Example:

```yaml
domains:
  POL:
    - "ase"
    - "RVT_1"
    - "RVT_2"
  GAG:
    - "Gag"
    - "zf"
    - "PTAP"
```

## `species`

Map of genome ID → scientific name. Keys drive `expand(..., genome=SPECIES)` targets in the Snakefile when `execution.use_species_dict: true`.

```yaml
species:
  HLartInt1A: "Artibeus intermedius"
```
