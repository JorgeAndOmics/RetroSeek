# Configuration reference

Every field in [`data/config/config.yaml`](../data/config/config.yaml) documented. Validation rules live in [`data/config/schema.yaml`](../data/config/schema.yaml) and are enforced by `validator.py::validation_run()` at pipeline start.

This page is the **canonical reference**: `config.yaml` itself is values-only, and `RetroSeek --config-help [KEY]` prints these field descriptions in the terminal by reading the tables below. Keep this file current and both the config and the CLI follow.

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
| `seed` | int ≥ 0 | `67` | Master RNG seed for reproducibility. Recorded in the run manifest so a run's stochastic steps (e.g. permutation tests) can be reproduced. |
| `identity_threshold` | int ≥ 0 | `0` | Minimum % identity for BLAST hits. `0` disables the filter. |
| `bitscore_threshold` | number ≥ 0 | `0` | Minimum bit score. `0` disables the filter. |
| `ltr_resize` | int ≥ 0 | `0` | Padding (bp) added to each LTR retrotransposon on both sides before overlap detection. |
| `ltr_flank_margin` | int ≥ 0 | `0` | Tolerance (bp) used when classifying flanking LTRs as left vs right. |
| `merge_option` | `virus` \| `label` | `virus` | How overlapping ranges group before `plyranges::reduce_ranges_directed`. **Strict enum** — typos fail validation. |
| `main_probes` | list of strings | `[POL, GAG, ENV, PRO]` | Probe names treated as *main* (as opposed to *accessory*). Semantically a set — duplicates ignored. Drives the `probe_type` column on plot dataframes and the `probe_category` attribute on GFF3 tracks. |
| `probe_min_length` | map (string → int) | `{ GAG: 200, POL: 400, ... }` | Per-probe minimum alignment length in residues. Ranges shorter than the probe-specific threshold are filtered out. |

### ERV-like assembly

Chains ≥2 *distinct* **main** probe loci from the **unreduced** `valid` tier into composite ERV-like candidates — candidate conserved full ERVs, or the longest recoverable fragment (e.g. GAG+POL when ENV is absent). Purely additive: the `valid` output is unchanged and isolated single-gene loci stay there. Output: `results/tracks/erv_like/{genome}.gff3` (a parent `erv_like` feature per candidate + its child `erv_like_member` loci carrying `Parent=`), a child-locus `.bed`, and the `{genome}.erv_like_loci` table.

| Key | Type | Default | Meaning |
|---|---|---|---|
| `erv_like.group_by` | `virus` \| `label` \| `none` | `virus` | Which loci may chain together. `virus`/`label` only chain probes sharing that attribute, so co-located different-group probes yield separate (possibly overlapping) candidates — **both retained**. `none` chains across all main loci in the window. **Strict enum** — typos fail validation. |
| `erv_like.max_join_distance` | int ≥ 0 | `1500` | Maximum gap (bp) between adjacent main-probe loci to still chain them. Inclusive (a gap exactly equal to this still joins). |
| `erv_like.require_canonical_order` | bool | `false` | When `true`, keep only candidates whose main probes occur in the order given by `main_probes` (forward on `+`, reversed on `-`, either on `*`). The `main_probes` list order *is* the canonical gene order. Dropped non-canonical candidates are still tallied (`erv_like_dropped_noncanonical` in `counts`) for the canonical-vs-rearranged plot. |
| `erv_like.completeness_threshold` | number 0–1 | `1.0` | `is_full` = (`n_main_present` / number of `main_probes`) ≥ this. `1.0` requires every main probe. |

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
| `list` | native multi-value (`CharacterList` in R; `list<string>` in parquet; comma-separated in GFF3 attributes) | Every unique contributor preserved — no information lost. Inflates plot row counts (entry explosion); the plot scripts emit a warning when `virus`/`label` use `list` or `concatenate`. |
| `concatenate` | single string joined by `concat_separator` | Pre-refactor behaviour. Human-readable in a text editor. Same plot-inflation caveat as `list`. |
| `best` | single value | Row with the highest `best_tiebreaker` wins. **Default** for `virus`/`label`. Tie-break is deterministic — see note below. |
| `majority` | single value | Mode (most frequent value). Ties broken by first alphabetical. |
| `first` | single value | Alphabetical first unique value. Deterministic but arbitrary. |
| `strict` | single value or `strict_marker` | Returns the value if all contributors agree, else the configured marker. |

#### `parameters.aggregation`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `aggregation.virus` | strategy | `best` | Strategy for `virus` column. |
| `aggregation.label` | strategy | `best` | Strategy for `label` column. |
| `aggregation.probe` | strategy | `list` | Strategy for `probe` column. Always a reduction grouping key, so single-valued in practice — strategy is effectively a no-op. |
| `aggregation.species` | strategy | `first` | Per-genome pipeline — all rows share one species. |
| `aggregation.best_tiebreaker` | `bitscore` \| `identity` \| `align_length` | `bitscore` | Column used to rank contributors when strategy is `best`. |
| `aggregation.concat_separator` | string | `"; "` | Separator used by `concatenate`. |
| `aggregation.strict_marker` | string | `"ambiguous"` | Value emitted by `strict` when contributors disagree. |

**Deterministic `best` tie-break.** When several contributors to a merged
range tie on `best_tiebreaker`, the winner is resolved by a fixed key chain
applied as a sort of the reduction input: `bitscore` desc → `query_coverage`
desc → `identity` desc → `evalue` asc → genomic position → `label` name. This
makes `best` reproducible across R / plyranges versions.

#### `parameters.solo_ltr_aggregation`

Separate block for solo-LTR probe-label propagation (produced by the LTR_retriever workstream).

| Key | Type | Default | Meaning |
|---|---|---|---|
| `solo_ltr_aggregation.probe` | strategy | `list` | Strategy for propagating probe labels from contributing ERVs to solo LTRs. |
| `solo_ltr_aggregation.best_tiebreaker` | `consensus_members` \| `bitscore` \| `identity` | `consensus_members` | Column used when strategy is `best`. `consensus_members` = number of ERVs that seeded the consensus family. |

#### Choosing a strategy

- **Default (`best`)** — recommended for `virus`/`label`. One deterministic "dominant" value per merged range based on alignment strength; clean single-value GFF3/parquet output and statistically meaningful plots.
- **`list`** — pick this if downstream code needs the full contributor set. Note: inflates plot row counts (entry explosion) — the plot scripts warn when `virus`/`label` use `list`.
- **`concatenate`** — pick this if your downstream tooling expects semicolon-delimited strings (e.g., grep-based inspection, legacy scripts). Same plot caveat as `list`.
- **`majority`** — pick this if you trust hit counts more than hit strength.
- **`strict`** — pick this when you want to flag ambiguity explicitly; useful for high-confidence result tables.
- **`first`** — mostly for testing/comparison; rarely the right production choice.

## `ltr_retriever`

Solo-LTR post-processing of LTRharvest output. See [`docs/solo_ltr.md`](solo_ltr.md) for the full mechanism (how LTR_retriever works end-to-end, what "family" means, how the pre-filter + label-propagation couplings with RetroSeek work, and the biology of solo LTRs) and [ADR-003](adr/ADR-003-ltr-retriever-pre-filter.md) for the pre-filter decision rationale.

| Key | Type | Default | Meaning |
|---|---|---|---|
| `substitution_rate` | number ≥ 0 | `1.3e-8` | bp substitutions per site per year used by LTR_retriever for age estimation. Mammals: `1.3e-8`; plants: `7e-9`. Does not affect solo-LTR detection sensitivity — only age annotations. |
| `min_ltr_similarity` | number 0–100 | `91` | LTR pair similarity floor (percent) for LTR_retriever's intact-ERV filter (`-miniden` flag). |
| `threads_per_genome` | int ≥ 1 | `4` | CPU threads per-genome LTR_retriever invocation. |
| `noanno` | bool | `true` | Skip LTR_retriever's internal TE-library annotation (`-noanno` flag). RetroSeek has its own probe-based classification. |
| `source_scn` | str (`retroviral` \| `full`) | `retroviral` | **Coupling A toggle.** Picks which SCN feeds LTR_retriever. `retroviral` (default) uses the prefilter-restricted SCN — rows overlapping `valid_ranges.gff3` — guaranteeing retroviral-only consensus families. `full` uses the unfiltered LTRharvest passthrough, useful for debugging or non-retroviral exploration. The prefilter rule always materialises both SCN files in `data/ltr_scn/` regardless of this setting. |
| `nearest_erv_max_distance` | int ≥ 0 | `10000` | Bp window for the solo-LTR → valid-ERV **nearest-ERV fallback** in Coupling B's label-propagation. Only used when the primary consensus-family path yields no labels for a given solo LTR. |

Related: `parameters.solo_ltr_aggregation` (already documented above under the `parameters` section) controls the strategy for summarising probe labels inherited from multiple contributing ERVs.

## `logging`

`level_styles` and `field_styles` are passed through to `coloredlogs`. See `coloredlogs.install()` documentation for accepted style dicts. Keys: `color`, `bold`, `background`.

## `plots`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `dpi` | int ≥ 1 | `300` | Output resolution. |
| `width` | int ≥ 1 | `15` | Plot width in inches. |
| `height` | int ≥ 1 | `12` | Plot height in inches. |
| `bitscore_x_scale` | str (`linear` \| `log10`) | `linear` | X-axis scale on density / raincloud bitscore plots. Use `log10` when the long-tail of low-bitscore hits crushes the lower modes. |
| `sankey_top_n` | int ≥ 1 or `null` | `null` | Long-tail handling for Sankey plots. `null` (default) shows every stratum. Set to a positive integer N to keep the top N strata per axis and fold the rest into a single labelled `Other (k)` stratum recording how many strata were collapsed. |
| `sankey_other_label` | str | `Other` | Label prefix for the bundled-tail stratum. The actual rendered label is `<prefix> (k)` where `k` is the number of folded strata. |
| `waffle_unit_hits` | int ≥ 1 | `1` | Number of ranges represented by one waffle square. Bump on huge inputs (e.g. `10` → "1 square = 10 ranges"). |
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
  example_genome_1: "Example species one"
```
