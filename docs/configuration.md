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

## `hotspot`

Per-genome detection of windows enriched for ERV integrations beyond chance. A single **Negative-Binomial GLM** models per-window hit counts with a mask-aware offset (`log(effective_bp)`) and an optional per-chromosome baseline; per-window upper-tail p-values are BH-adjusted, thresholded, and merged into hotspot regions. The model is **deterministic** — identical inputs give identical results — so the run is reproducible by construction; the global `parameters.seed` is set and recorded in the manifest for provenance only (the core result does not draw on it). Outputs: per-window `{genome}.csv`/`.parquet`, a `{genome}.manifest.yaml`, merged-region `{genome}.gff3`/`.bed` tracks, and Manhattan/karyotype/Q-Q/summary PDFs.

| Key | Type | Default | Meaning |
|---|---|---|---|
| `hotspot.input` | `valid` \| `original` | `valid` | Which upstream track tier feeds detection. `valid` = domain-validated reduced loci (the `{genome}_reduced.gff3` track); `original` = raw unvalidated `tblastn` hits. **Strict enum** — typos fail validation. **Power note**: window enrichment needs density. The curated `valid` (default) tier is sparse — expect few/no hotspots and occasional NB convergence failures. Switch to `original` (with `strata_by_chromosome: false`) for genome-wide hotspot calling. |
| `hotspot.group_split` | bool | `false` | When `true`, fit one NB model per retrovirus genus (`mcols$label`); when `false`, pool all hits into one `Ungrouped` model. |
| `hotspot.window_size` | int ≥ 1 | `500000` | Tile width (bp). Sets resolution and the number of windows tested. Defaults to 500 kb because finer windows over a multi-Gb genome are ~99% empty, which collapses the NB dispersion (`theta`) and destroys power; 10 kb yields zero calls on real assemblies. |
| `hotspot.mask_size` | int ≥ 0 | `20` | Length of the all-N run treated as an unsequenceable gap and subtracted from a window's callable `effective_bp`. `0` disables masking. |
| `hotspot.mask_mismatch` | int ≥ 0 | `3` | Non-N bases tolerated inside the N-mask motif. |
| `hotspot.pvalue_threshold` | number 0–1 | `0.05` | BH-adjusted q-value cutoff for calling a window significant. |
| `hotspot.min_hits` | int ≥ 0 | `2` | Minimum total hits in a merged hotspot region (applied after merging). |
| `hotspot.merge_gap` | int ≥ -1 | `0` | bp gap allowed when merging adjacent significant windows. `0` = merge strictly adjacent windows; `>0` = bridge gaps up to that size; `-1` = no merging (each significant window stays its own region). |
| `hotspot.strata_by_chromosome` | bool | `true` | Include chromosome as an NB covariate so each chromosome gets its own baseline rate. The right value is **tier-dependent**: defaults `true` to match the sparse default `valid` tier, where the covariate is needed for the NB to converge (without it, fits fail on low-count genomes). **Set `false` when switching to a dense tier like `original`** on fragmented scaffold-level assemblies, where per-contig baselines instead absorb local enrichment and drain power. |
| `hotspot.unplaced_min_factor` | int ≥ 1 | `10` | Scaffolds shorter than `this × window_size` are pooled into a single `Unplaced` stratum (avoids unstable per-scaffold coefficients). |

## `parameters`

### Thresholds

| Key | Type | Default | Meaning |
|---|---|---|---|
| `seed` | int ≥ 0 | `67` | Master RNG seed for reproducibility. Recorded in the run manifest for provenance and consumed by any stochastic step. (Hotspot detection is deterministic and does not draw on it.) |
| `identity_threshold` | int ≥ 0 | `0` | Minimum % identity for BLAST hits. `0` disables the filter. |
| `bitscore_threshold` | number ≥ 0 | `0` | Minimum bit score. `0` disables the filter. |
| `ltr_resize` | int ≥ 0 | `0` | Padding (bp) added to each LTR retrotransposon on both sides before overlap detection. |
| `ltr_flank_margin` | int ≥ 0 | `0` | Tolerance (bp) used when classifying flanking LTRs as left vs right. |
| `merge_option` | `virus` \| `label` | `virus` | How overlapping ranges group before `plyranges::reduce_ranges_directed`. **Strict enum** — typos fail validation. |
| `hit_domain_mode` | `membership` \| `positional` | `membership` | How the per-hit `domain_hit_class` on LTR-flanked hits is decided. `membership`: a hit is `substring_match` when its own gene has a config-matched (`domains`) Pfam domain **anywhere in its enclosing LTR element** (co-occurrence; cheap) — values `substring_match` / `no_substring_match`. `positional`: `substring_match` only when the hit **physically overlaps** a config-matched domain of its gene (co-localization; stronger), adding a `non_domain` level for hits overlapping no domain. Orthogonal to the per-provirus `domain_tier`, which is always element-wise. **Strict enum.** |
| `main_probes` | list of strings | `[POL, GAG, ENV, PRO]` | Probe names treated as *main* (as opposed to *accessory*). Semantically a set — duplicates ignored. Drives the `probe_type` column on plot dataframes and the `probe_category` attribute on GFF3 tracks. |
| `probe_min_length` | map (string → int) | `{ GAG: 200, POL: 400, ... }` | Per-probe minimum alignment length in residues. Ranges shorter than the probe-specific threshold are filtered out. |

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

## `classification`

Per-locus ERV taxonomic classification — turns each valid LTR-element locus into a calibrated **taxon call** (`taxon_call` + `rank` + confidence + mosaic flag + ERV class) from the locus's own marker sequence, instead of transferring the best-bitscore probe label. The classification is **rank-agnostic** (ADR-008): the *axis* — the taxa a locus can resolve to — is declared (see `reference_taxa`), at whatever rank, so a locus resolves to that rank when its evidence lands on an axis taxon, or backs off to an honest higher rank (`rank`) otherwise. Each gene is classified independently against a pinned, taxon-comprehensive reference: POL/GAG by phylogenetic placement (MAFFT → EPA-ng → gappa) when a tree resolves an axis taxon, weighted-LCA otherwise, presence-diagnostic genes (e.g. REX/TAX) by presence. The per-gene calls are then combined into a locus call and a mosaic composition. The reference is built once by the `taxonomy_reference*` rules (`make reference`); see [`docs/taxonomy_classification/`](taxonomy_classification/) and the ADRs for the design. Reuses `parameters.seed` (placement/tree determinism), `parameters.main_probes` (gene reliability order + mosaic gene set), and `execution.entrez_email` (reference build).

| Key | Type | Default | Meaning |
|---|---|---|---|
| `enable` | bool | `true` | Master switch for the classification stage. When `false`, the `taxonomy_classify` target produces nothing and the pipeline keeps the legacy probe-label provenance only. |
| `placement_genes` | list of str | `[POL]` | Genes classified by phylogenetic placement onto a per-gene reference tree; every other gene uses weighted-LCA. `POL` is the reliable tree; `GAG` is shipped but opt-in (its reference alignment is low-identity ~19.7%, so its placements are low-confidence). Genes here must have a built tree package under `data/taxonomy_reference/trees/`. |
| `search` | str (`blastx`) | `blastx` | Translated-search engine mapping each locus marker region to reference proteins. `blastx` reuses the BLAST+ already in the env (no extra dependency). |
| `evalue` | number ≥ 0 | `0.001` | blastx e-value cutoff for marker → reference hits. |
| `top_percent` | number 0–1 | `0.1` | Weighted-LCA bitscore band: hits within this fraction of the best bitscore per marker vote on the lowest-common-ancestor call. Smaller = stricter (fewer, higher-confidence ancestors). |
| `min_orf` | int ≥ 0 | `30` | Minimum translated marker length (amino acids) for a region to be eligible for phylogenetic placement; shorter markers fall back to weighted-LCA. |
| `confidence_min` | number 0–1 | `0.5` | Confidence floor for the high/low confidence tag. A locus whose call confidence is **below** this value is tagged `LC` (low confidence) in the `confidence_tag` column of the loci/fragments tables; at or above it is `HC`. The threshold is inclusive (`conf == confidence_min` ⇒ `HC`) and applies to every method (placement, weighted-LCA, presence). Raise it to flag more marginal calls. |
| `structure_full_min` | number 0–1 | `1.0` | Minimum gene completeness (fraction of `main_probes` present) for a locus to be catalogued as a **`full`** ERV in the `structure_class` column. `1.0` requires every main gene. A single-main-gene locus is always `gene`; a multi-gene locus present but below this floor is `partial`. Deliberately gene-content only — flanking-LTR evidence stays in the anchoring axis (`source`) and the solo-LTR module, not here. |
| `reference_taxa` | list of str | `[]` | The **classification axis** (ADR-008): the taxa — at **any** rank (genus `Lentivirus`, family `Bornaviridae`, …) — the reference is built at and that a locus can resolve to as a first-class `taxon_call`. Empty or absent derives the axis from the distinct probeset `Label` values, so `Label` seeds the classifier; setting an explicit list decouples the classifier from the probeset. Changing it requires rebuilding the reference (`make reference`). |

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
| `per_stratum` | number ≥ 0 | `0.18` | Inches of canvas added per category past the base canvas on per-species panels. `width`/`height` size a small study; beyond that the canvas grows by this much per extra genome so labels keep their room. Raise it if ticks still crowd at your genome count; `0` disables growth (fixed canvas). |
| `max_dim` | number ≥ 1 | `60` | Hard ceiling in inches for a grown canvas. At 300 dpi, 60 in ≈ 18,000 px — the practical PNG limit. Panels that would exceed it are clamped rather than failing to render. |

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
