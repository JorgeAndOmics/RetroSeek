# Configuration reference

Every field in [`data/config/config.yaml`](../data/config/config.yaml) documented. Validation rules live in [`data/config/schema.yaml`](../data/config/schema.yaml) and are enforced by `validator.py::validation_run()` at pipeline start.

This page is the **canonical reference**: `config.yaml` itself is values-only, and `RetroSeek --config-help [KEY]` prints these field descriptions in the terminal by reading the tables below. Keep this file current and both the config and the CLI follow.

## `blast`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `e_value` | number >= 0 | `0.01` | E-value threshold for `tblastn`. Lower = stricter. |
| `optional_parameters` | string | `""` | Extra flags appended to every `tblastn` invocation (e.g., `-word_size 3`). |

## `genome_tools`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `suffix_array_parts` | int | `30` | `gt suffixerator -parts` value. Higher = lower memory, more temp files. Tune to 20-40 for mammalian genomes. |
| `suffixerator_optional_parameters` | string | `""` | Extra flags for `gt suffixerator`. |
| `ltrharvest_optional_parameters` | string | `""` | Extra flags for `gt ltrharvest` (e.g., `-minlenltr 100 -maxlenltr 1000`). |
| `ltrdigest_optional_parameters` | string | `""` | Extra flags for `gt ltrdigest`. |

## `hotspot`

Per-genome detection of windows enriched for ERV integrations beyond chance. A single **Negative-Binomial GLM** models per-window hit counts with a mask-aware offset (`log(effective_bp)`) and an optional per-chromosome baseline; per-window upper-tail p-values are BH-adjusted, thresholded, and merged into hotspot regions. The model is **deterministic** - identical inputs give identical results - so the run is reproducible by construction; the global `parameters.seed` is set and recorded in the manifest for provenance only (the core result does not draw on it). Outputs: per-window `{genome}.csv`/`.parquet`, a `{genome}.manifest.yaml`, merged-region `{genome}.gff3`/`.bed` tracks, and Manhattan/karyotype/Q-Q/summary PDFs.

| Key | Type | Default | Meaning |
|---|---|---|---|
| `hotspot.input` | `catalog` \| `original` | `catalog` | Which event set feeds detection (ADR-012). `catalog` counts **one row per integration event** from the authoritative non-overlapping per-locus assembly, carrying `taxon_call` / `segment` / `structure_class` / `source`; this is the defensible unit for an integration hotspot, and it makes hotspot a consumer of the classification stage (a stale `catalog.csv` pulls `--classify` in). `original` counts raw tBLASTn hits: several features per multi-gene provirus, so counts are weighted by gene content, but the extra density is what gives the NB model power on sparse assemblies. The former `valid` tier is **retired** - it was the same events as `catalog`, multi-counted. **Strict enum.** |
| `hotspot.group_by` | `segment` \| `taxon_call` \| `none` | `segment` | Axis for the per-group NB models. `segment` / `taxon_call` are the calibrated taxonomic calls (ADR-008/011), so grouping is rank-agnostic and follows `classification.segment_rank`; `none` pools every event into one model. Replaces the old `group_split` bool, which split by the legacy probe-derived `label` the redesign superseded (`group_split: true` becomes `group_by: segment`). When the chosen column is absent - notably under `input: original`, where raw hits carry no taxonomic call - the run logs a warning and pools rather than failing. **Strict enum.** |
| `hotspot.source` | `ltr-flanked` \| `orphan` \| `both` | `ltr-flanked` | Which catalog tier is counted. `ltr-flanked` (default) restricts detection to LTR-confirmed proviruses, so a hotspot means 'where do complete integrations cluster'. `orphan` asks the same of the recovered non-LTR tier, and `both` counts the whole catalog - denser, which helps NB convergence, at the cost of mixing LTR-confirmed proviruses with single-gene orphan fragments. Ignored under `input: original`. **Strict enum.** |
| `hotspot.window_size` | int >= 1 | `500000` | Tile width (bp). Sets resolution and the number of windows tested. Defaults to 500 kb because finer windows over a multi-Gb genome are ~99% empty, which collapses the NB dispersion (`theta`) and destroys power; 10 kb yields zero calls on real assemblies. |
| `hotspot.mask_size` | int >= 0 | `20` | Length of the all-N run treated as an unsequenceable gap and subtracted from a window's callable `effective_bp`. `0` disables masking. |
| `hotspot.mask_mismatch` | int >= 0 | `3` | Non-N bases tolerated inside the N-mask motif. |
| `hotspot.pvalue_threshold` | number 0-1 | `0.05` | BH-adjusted q-value cutoff for calling a window significant. |
| `hotspot.min_hits` | int >= 0 | `2` | Minimum total hits in a merged hotspot region (applied after merging). |
| `hotspot.merge_gap` | int >= -1 | `0` | bp gap allowed when merging adjacent significant windows. `0` = merge strictly adjacent windows; `>0` = bridge gaps up to that size; `-1` = no merging (each significant window stays its own region). |
| `hotspot.strata_by_chromosome` | bool | `true` | Include chromosome as an NB covariate so each chromosome gets its own baseline rate. The right value is **tier-dependent**: defaults `true` to match the sparse default `catalog` tier, where the covariate is needed for the NB to converge (without it, fits fail on low-count genomes). **Set `false` when switching to a dense tier like `original`** on fragmented scaffold-level assemblies, where per-contig baselines instead absorb local enrichment and drain power. |
| `hotspot.unplaced_min_factor` | int >= 1 | `10` | Scaffolds shorter than `this x window_size` are pooled into a single `Unplaced` stratum (avoids unstable per-scaffold coefficients). |

## `parameters`

### Thresholds

| Key | Type | Default | Meaning |
|---|---|---|---|
| `seed` | int >= 0 | `67` | Master RNG seed for reproducibility. Recorded in the run manifest for provenance and consumed by any stochastic step. (Hotspot detection is deterministic and does not draw on it.) |
| `identity_threshold` | int >= 0 | `0` | Minimum % identity for BLAST hits. `0` disables the filter. |
| `bitscore_threshold` | number >= 0 | `0` | Minimum bit score. `0` disables the filter. |
| `ltr_resize` | int >= 0 | `0` | Padding (bp) added to each LTR retrotransposon on both sides before overlap detection. |
| `ltr_flank_margin` | int >= 0 | `0` | Tolerance (bp) used when classifying flanking LTRs as left vs right. |
| `merge_option` | `virus` \| `label` | `virus` | How overlapping ranges group before `plyranges::reduce_ranges_directed`. **Strict enum** - typos fail validation. |
| `main_probes` | **ordered** list of strings | `[POL, GAG, ENV, PRO]` | Probe names treated as *main* (as opposed to *accessory*). **The order is meaningful and does four jobs at once** - see "How `main_probes` is used" below. Drives the `probe_type` column on plot dataframes and the `probe_category` attribute on GFF3 tracks. Duplicates are ignored, but unlike a set the sequence is read, so reordering this list changes published columns. |
| `probe_min_length` | map (string -> int) | `{ GAG: 200, POL: 400, ... }` | Per-probe minimum alignment length in residues. Ranges shorter than the probe-specific threshold are filtered out. |

#### How `main_probes` is used

One ordered list drives four separate things. They are listed here because they
do not all want the same order, and the list is a single knob: whichever order
you set governs all four.

1. **Main vs accessory membership** (order irrelevant). The set alone decides
   `probe_type` / `probe_category`.
2. **The mosaic gene set** (order irrelevant). `is_mosaic` is evaluated over
   confident calls on these genes only, so accessory probes such as `P15E` or
   `PR160` can never manufacture a false recombination signal.
3. **Gene reliability ranking** (order used, first = most reliable). When a
   locus's genes disagree, the locus call is taken from the highest-ranked gene
   - after preferring `placement` over `lca`, and before confidence breaks
   remaining ties. `POL` first is the usual choice: it is the most conserved
   marker and the only one with a reliable placement tree.
4. **Expected canonical gene order** (order used). `canonical_order` (and
   `is_canonical` on the ERV-like table) is `true` when a locus's genes, sorted
   by genomic start, read as this list or its exact reverse. The reverse is
   accepted because a minus-strand provirus reads backwards.

**The tension is between (3) and (4).** Reliability wants `POL` first; the
retroviral genome is `5'-gag-pro-pol-env-3'`. Setting `[POL, GAG, ENV]` gives
the best taxon calls, and makes `canonical_order` report `false` for a
structurally textbook `gag -> pol -> env` provirus, because that order is
neither the list nor its reverse. Setting `[GAG, PRO, POL, ENV]` makes
`canonical_order` biologically literal and demotes `POL` in tie-breaking.

Loci carrying two or fewer main genes are unaffected either way: any two-element
order matches either the list or its reverse.

Pick the order for the column you intend to read, and note the choice alongside
the results. `completeness` (and therefore `structure_full_min` and
`structure_class`) uses only the list's **length** as denominator, so it is
insensitive to order.

### Pair detection

| Key | Type | Default | Meaning |
|---|---|---|---|
| `probe_to_pair` | string | `"ENV"` | Probe name used as the anchor for `pair_detector.R`. |
| `pair_max_gap` | number >= 0 | `300000` | Maximum bp distance between paired probes. |

### Aggregation strategies

When overlapping ranges are collapsed via `plyranges::reduce_ranges_directed`, their metadata (virus, label, probe, species) must be reduced to a single value *per merged range*. The strategy is configurable per field.

#### Vocabulary

| Strategy | Output shape | Semantics |
|---|---|---|
| `list` | native multi-value (`CharacterList` in R; `list<string>` in parquet; comma-separated in GFF3 attributes) | Every unique contributor preserved - no information lost. Inflates plot row counts (entry explosion); the plot scripts emit a warning when `virus`/`label` use `list` or `concatenate`. |
| `concatenate` | single string joined by `concat_separator` | Pre-refactor behaviour. Human-readable in a text editor. Same plot-inflation caveat as `list`. |
| `best` | single value | Row with the highest `best_tiebreaker` wins. **Default** for `virus`/`label`. Tie-break is deterministic - see note below. |
| `majority` | single value | Mode (most frequent value). Ties broken by first alphabetical. |
| `first` | single value | Alphabetical first unique value. Deterministic but arbitrary. |
| `strict` | single value or `strict_marker` | Returns the value if all contributors agree, else the configured marker. |

#### `parameters.aggregation`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `aggregation.virus` | strategy | `best` | Strategy for `virus` column. |
| `aggregation.label` | strategy | `best` | Strategy for `label` column. |
| `aggregation.probe` | strategy | `list` | Strategy for `probe` column. Always a reduction grouping key, so single-valued in practice - strategy is effectively a no-op. |
| `aggregation.species` | strategy | `first` | Per-genome pipeline - all rows share one species. |
| `aggregation.best_tiebreaker` | `bitscore` \| `identity` \| `align_length` | `bitscore` | Column used to rank contributors when strategy is `best`. |
| `aggregation.concat_separator` | string | `"; "` | Separator used by `concatenate`. |
| `aggregation.strict_marker` | string | `"ambiguous"` | Value emitted by `strict` when contributors disagree. |

**Deterministic `best` tie-break.** When several contributors to a merged
range tie on `best_tiebreaker`, the winner is resolved by a fixed key chain
applied as a sort of the reduction input: `bitscore` desc -> `query_coverage`
desc -> `identity` desc -> `evalue` asc -> genomic position -> `label` name. This
makes `best` reproducible across R / plyranges versions.

## `solo_ltr`

Native solo-LTR detection (ADR-017). A solo LTR is the single LTR left behind when a provirus's two LTRs recombine homologously and excise everything between them; each one marks an ancestral integration whose provirus is gone. The detector uses the LTR arms of ERV-bearing elements as `blastn` bait and then subtracts the hits that are something else: hits overlapping an intact LTRharvest element, and hits sitting close to an orphan locus (a monoLTR beside surviving coding sequence, so the provirus is damaged rather than excised). See [ADR-017](adr/ADR-017-native-solo-ltr-detection.md) for the method, its calibration against a length-matched random-window null, and what it deliberately gives up.

Every threshold the method depends on is here, and the scripts take them as required arguments with no defaults of their own, so this table is the single source of truth.

| Key | Type | Default | Meaning |
|---|---|---|---|
| `solo_ltr.min_bait_length` | int >= 1 | `300` | Minimum length (bp) of a bait LTR arm. Arms run 100 to 1000 bp, and 80% coverage of a 102 bp arm is not evidence of anything, so short arms are dropped rather than trusted. |
| `solo_ltr.min_hit_length` | int >= 1 | `300` | Minimum alignment length (bp) of an accepted hit. Together with `min_bait_length` this is the near-full-length requirement that took the solo/intact ratio from 492:1 to 27.6:1 on Desmodus. |
| `solo_ltr.min_identity` | number 0-100 | `95` | Minimum percent identity between bait arm and hit. This is an **age filter**: it keeps solos whose family still has a close modern relative, and under-reports ancient ones. |
| `solo_ltr.min_coverage` | number >= 0 | `0.8` | Lower bound on alignment length divided by bait arm length (Ou and Jiang's published rule). |
| `solo_ltr.max_coverage` | number >= 0 | `1.2` | Upper bound on the same ratio, so an alignment running far past the bait is rejected. |
| `solo_ltr.min_alignment_length` | int >= 1 | `80` | Absolute minimum alignment length (bp), also from Ou and Jiang. Redundant with `min_hit_length` at default settings, kept because it is a separate published criterion. |
| `solo_ltr.orphan_pad` | int >= 0 | `10000` | Distance (bp) within which a candidate counts as a **monoLTR at an orphan** rather than a solo. A proviral distance: close enough that the surviving coding sequence plausibly belongs to the same element. |
| `solo_ltr.merge_gap` | int >= 0 | `0` | Gap tolerated when merging overlapping hits into one candidate locus. `0` merges only touching or overlapping hits. |
| `solo_ltr.blast_task` | str (`blastn` \| `megablast` \| `dc-megablast` \| `blastn-short`) | `dc-megablast` | `blastn -task`. Discontiguous megablast is the sensitive-but-affordable setting for diverged nucleotide matches. |
| `solo_ltr.blast_evalue` | number >= 0 | `1.0e-5` | `blastn -evalue`. Generous on purpose: the acceptance criteria above do the real filtering, not the E-value. |
| `solo_ltr.max_target_seqs` | int >= 1 | `5000` | `blastn -max_target_seqs`. High because one LTR family can have thousands of genomic copies and truncating the hit list would silently lose solos. |
| `solo_ltr.blast_threads` | int >= 1 or `null` | `null` | Threads per solo-LTR `blastn` job. `null` gives each job every core, so genomes are searched one at a time. With many genomes, set it to a fraction of the cores (for example 4 on a 64-core machine) so several genomes run in parallel; `blastn` scales poorly past a few threads. Results do not depend on it. |
| `solo_ltr.calibration_lengths` | list of int >= 0 | `[0, 300, 400, 500]` | Length thresholds swept by the per-genome calibration plot. |
| `solo_ltr.calibration_identities` | list of number 0-100 | `[70, 80, 85, 90, 95, 97, 99]` | Identity thresholds swept by the same plot. |
| `solo_ltr.tree.enable` | bool | `true` | Build the LTR phylogeny for this stage. |
| `solo_ltr.tree.n_element_tips` | int >= 0 | `300` | Elements sampled as flanking-arm tips, both arms of each kept. Capped because the clustering statistic compares against class abundance and saturates when one class dominates: using every arm made the Mus tree 96.5% flanking and collapsed its enrichment to 1.03x. |
| `solo_ltr.tree.family_max_distance` | number >= 0 | `0.2` | LTR families are cut from the tree as maximal clades whose largest tip-to-tip distance (substitutions/site) is at most this. 0.2 is the 80-80-80 transposable-element family convention, and where family counts stop tracking the cut on the model 5; below about 0.1 the tree shatters into pairs and the number of families without an intact member becomes an artefact. |
| `solo_ltr.tree.family_panels_per_kind` | int >= 0 | `3` | Families of each kind (no intact member, and with one) drawn as their own subtree, largest by solo count. |
| `solo_ltr.tree.n_solo_tips` | int >= 0 | `200` | Solos sampled as tree tips. A tree over every solo would be neither computable nor readable, so the sample is seeded and reported. |
| `solo_ltr.tree.n_mono_tips` | int >= 0 | `200` | MonoLTR-at-orphan candidates sampled as tree tips. |
| `solo_ltr.tree.permutations` | int >= 0 | `20` | Label permutations for the same-class-sister null, which is what makes "solos cluster with solos" a measurement rather than an impression. |
| `solo_ltr.tree.seed` | int >= 0 | `67` | RNG seed for tip sampling and permutations, so the tree and its statistics are reproducible. |
| `solo_ltr.tree.model` | str | `GTR+G` | IQ-TREE substitution model. |
| `solo_ltr.tree.fast` | bool | `true` | Pass IQ-TREE `-fast` (two iterations, no bootstrap). The tree is **exploratory evidence** and is never an input to classification, so speed beats support values here. |

## `placement`

Phylogenetic-placement figures built from the evidence EPA-ng already produces.
See [ADR-014](adr/ADR-014-publishing-placement-evidence-and-cophylogeny.md).

| Key | Type | Default | Meaning |
|---|---|---|---|
| `mass_norm` | str (`absolute` \| `relative`) | `absolute` | Colour scale for the per-genome heat-trees. `absolute` keeps raw placement mass, so a genome with more ERVs reads hotter; `relative` normalises within each sample, which is what you want when comparing genomes of very different ERV load. |

## `classification`

Per-locus ERV taxonomic classification - turns each valid LTR-element locus into a calibrated **taxon call** (`taxon_call` + `rank` + confidence + mosaic flag + ERV class) from the locus's own marker sequence, instead of transferring the best-bitscore probe label. The classification is **rank-agnostic** (ADR-008): the *axis* - the taxa a locus can resolve to - is declared (see `reference_taxa`), at whatever rank, so a locus resolves to that rank when its evidence lands on an axis taxon, or backs off to an honest higher rank (`rank`) otherwise. Each gene is classified independently against a pinned, taxon-comprehensive reference: POL/GAG/ENV by phylogenetic placement (MAFFT -> EPA-ng -> gappa) when a tree resolves an axis taxon, weighted-LCA otherwise, presence-diagnostic genes (e.g. REX/TAX) by presence. The per-gene calls are then combined into a locus call and a mosaic composition. The reference is built once by the `taxonomy_reference*` rules (`make reference`); see [`docs/taxonomy_classification/`](taxonomy_classification/) and the ADRs for the design. Reuses `parameters.seed` (placement/tree determinism), `parameters.main_probes` (gene reliability order + mosaic gene set + expected canonical gene order - see [How `main_probes` is used](#how-main_probes-is-used)), and `execution.entrez_email` (reference build).

| Key | Type | Default | Meaning |
|---|---|---|---|
| `enable` | bool | `true` | Master switch for the classification stage. When `false`, the `taxonomy_classify` target produces nothing and the pipeline keeps the legacy probe-label provenance only. |
| `placement_genes` | list of str | `[POL, GAG, ENV]` | Genes classified by phylogenetic placement onto a per-gene reference tree; every other gene uses weighted-LCA. Genes listed here must have a built tree package under `data/taxonomy_reference/trees/` (`--build-reference`). **Two measurements disagree about GAG and ENV, and the second is the one that matters.** Reference-alignment identity flags both as weak - the builder reports POL 25.9% (OK), GAG 19.7% and ENV 19.1% (LOW), and all-pairs medians of 23.2 / 16.7 / 16.2 confirm POL is genuinely ~7 points less divergent. But bootstrap support on the resulting trees is comparable or better: median UFBoot 100 / 100 / 99, nodes >=95 at 65.6% / 73.8% / 69.5%, and nodes below 70 at 9.8% / 12.3% / **5.1%**. More divergence yields more informative sites, so a lower-identity alignment can still resolve deep splits; the identity flag is a heuristic proxy, not a verdict on the tree. The practical case for including them is coverage: **~29% of catalog loci carry GAG or ENV but no POL**, and with `[POL]` alone those can never be placed at all. Treat GAG/ENV placements as secondary evidence: low identity leaves them more exposed to systematic error (e.g. long-branch attraction) than bootstrap support alone reveals. Note the ranking consequence - the locus call prefers `placement` over `lca` *before* it applies `parameters.main_probes` order, so a GAG or ENV **placement** now outranks a POL **weighted-LCA** call, where previously POL always spoke for the locus. |
| `search` | str (`blastx`) | `blastx` | Translated-search engine mapping each locus marker region to reference proteins. `blastx` reuses the BLAST+ already in the env (no extra dependency). |
| `evalue` | number >= 0 | `0.001` | blastx e-value cutoff for marker -> reference hits. |
| `top_percent` | number 0-1 | `0.1` | Weighted-LCA bitscore band: hits within this fraction of the best bitscore per marker vote on the lowest-common-ancestor call. Smaller = stricter (fewer, higher-confidence ancestors). |
| `min_orf` | int >= 0 | `30` | Minimum translated marker length (amino acids) for a region to be eligible for phylogenetic placement; shorter markers fall back to weighted-LCA. |
| `confidence_min` | number 0-1 | `0.5` | Confidence floor for the high/low confidence tag. A locus whose call confidence is **below** this value is tagged `LC` (low confidence) in the `confidence_tag` column of the loci/orphans tables; at or above it is `HC`. The threshold is inclusive (`conf == confidence_min` => `HC`) and applies to every method (placement, weighted-LCA, presence). Raise it to flag more marginal calls. |
| `structure_full_min` | number 0-1 | `1.0` | Minimum gene completeness (fraction of `main_probes` present) for a locus to be catalogued as a **`full`** ERV in the `structure_class` column. `1.0` requires every main gene. A single-main-gene locus is always `gene`; a multi-gene locus present but below this floor is `partial`. Deliberately gene-content only - flanking-LTR evidence stays in the anchoring axis (`source`) and the solo-LTR module, not here. |
| `segment_rank` | str | `genus` | Taxonomic rank each `taxon_call` is rolled up to for the by-segment stage (ADR-011). Any NCBI rank (`genus`, `subfamily`, `family`, ...) - the roll-up walks the reference `taxonomy.tsv` hierarchy, so no taxon name is ever hard-coded and the pipeline stays rank-agnostic. A locus whose call is *coarser* than this rank (e.g. `Retroviridae` when segmenting by genus) becomes `unassigned_at_<rank>` rather than being given precision its evidence does not support. |
| `reference_taxa` | list of str | `[]` | The **classification axis** (ADR-008): the taxa - at **any** rank (genus `Lentivirus`, family `Bornaviridae`, ...) - the reference is built at and that a locus can resolve to as a first-class `taxon_call`. Empty or absent derives the axis from the distinct probeset `Label` values, so `Label` seeds the classifier; setting an explicit list decouples the classifier from the probeset. Changing it requires rebuilding the reference (`make reference`). |

## `logging`

`level_styles` and `field_styles` are passed through to `coloredlogs`. See `coloredlogs.install()` documentation for accepted style dicts. Keys: `color`, `bold`, `background`.

## `plots`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `segment_panel` | `full` \| `curated` \| `none` | `full` | How many pages `--segment` draws in each segment's PDF, `results/plots/classification/segments/by_<rank>/<segment>.pdf`. `full` draws every taxonomy and structure page that carries within-segment signal (27). Three of the 30 are excluded automatically because they are degenerate for a single segment: `erv_class_composition` (ERV class is a function of genus, so it is constant within one), and `taxon_confidence_tree` / `taxon_tier_tree` (they draw the taxonomy cladogram, which collapses to a single tip). `curated` draws the 3-page legacy subset: `taxon_composition`, `confidence_gradient`, `structure_class_composition`. `none` writes the per-segment tables and skips the figures. Cost scales as segments x pages and every page has a row per genome, so `curated` is the escape hatch for a large study. **Strict enum.** |
| `bitscore_x_scale` | str (`linear` \| `log10`) | `linear` | X-axis scale on density / raincloud bitscore plots. Use `log10` when the long-tail of low-bitscore hits crushes the lower modes. |
| `sankey_top_n` | int >= 1 or `null` | `null` | Long-tail handling for Sankey plots. `null` (default) shows every stratum. Set to a positive integer N to keep the top N strata per axis and fold the rest into a single labelled `Other (k)` stratum recording how many strata were collapsed. |
| `sankey_other_label` | str | `Other` | Label prefix for the bundled-tail stratum. The actual rendered label is `<prefix> (k)` where `k` is the number of folded strata. |
| `waffle_unit_hits` | int >= 1 | `1` | Number of ranges represented by one waffle square. Bump on huge inputs (e.g. `10` -> "1 square = 10 ranges"). |
| `circle_plot_bitscore_threshold` | number >= 0 | `0` | Bit-score cutoff for circle-plot display. The circle-plot stage is currently broken (see README). |
| `per_stratum` | number >= 0 | `0.18` | Inches added to a stage PDF's page height per genome past 20, so each genome's row keeps its room in a large study (pages are A4 landscape otherwise; see [visual_style.md](visual_style.md)). `0` keeps every page A4. |

## `execution`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `num_cores` | int >= 1 | `64` | Snakemake parallelism cap. |
| `use_species_dict` | bool | `true` | If `true`, `SPECIES` comes from the `species:` map below. If `false`, derive by scanning `SPECIES_DB` for `*.fa` files. |
| `retrieval_time_lag` | number >= 0 | `0.3` | Seconds between Entrez API calls. Required > 0 for NCBI politeness. |
| `max_retrieval_attempts` | int >= 1 | `9` | Entrez retry count. |
| `max_threadpool_workers` | int >= 1 or `null` | `null` | Python `ThreadPoolExecutor` size. `null` = CPU default. |
| `entrez_email` | string | - | **Required.** Contact email for NCBI Entrez API. NCBI ToS requirement. |

## `input`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `probe_csv` | string | - | **Required.** Path to the probe metadata CSV; relative paths resolve against the repo root. Columns expected: `Label, Name, Abbreviation, Probe, Accession`. |
| `species_tree` | str | `''` (none) | Path to a Newick file of the host-species phylogeny (ADR-011). Used to order and annotate the species panels; empty means no species tree and those plots render an explanatory placeholder instead. Tip labels are matched to the `species:` display names, ignoring case and `_` vs space, and any species not found in the tree (or tip not found in the study) is reported in the log rather than silently dropped. Pin the file in your repo/data dir for reproducibility - e.g. a dated TimeTree export. |

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

## `input.pfam_domain_classes`

Path to the curated Pfam domain class table. Default
`data/config/pfam_domain_classes.tsv`.

Three columns, one row per Pfam family:

```
pfam_acc    pfam_name       class
PF00665     rve             retroviral_diagnostic
PF00078     RVT_1           retroelement_shared
PF02994     Transposase_22  non_ltr
PF03184     DDE_1           dna_transposon
PF00098     zf-CCHC         other
```

`pfam_acc` is the key and is unversioned (`PF00665`, not `PF00665.33`). Pfam
guarantees accession stability but reserves the right to rename families, so an
accession-keyed table survives a Pfam upgrade where a name-keyed one would not.
`pfam_name` is kept because `gt ltrdigest` records only names.

The five classes, ordered most to least informative about retroviral identity:

| class | meaning |
|---|---|
| `retroviral_diagnostic` | essentially exclusive to retroviruses and close relatives (`rve`, `RVP`, `TLV_coat`, `Gag_p24`, `GP41`) |
| `retroelement_shared` | reverse-transcribing element, but shared with LINEs and others (`RVT_1`, `RNase_H`) |
| `non_ltr` | LINE/L1 machinery, evidence *against* a retroviral origin (`Transposase_22`, `ORF2p_C`) |
| `dna_transposon` | cut-and-paste transposon (`DDE_1`, `Dimer_Tnp_hAT`) |
| `other` | a real domain, but host housekeeping; says nothing about element type |

The first two classes make a locus `domain_selected`; any other class makes it
`domain_unlisted`; no domain at all is `non_domain`. A family absent from the
table is treated as `other` and is still reported in `domain_names`, never
dropped.

This table replaces the former `domains` block, a map of probe name to
case-insensitive name regexes. That mechanism was retired in ADR-015: measured
across the model genomes it missed 43.8% of the retroviral-diagnostic signal
(`rve`, `RVP`, `IN_DBD_C`, `MLVIN_C` and `GP41` matched no pattern and were
discarded) while its POL pattern `ase` captured `Transposase_22`, an L1 ORF1p
domain, 29,081 times.

**Where it applies.** Exactly one place: the domain scan, which reads Pfam
accessions and works at locus grain, for both tiers (ADR-016). The `ranges` stage
counts LTRdigest's domains (`n_domains_total`, `element_domains`) but deliberately
does not classify them, so there is one `domain_tier` in the project and it means
one thing.

Regenerating the subset is automatic: `pfam_subset_builder` reruns whenever this
table or `Pfam-A.hmm` changes, and never per genome or per run.

## `species`

Map of genome ID -> scientific name. Keys drive `expand(..., genome=SPECIES)` targets in the Snakefile when `execution.use_species_dict: true`.

```yaml
species:
  example_genome_1: "Example species one"
```
