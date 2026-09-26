# Changelog

All notable changes to RetroSeek are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Code-quality gates**: Qlty (`.qlty/qlty.toml`) measures function complexity
  against a limit of 8, and every function in `workflow/scripts` (Python and R)
  now meets it. `make check` also runs lintr on the R sources (`.lintr.R`, zero
  findings) and ruff's Google-style docstring rules, and `make test-py` enforces
  a coverage floor (68%). Every refactor behind this was checked against the
  model genomes' outputs: identical tables, tracks and rendered figures.
- The hotspot detector has an end-to-end test on a toy genome, replacing a test
  that was always skipped.
- **Phylogenetic placement evidence is now published** (ADR-014). EPA-ng places
  every locus onto the retroviral reference tree on each run; the resulting
  `.jplace` was written into `data/tmp/`, documented as cleared between runs.
  It is the evidence behind every `taxon_call` and the format iTOL and gappa
  read, so it now lands in `results/tracks/taxonomy/placements/`.
- Per-genome **heat-trees** (SVG / Newick / Nexus) showing where a genome's ERV
  load concentrates on the phylogeny, plus **EDPL** placement-uncertainty tables
  - an axis independent of the existing `confidence` column - and LWR
  histograms. New `--placement-trees` stage.
- **Co-phylogeny**: a tree of host genomes built from their ERV placement
  distributions, compared against the host phylogeny by bipartition, with the
  KRD distance matrix behind it. On the model 5 the LTR-flanked tier is
  discordant (RF = 4, 0 of 2 splits shared).
- Tree-ordered composition panels: `species_composition_tree` and
  `taxon_tier_tree`.
- `tree_layout.py` warns when a supplied species tree has uninformative branch
  lengths, so a cladogram is not mistaken for a timetree.

### Changed

- biopython pinned at 1.88 (was 1.87), verified by rerunning every downstream
  stage on the model genomes and comparing the tables: all identical. The one
  visible change is in `{genome}.solos.treefile`, the pruned solo-LTR tree,
  whose branch lengths are now written at full precision (`0.11677255`) where
  1.87's Newick writer rounded to five decimals; no stage reads that file.

### Fixed

- **BLAST hits could overwrite each other.** Each hit was keyed by its accession
  and six random characters; a repeated draw on one chromosome replaced the
  earlier hit without a word (about 2 hits expected lost across the model
  genomes, nearly all in *Mus musculus*). Identifiers are now unique per genome.
- A gappa table whose header lacks `name`, `taxopath` and `aLWR`/`LWR` now stops
  the classification instead of being read from the wrong columns; unreadable
  confidences are reported instead of silently read as 0.
- Unreadable or short rows in a solo-LTR list are reported instead of dropped in
  silence.
- LTR-flanked annotation no longer crashes when a genome has candidate hits but no
  LTR element.
- Probe categories (main, accessory, mixed) are computed for the whole column at
  once: 300 times faster on list-aggregated columns, same values.
- Input validation no longer aborts unattended runs. `validate_ncbi_key` and
  `green_light` called bare `input()`, so a run with no terminal attached (CI, a
  scheduler, `nohup`) died with `EOFError` before Snakemake started. Both now use
  `validator.ask`, which falls back to the default answer. A failed validation
  still refuses to proceed.

- `PATH_DICT["CONFIG_DIR"]` now resolves from the repository root instead of
  `root.data_root_folder`. The two files read from it - `schema.yaml`
  (`validator.py`) and `erv_class.tsv` (the `taxonomy_reference` rule) - ship with
  the code, so pointing the data root outside the repo made validation fail with
  `FileNotFoundError: <data_root>/config/schema.yaml`. Guarded by
  `tests/unit/test_defaults_paths.py`.
- Corrected the `--blast` target in the CLI reference (`docs/usage.md`): the flag
  builds the `blast_pkl2parquet` checkpoint, not `full_genome_blaster`.
- Corrected the conda environment name in the test-fixture docs (`RetroSeek`, not
  `retroseek`), which fails on case-sensitive filesystems.

### Removed

- **Circle plots** (`--generate-circle-plots`, `circle_plot_generator`): the stage
  had not run since April 2026 (it read per-locus scores a refactor replaced) and
  was never brought into the house style. Its config key
  `plots.circle_plot_bitscore_threshold` is retired, and `bioconductor-ggbio`, used
  only by it, leaves the environment.

## [1.1.1] - 2026-05-27

### Added

- GitHub Actions CI (`.github/workflows/ci.yml`): lint, format-check, type-check,
  fast Python tests, and the R `testthat` suite - mirroring `make check`.
- Reproducible anonymized demo-figure generator
  (`workflow/scripts/demo_figures.R`) that rebuilds the README figures from real
  output with neutral placeholder labels.
- Unit tests for species segmentation (`segment_by_probe`) and probe-pair
  detection (`find_pairs`), extracted into pure, sourced modules.

### Changed

- Documented the branching model (short-lived `feat/*` / `fix/*` branches off
  `main`, merged via PR once CI is green), replacing the retired
  `Experimental -> main` flow.
- README quick-start, screenshots (now anonymized demo figures), and signposting
  (CI badge, CHANGELOG link).
- `hotspot_detector` and `circle_plot_generator` are now explicitly marked
  **experimental** (honest `skip()` test scaffolds instead of fake-passing stubs).
- Pinned directly-imported R packages (`scales`, `IRanges`, `S4Vectors`)
  explicitly in `environment.yml` (previously present only transitively).

### Fixed

- Corrected the stale `enERVate` clone URL and `conda activate` env name.
- Replaced personal email addresses and machine-specific `/mnt/v` paths in
  tracked files with placeholders.
- `full_genome_blaster` now serializes an empty table for genomes with zero
  tBLASTn hits (previously `None`, which the converter could not load).

### Removed

- Stale `RetroSeek.yaml` conda export (`environment.yml` is the single env spec).

## [1.1.0] - 2025-06-20

### Added

- Solo-LTR detection via LTR_retriever, pre-filtered to retroviral candidates,
  with probe-label propagation and per-family solo/intact ratios.
- ERV-like composite candidate assembly and a dedicated plotting panel.
- Configurable metadata aggregation across merged ranges (list / concatenate /
  best / majority / first / strict).
- Expanded provirus and stage plotting panels.

## [1.0.1] - 2025-04-07

### Added

- Initial public release: end-to-end Snakemake pipeline for ERV-integration
  detection - genome acquisition (NCBI Datasets), BLAST+ homology search,
  LTRharvest / LTRdigest discovery, R-based range analysis, and plotting.

[Unreleased]: https://github.com/JorgeAndOmics/RetroSeek/compare/v1.1.1...HEAD
[1.1.1]: https://github.com/JorgeAndOmics/RetroSeek/compare/v1.1.0...v1.1.1
[1.1.0]: https://github.com/JorgeAndOmics/RetroSeek/compare/v1.0.1...v1.1.0
[1.0.1]: https://github.com/JorgeAndOmics/RetroSeek/releases/tag/v1.0.1
