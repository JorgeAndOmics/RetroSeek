# Changelog

All notable changes to RetroSeek are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

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

### Fixed

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
