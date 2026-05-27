# Changelog

All notable changes to RetroSeek are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- GitHub Actions CI (`.github/workflows/ci.yml`): lint, format-check, type-check,
  fast Python tests, and the R `testthat` suite — mirroring `make check`.
- Reproducible anonymized demo-figure generator
  (`workflow/scripts/demo_figures.R`) that rebuilds the README figures from real
  output with neutral placeholder labels.
- Unit tests for species segmentation (`segment_by_probe`) and probe-pair
  detection (`find_pairs`), extracted into pure, sourced modules.

### Changed

- Documented the branching model (short-lived `feat/*` / `fix/*` branches off
  `main`, merged via PR once CI is green), replacing the retired
  `Experimental → main` flow.
- README quick-start, screenshots (now anonymized demo figures), and signposting
  (CI badge, CHANGELOG link).
- `hotspot_detector` and `circle_plot_generator` are now explicitly marked
  **experimental** (honest `skip()` test scaffolds instead of fake-passing stubs).

### Fixed

- Corrected the stale `enERVate` clone URL and `conda activate` env name.
- Replaced personal email addresses and machine-specific `/mnt/v` paths in
  tracked files with placeholders.
- `full_genome_blaster` now serializes an empty table for genomes with zero
  tBLASTn hits (previously `None`, which the converter could not load).

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
  detection — genome acquisition (NCBI Datasets), BLAST+ homology search,
  LTRharvest / LTRdigest discovery, R-based range analysis, and plotting.

[Unreleased]: https://github.com/JorgeAndOmics/RetroSeek/compare/v1.1.0...HEAD
[1.1.0]: https://github.com/JorgeAndOmics/RetroSeek/compare/v1.0.1...v1.1.0
[1.0.1]: https://github.com/JorgeAndOmics/RetroSeek/releases/tag/v1.0.1
