# Changelog

All notable changes to RetroSeek are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Solo-LTR detection now runs** (ADR-013). The stage was fully scaffolded but had
  never produced output. Solo LTRs are the single-LTR remnants of proviruses whose
  flanking LTRs recombined away; in mammals they outnumber intact ERVs by one to
  two orders of magnitude, so this is the catalog's most numerous tier.
  - New `ltr_scn_from_gff3` rule reconstructs LTRharvest's `.scn` from its GFF3,
    which carries every SCN field. Verified byte-identical on *Desmodus rotundus*
    (9,893/9,893 data rows). This unblocks *Antrozous pallidus*, whose SCN is gone
    and whose suffix-array files are all zero bytes, without a ~20 GB index rebuild.
  - Solo LTRs inherit `taxon_call`, `rank`, `segment` and `erv_class` from the
    classified locus their LTR library entry came from - sequence homology rather
    than proximity - because LTR_retriever names library sequences by genomic
    coordinate. A nearest-locus fallback remains for names without coordinates and
    is recorded distinctly in `label_source`.
  - Solos join `catalog.csv` as a third tier (`source=solo-ltr`) under
    `classification.include_solo_ltr`, with `reconcile_catalog` generalised to
    `ltr-flanked > solo-ltr > orphan` precedence.
  - New `ltr_retriever.group_by` for the solo/intact ratio table, following
    ADR-012's grouping vocabulary.

### Fixed

- **Solo LTRs were being read from the wrong file.** The integrator consumed
  `nmtf.pass.list`, which holds *intact* LTR-RTs whose termini lack the canonical
  TGCA motif ("Non-TGCA LTR-RTs" in LTR_retriever's own banner), not solo LTRs.
  Solos now come from `solo_finder.pl` driven off the whole-genome RepeatMasker
  annotation, as LTR_retriever intends.
- **`ltr_retriever.noanno` made solo detection impossible** and is retired.
  `-noanno` suppresses the whole-genome annotation, which produces the only file
  solo detection can read, so the stage could never have found a solo. The runner
  never passes the flag.
- The LTR_retriever pre-filter treated LTRharvest SCN coordinates as 0-based and
  shifted GFF3 starts by `- 1`, widening every valid interval by one base. Both
  formats are 1-based closed, as the byte-identical SCN reconstruction proves. The
  corrected filter reproduces the catalog's locus counts exactly (406 for
  *Desmodus*, 905 for *Antrozous*).
- The pre-filter no longer reads the suffix array's `.des` file, taking the
  `seq-nr` to chromosome mapping from the LTRharvest GFF3 instead. `.des` is zero
  bytes for *Antrozous*.

### Removed

- `parameters.solo_ltr_aggregation` and `ltr_retriever.noanno`. The first
  propagated probe labels onto solo LTRs, a vocabulary ADR-007/008 replaced with
  `taxon_call`; `nearest_erv_max_distance` is renamed `nearest_locus_max_distance`
  for the same reason.

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
