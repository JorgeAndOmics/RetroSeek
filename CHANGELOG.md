# Changelog

All notable changes to RetroSeek are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **LTR families** (ADR-023). The solo-LTR stage groups its bait arms into
  families with `cd-hit-est`: arms at least `solo_ltr.families.identity` (0.8 by
  default) identical over the whole shorter arm. New tables per genome: every arm
  with its family (`Mmus_F001`...), a family summary (size, genus mix and purity,
  median arm-pair similarity as an age signal) and a family-by-genus table. On the
  model genomes: 820 families, 97% of elements with both arms in one family.
  `cd-hit` joins the environment. Every solo LTR carries its seed arm's family
  (`ltr_family` on the solo table and track), and a new table gives each family's
  solos per intact element.
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

- **Three gene lists instead of one** (ADR-022). `parameters.main_probes` now only
  says which genes are main. Two new settings take over its other jobs:
  `classification.gene_priority` (which gene's call wins the locus call) and
  `parameters.gene_order` (the genes 5' to 3'). A config without them falls back
  to `main_probes` and gets a warning at launch; its taxon calls do not change.
- **`canonical_order` means what it says.** It is checked against `gene_order`
  along the locus's strand and has three answers: `True`, `False`, or blank when
  a locus has fewer than two of the listed genes. With a POL-first `main_probes`
  it used to read `False` for nearly every full provirus. The gene-order page
  now draws only the loci that could be checked.
- The launcher stops on a repeated or lower-case name in the gene lists, and
  warns about a listed name that is not a probe. The ranges manifest records
  `main_probes`.
- The README's demo figures were regenerated from the current model-genome run.
  Retroviral genera now show under their real names and house colours; host
  species and provirus names stay anonymised.
- Environment: filelock 4.0.4 (was 3.32.6), libglib 2.90.0 (2.88.3), virtualenv
  21.13.0 (21.7.9), yq 4.3.0 (4.1.2), and later tqdm 4.70.1, sqlalchemy 2.0.54, idna
  3.20, pyparsing 3.3.3 and platformdirs 4.12.2, all indirect dependencies that a
  fresh `make env` already resolves to. Verified by a full downstream rerun on the model
  genomes against the previous one.
- biopython pinned at 1.88 (was 1.87), verified by rerunning every downstream
  stage on the model genomes and comparing the tables: all identical. The one
  visible change is in `{genome}.solos.treefile`, the pruned solo-LTR tree,
  whose branch lengths are now written at full precision (`0.11677255`) where
  1.87's Newick writer rounded to five decimals; no stage reads that file.
- The ten longest functions (tree views, tree build, hit tabulation, domain scan,
  solo finder and annotator, classifier entry, counts, probe CSV check) are split
  into named steps. Checked old against new: identical outputs on the model
  genomes, or identical tool commands where the tools are slow.

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
- **The launcher no longer renames genome files.** Before every run it renamed
  `.fasta`, `.fna` and `.fas` files in the genome folder to `.fa`, replacing an
  existing `.fa` of the same name. The `genome_fasta_normalizer_setup` rule already
  gives each genome a `{genome}.fa` link without touching the source files, and now
  accepts `.fas` too; the launcher's validation reads the same file that rule links
  to. Extensions are matched in lower case: rename a `.FNA` to `.fna`.
- The example config in `tests/fixtures/` validated no more: it still carried the
  retired `domains:` block. A test now checks both shipped configs against the
  schema.
- Every reader of the pipeline's own GFF3 tracks (solo-LTR finder, bait
  builder, classifier) follows one rule: a damaged row stops the job and names
  its line, and a `##FASTA` section ends the features. Some readers used to skip
  such rows without a word (losing an element, an arm or a locus), one failed
  without saying where.
- The classifier read GFF3 attributes with a pattern that also matched inside a
  longer key (`probe=` inside `subprobe=`). All GFF3 readers now share one
  attribute parser that matches whole keys; the pipeline's own tracks read the
  same as before.
- A BLAST hit with a blank or missing sequence name stops the hit parser with a
  message naming the genome, instead of an IndexError (blank) or an empty
  chromosome name that overlapped nothing downstream (missing).
- A `{genome}.fa` link pointing at a different file from the `.fna`, `.fasta`,
  `.ffn` or `.fas` beside it now stops the launcher's preflight as ambiguous, like
  two such files without a `.fa`; it used to be used as it was. A link with no such
  file beside it (a genome stored elsewhere) is still honoured.
- A tie for a hotspot's dominant lineage is broken in byte order, the same on
  every machine; it followed the machine's locale. The model genomes' lineage
  names sort alike either way, so their tables do not change.
- The hotspot stage's notice about FASTA headers with no name now follows the
  console line contract (`time level step genome | message`) instead of a bare
  R message.
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
- `pfam_name_to_acc.tsv`, a Pfam name to accession table the Pfam subset step
  wrote and no stage read. The step reruns once, with the domain scan after it.

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
