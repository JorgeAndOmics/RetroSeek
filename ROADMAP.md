# RetroSeek roadmap

Feature maturity for RetroSeek's pipeline stages.

- **Stable** — tested and used in production runs.
- **Experimental** — functional but with limited/no automated test coverage; the
  interface or output may change, and edge cases are under-validated.
- **Planned** — not yet implemented.

| Stage / feature | CLI flag | Status | Notes |
|---|---|---|---|
| Genome acquisition + indexing | `--download-genomes`, `--blast-dbs`, `--suffix-arrays` | Stable | |
| LTR discovery (LTRharvest / LTRdigest) | `--ltr-candidates`, `--ltr-domains` | Stable | |
| Probe extraction + tBLASTn | `--probe-extractor`, `--blast` | Stable | |
| Range analysis (integration → GFF3 tracks + tables) | `--ranges-analysis` | Stable | Core R layer; well covered by testthat |
| ERV-like candidate assembly | part of `--ranges-analysis` | Stable | |
| Species segmentation (main / accessory) | part of the BLAST → tables stage | Stable | Unit-tested (`segment_by_probe`) |
| Solo-LTR detection | `--solo-ltr-detection` | Stable | See [`docs/solo_ltr.md`](docs/solo_ltr.md) |
| Provirus / stage plots | `--generate-global-plots` | Stable | |
| Probe-pair detection | `--pair-detection` | Stable | Unit-tested (`find_pairs`) |
| Hotspot detection | `--hotspot-detection` | **Experimental** | regioneR permutation testing — stochastic and slow; not yet unit-tested |
| Per-genome circle plots | `--generate-circle-plots` | **Experimental** | ggbio Circos-style rendering — visual output; not yet unit-tested |

## Planned / under consideration

- **Promote the experimental rules to Stable** by extracting their pure logic
  into testable helpers (hotspot window-tiling; circle-plot karyotype builder),
  mirroring the `species_segmenter` / `pair_detector` refactors.
- **A wildcard-free default target** (`rule all`) so `make test-snakemake`
  (a bare `snakemake -n`) can resolve the DAG without an explicit target.
