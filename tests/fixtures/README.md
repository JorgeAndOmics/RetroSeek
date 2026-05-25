# Test fixtures

Miniature inputs for integration and end-to-end tests. Keep everything here:

- Small enough to run in < 30 s (target: ≤ a few KB per file).
- Self-contained — no external downloads during the test run.
- Deterministic — no randomness, no timestamps, no environment-dependent output.

## Toy genome generator

`build_toy_genomes.py` generates a deterministic, seeded set of synthetic ERV-containing genomes for rapid iteration. Default output: `/mnt/v/databases/toy-genomes/`.

```bash
conda activate enERVate
python tests/fixtures/build_toy_genomes.py              # default: /mnt/v/databases/toy-genomes, seed 1337
python tests/fixtures/build_toy_genomes.py --force      # overwrite existing
python tests/fixtures/build_toy_genomes.py --dest /tmp  # elsewhere
python tests/fixtures/build_toy_genomes.py --seed 42    # different seed
```

### What it produces

- **3 synthetic "species"** — `Toyus_simplex`, `Toyus_complex`, `Toyus_sparsus` — at 30–60 kb × 2–3 scaffolds each.
- **Planted LTR retrotransposons**, 1–3 per scaffold, with direct-repeat flanks that LTRharvest reliably detects.
- **Protein-coding ORFs** back-translated from real retroviral fragments (GAG, POL, ENV, PRO) inside each ERV so tBLASTn produces non-trivial hits.
- **`toy_probes.csv`** — canonical-schema probe CSV (`Label, Name, Abbreviation, Probe, Accession`) pointing at the same proteins.
- **`MANIFEST.tsv`** — per-species size + sha256, for determinism verification.

### Why outside the repo

The genomes live at `/mnt/v/databases/toy-genomes/` — outside the repo, regeneratable anytime from the committed generator script. Keeps the repo small and version-controlled separately from runtime data.

## Conventions

- Synthetic sequences where possible (avoid redistributing large third-party data).
- Fasta IDs should be short and recognisable: `chr1_toy`, `chr2_toy`.
- Name fixture files descriptively so their purpose is obvious in test output.

## Adding a new fixture

1. Add the file under the appropriate subdirectory.
2. Reference it from at least one test.
3. Document any non-obvious property (e.g., "contains one exact LTR match for probe GAG at 500–1500 bp").
