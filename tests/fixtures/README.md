# Test fixtures

Miniature inputs for integration and end-to-end tests. Keep everything here:

- Small enough to run in < 30 s (target: ≤ a few KB per file).
- Self-contained — no external downloads during the test run.
- Deterministic — no randomness, no timestamps, no environment-dependent output.

## Suggested layout

```
fixtures/
├── genomes/            # Tiny synthetic or trimmed FASTA files
│   └── toy_genome.fa
├── probes.csv          # Minimal probe metadata CSV
├── config.yaml         # Shrunken config pointing at the files above
└── expected/           # Reference outputs for comparison (optional)
```

## Conventions

- Synthetic sequences where possible (avoid redistributing large third-party data).
- Fasta IDs should be short and recognisable: `chr1_toy`, `chr2_toy`.
- Name fixture files descriptively so their purpose is obvious in test output.

## Adding a new fixture

1. Add the file under the appropriate subdirectory.
2. Reference it from at least one test.
3. Document any non-obvious property (e.g., "contains one exact LTR match for probe GAG at 500–1500 bp").
