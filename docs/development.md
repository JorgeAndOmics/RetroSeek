# Development

Contributor guide for RetroSeek.

## Environment

```bash
make env              # Create the conda/mamba env from data/config/environment.yml
conda activate retroseek
make env-update       # Update env in place (after pulling changes)
```

The env provides Python 3.10, R 4.3, Snakemake 8, BLAST+, GenomeTools, NCBI Datasets, Bioconductor, and every dev tool (ruff, mypy, pytest, pre-commit, lintr, styler).

## Branching

All work happens on the **`Experimental`** branch. `main` is the stable trunk; do not commit to it directly. Merges from `Experimental` to `main` happen only after explicit review.

```bash
git switch Experimental          # default working branch
git pull                         # stay current
# ... work ...
git push                         # pushes to origin/Experimental
```

## Test-driven development

Three layers:

1. **Python unit** (`pytest`): `tests/unit/` for pure functions in `workflow/scripts/*.py`.
2. **R unit** (`testthat`): `workflow/tests/` for R transforms. Run via `make test-r`.
3. **Snakemake integration**: `tests/fixtures/` holds miniature genomes + config for `snakemake --dry-run` (and full end-to-end runs on tiny data). Run via `make test-snakemake`.

### Cycle

1. Add a failing test that expresses the new behaviour.
2. Run `make test-py` (or `make test-r`) — confirm red.
3. Write the minimum code to turn it green.
4. Run `make check` — all gates must pass.
5. Refactor for SOLID / DRY while keeping tests green.
6. Update [`docs/`](.) and [`README.md`](../README.md) as relevant.
7. Commit.

### Coverage

- Python: target ≥ 80 % line coverage on `workflow/scripts/` (growing from 0 %).
- R: "every non-trivial helper has at least one test."
- Integration: every Snakefile rule exercised at least once by fixtures.

### Markers

- `@pytest.mark.slow` — deselect with `-m 'not slow'`.
- `@pytest.mark.integration` — requires BLAST / GenomeTools.
- `@pytest.mark.network` — requires internet.

## Code style

- **Python**: ruff (lint + format), mypy strict, Google-style docstrings. Configured in [`pyproject.toml`](../pyproject.toml). Target version 3.10.
- **R**: tidyverse style via `styler`, lint via `lintr`. Both installed in the env.
- **Bash**: POSIX-compatible where possible; `shellcheck` if available.

## Commit messages

[Conventional Commits](https://www.conventionalcommits.org/):

```
<type>(<scope>): <summary>

[optional body explaining *why*]
```

Types: `feat`, `fix`, `refactor`, `test`, `docs`, `chore`, `perf`, `style`.

Examples:
- `feat(ranges_analysis): add reduced-range fallback when no domain match`
- `fix(seq_utils): handle empty hit_def gracefully`
- `refactor(defaults): consolidate PATH_DICT construction`

Use the imperative mood. Reference issues where relevant. No co-author trailers.

## Hooks

Minimal hygiene via `pre-commit` (configured in [`.pre-commit-config.yaml`](../.pre-commit-config.yaml)):

```bash
pre-commit install
```

Hooks: trailing whitespace, end-of-file, YAML/JSON validity, mixed line endings, large-file guard. Lint, typecheck, and tests are **not** in hooks — run them via `make check`.

## Before every commit

```bash
make check                          # lint + format-check + typecheck + tests
```

Mentally walk through:

- On `Experimental`? (`git branch --show-current`)
- Tests updated first (TDD)?
- Docs synchronised? ([`docs/architecture.md`](architecture.md), [`docs/usage.md`](usage.md), [`docs/configuration.md`](configuration.md), [`docs/solo_ltr.md`](solo_ltr.md) if the LTR chain was touched, this file, [`docs/adr/`](adr/), [`README.md`](../README.md))
- Conventional Commit message?

## Snakemake patterns to know

A few non-obvious patterns are used in [`workflow/Snakefile`](../workflow/Snakefile) that contributors should understand before editing rules:

- **Checkpoint-gated aggregates.** `blast_pkl2parquet` is a Snakemake `checkpoint`, and the `species_with_hits(wildcards)` function at the top of the Snakefile reads its parquet at DAG-evaluation time. Aggregate rules that previously used a parse-time `SPECIES_POST` constant now use `lambda wildcards: expand(..., genome=species_with_hits(wildcards))`. This makes the DAG a pure function of config + inputs rather than of prior-run disk state. See [ADR-004](adr/ADR-004-species-post-checkpoint.md).

- **Global `wildcard_constraints` on `{genome}`.** Pinned to the exact configured species list via `"|".join(re.escape(s) for s in SPECIES)` near the top of the Snakefile. Prevents the default `.+` regex from greedily absorbing suffixes like `_retroviral` into the wildcard and producing `AmbiguousRuleException` at runtime when two rules write into the same directory with overlapping filename patterns. When adding a new rule whose output shares a directory with another rule's output, prefer unambiguous filename prefixes *and* rely on the constraint — defense in depth.

- **Trap-backed heartbeat on silent long-running rules.** `ltr_index_generator_setup` (suffixerator) and `ltr_harvester_setup` (ltrharvest) wrap their shell commands in a background subshell that emits `[heartbeat:<rule>:<genome>] still running at Nm elapsed` to stderr every 60 s, with `trap '... EXIT'` for cleanup. Adopt the same pattern for any new rule whose primary tool writes all output to stdout (and therefore leaves stderr silent). Do *not* add heartbeats to rules whose tool already emits chatty stderr (e.g. `gt ltrdigest -v`, anything using Python `logging` / `coloredlogs` / `tqdm`) — the cost is zero but the duplication is noise.

## Architectural decision records

Create a new ADR under [`docs/adr/`](adr/) whenever you:

- Pick or replace a tool / library.
- Restructure a significant piece of the pipeline.
- Lock in a default whose reversal someone might later consider.

Use [`docs/adr/ADR-000-template.md`](adr/ADR-000-template.md) as the starting point. Number sequentially.

## Getting help

Open an issue on the project's GitHub repository.
