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

All work happens on the **`Experimental`** branch. `master` is the stable reference; do not commit to it directly. Merges from `Experimental` to `master` happen only after explicit review.

```bash
git switch Experimental          # default working branch
git pull                         # stay current
# ... work ...
git push                         # pushes to origin/Experimental
```

Note: if `git-lfs` is not installed locally, pass `-c core.hooksPath=/dev/null` to git operations that trigger LFS hooks (pull, push, checkout).

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
- Docs synchronised? ([`docs/architecture.md`](architecture.md), [`docs/usage.md`](usage.md), this file, [`docs/adr/`](adr/), [`README.md`](../README.md))
- Conventional Commit message?

## Architectural decision records

Create a new ADR under [`docs/adr/`](adr/) whenever you:

- Pick or replace a tool / library.
- Restructure a significant piece of the pipeline.
- Lock in a default whose reversal someone might later consider.

Use [`docs/adr/ADR-000-template.md`](adr/ADR-000-template.md) as the starting point. Number sequentially.

## Getting help

Open an issue or email `jgonzlez@tcd.ie`.
