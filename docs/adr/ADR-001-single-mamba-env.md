# ADR-001: Single conda/mamba env instead of per-rule `--use-conda`

- **Status**: Accepted
- **Date**: 2026-04-23
- **Deciders**: Jorge González García

## Context

RetroSeek spans Python 3.10, R 4.3, Snakemake 8, external bio tools (BLAST+, GenomeTools, NCBI Datasets CLI), and a large Bioconductor + tidyverse R stack. Snakemake offers two idiomatic environment-management strategies for bioinformatics workflows:

1. One env per rule via `--use-conda`, each with a minimal `environment.yaml` declared in the rule.
2. One env for the whole pipeline, activated before invoking Snakemake.

Both are common. The decision has consequences for installation time, mental load, reproducibility, and flexibility.

## Decision

Use a **single conda/mamba environment** for the whole pipeline, pinned in [`data/config/environment.yml`](../../data/config/environment.yml). `--use-conda` is **not** used.

## Consequences

- **Positive**:
  - One `mamba env create` on install; one activation before running. Low friction for new contributors.
  - Single source of truth for all versions — one file to snapshot for reproducibility.
  - No mamba solve per rule — faster iteration once the env exists.
  - Simpler CI (no nested env management).
- **Negative**:
  - Heavier env (all deps loaded even if a run only touches a subset).
  - If a future rule needs a dependency that conflicts with an existing one, we must either drop the conflict (preferred) or refactor to per-rule envs (see revisit trigger).
  - Less alignment with the bioinformatics-Snakemake norm of per-rule envs.
- **Neutral**:
  - Reproducibility is equivalent (pinning the env fully matches pinning per-rule envs), provided we snapshot `conda list --explicit` alongside runs.

## Alternatives considered

- **Per-rule `--use-conda` envs**: idiomatic for Snakemake, better isolation. Rejected for now: our dependency set has no currently known conflicts, and the single-env model keeps cognitive load lower for a small team.
- **Pixi (unified conda + PyPI lockfile)**: modern tooling, growing adoption. Rejected for now: less battle-tested in bio workflows; can revisit as it matures.
- **Docker / Singularity container**: maximum reproducibility. Rejected for now: heavier setup; slower iteration loop; users on HPC may not have permission to build/run containers.

## Revisit trigger

- A new dependency is added that genuinely conflicts with an existing one in the env.
- The env grows so large that solves or installs become routinely painful (> 15 minutes on modern hardware).
- A user community emerges that needs container-based distribution for publication / HPC policy reasons.

## References

- Snakemake software-deployment docs: https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html
- `data/config/environment.yml` — the pinned env.
