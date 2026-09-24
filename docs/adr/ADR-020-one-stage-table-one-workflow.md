# ADR-020: One stage table, one workflow per run, and a guard on the heavy searches

- **Status**: Accepted
- **Date**: 2026-09-23
- **Deciders**: Jorge González García
- **Builds on**: [ADR-019](ADR-019-pinned-pfam-release.md) (the Pfam trap the guard catches)

## Context

The launcher had grown by accretion:

- **Two parsers.** The root `./RetroSeek` and `workflow/scripts/RetroSeek.py` each
  declared every flag, with help texts that had drifted apart.
- **Every run looked like a success.** The root shim dropped the pipeline's exit
  code, so a failed stage or a refused validation still returned 0.
- **One Snakemake per flag, in a fixed order that was not the dependency order**
  (`--placement-trees` and `--segment` before `--classify`). A dry run with several
  flags printed several overlapping job tables, and a forced rule reran once per
  flag.
- **No protection for the heavy searches.** LTRdigest costs about a day per large
  genome, and harmless-looking changes (a redownloaded Pfam file, a reworded rule)
  wake it for every genome. The run protocol said "dry-run first and stop if a heavy
  rule appears", but only a careful human enforced it.
- **The element and orphan domain scan had no flag.** It ran only as a hidden
  dependency of `--classify`, next to the unrelated LTRdigest flag `--ltr-domains`.

## Decision

- **One table of stages**, `workflow/scripts/stages.py`. Each row holds a flag, its
  phase, its targets, its help, the heavy rules it exists to run and the tools it
  calls. The parser (grouped by phase), the run order, the guard and the tool check
  are all built from it. The module imports nothing from the pipeline, so `-h` works
  without a config.
- **A thin root shim.** It resolves `--configfile`, handles `--config-help`, runs
  the pipeline with the same interpreter and returns its exit code.
- **One Snakemake call per run**, with every requested target: one DAG, one dry-run
  table, parallelism across stages. `--keep-going` is the default so one failed
  genome does not stop the rest overnight; `--stop-on-error` turns it off.
- **A guard.** Before a real run the launcher captures a dry run. If a heavy rule
  would run and no requested stage owns it, nothing runs; the message names each
  rule, Snakemake's reason and the usual fix. `--allow-heavy` overrides.
- **Fast checks always, slow checks optional.** The schema, the tools of the chosen
  stages and the Pfam library are checked on every run; `-skp` now skips only the
  NCBI lookups and the prompts.
- **New flags:** `--domain-scan` (the symmetric domain scan, after ranges) and the
  preset `--downstream` (every Analysis and Figures stage except circle plots).

## Consequences

- Positive:
  - A failed run fails the command, so scripts and schedulers can trust `$?`.
  - The heavy searches cannot rerun by accident.
  - `-h` reads as the pipeline's phases, and a dry run shows one table.
  - A new stage is one table row; its help, order, guard and tool check follow.
- Negative:
  - Every real run starts with a dry run (seconds on five genomes, up to about a
    minute on a hundred).
  - The guard reads Snakemake's text output (the "Job stats" table and "reason:"
    lines). A future Snakemake that changes them makes the guard see no jobs; the
    unit tests, built on a real 9.24 dry run, catch that.
  - `-skp` no longer skips the schema check; a config that used to pass unchecked
    may now stop.
- Neutral:
  - Snakemake options still pass through unchanged. `--unlock` and
    `--cleanup-metadata` bypass the guard because they run no jobs.

## Alternatives considered

- **Subcommands** (`./RetroSeek setup | discover | analyse`): self-documenting, but
  changes every existing command, doc and runbook for the same grouping the help
  sections now give.
- **Keep one Snakemake per flag and only regroup the help**: leaves the repeated
  dry-run tables and the out-of-order stages.
- **Protect heavy outputs with Snakemake's `ancient()`**: stops timestamp reruns,
  but not reruns triggered by a changed rule or by an upstream job that reruns.

## Revisit trigger

- A Snakemake release that changes its dry-run text, or offers a structured
  (JSON) dry-run report the guard could read instead.
- A cluster profile that needs per-stage Snakemake calls.

## References

- `workflow/scripts/stages.py`, `workflow/scripts/guard.py`, `RetroSeek`
- [usage.md, CLI](../usage.md#cli)
