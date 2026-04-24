# ADR-004: `SPECIES_POST` → Snakemake checkpoint + runtime `species_with_hits(wildcards)` resolver

- **Status**: Accepted
- **Date**: 2026-04-22
- **Deciders**: Jorge González García

## Context

Historically, RetroSeek's Snakefile computed a parse-time constant `SPECIES_POST` by reading `full_genome_blast.parquet` off disk:

```python
SPECIES_POST = retrieve_species_list(...)    # evaluated when Snakemake parses the Snakefile
```

Aggregate rules (`ranges_analysis`, `hotspot_detector`, `pair_detector`, plot generators) then expanded over `SPECIES_POST` to limit downstream work to species that actually produced tBLASTn hits:

```python
rule ranges_analysis:
    input:
        expand(..., genome=SPECIES_POST)
```

This introduced a **disk-state dependency in DAG evaluation**:

- On a **fresh run**, `full_genome_blast.parquet` didn't yet exist, so `retrieve_species_list()` fell back to the full `SPECIES` list. The DAG included every configured species' downstream rules — work that would be wasted if some species had zero hits.
- On a **second run** after the BLAST stage had produced the parquet, `SPECIES_POST` reflected the hit subset, and downstream rules correctly scoped to those species.
- On a **third run after probe changes** but without clearing the parquet, `SPECIES_POST` still reflected the stale prior-run list, so the new-probe hit set wouldn't be applied.

In short: **the DAG was a function of the filesystem, not of config + inputs**. Snakemake's reproducibility guarantees assume the opposite.

A secondary issue: the pattern baked the filesystem dependency into Snakefile parsing, which happens before any rule runs. There was no way for the DAG to actually *wait* for the parquet to be produced — it either existed at parse time or it didn't.

## Decision

Convert `blast_pkl2parquet` from a regular rule to a **Snakemake checkpoint**. Introduce a runtime resolver `species_with_hits(wildcards)` that reads the checkpoint's parquet output at DAG-evaluation time. Every aggregate rule that used to use `expand(..., genome=SPECIES_POST)` now uses:

```python
lambda wildcards: expand(..., genome=species_with_hits(wildcards))
```

The `checkpoints.blast_pkl2parquet.get().output.parquet_full` call inside the resolver forces Snakemake to block DAG evaluation of any dependent rule until the checkpoint has executed, at which point the parquet is guaranteed to exist and can be read authoritatively.

One coupled change: `species_segmenter_setup` now declares its per-species outputs via `genome=SPECIES` (the full configured list) rather than `genome=SPECIES_POST`. The R script writes an empty-schema parquet for any configured species without hits, so the declared outputs always materialise. Downstream rules demand only the hit subset via the runtime resolver.

## Consequences

- Positive:
  - **The DAG is now a pure function of config + inputs**, not of prior-run filesystem state. Same config + same inputs → deterministic DAG on first invocation.
  - The hit-filter optimisation (skip downstream work for zero-hit species) now works on first invocation rather than only on re-runs.
  - The parquet is explicitly the source of truth — no silent fallback to "run everything" when the file is absent.
  - Rule graphs are explicit about their dependency on the BLAST checkpoint; no hidden parse-time disk reads.

- Negative:
  - Checkpoints carry a small runtime overhead — Snakemake has to pause the scheduler, evaluate the checkpoint, then re-plan the DAG. Measured overhead on 5-genome mammalian runs: < 5 seconds.
  - Rule `input:` blocks for dependent aggregates are wrapped in a `lambda wildcards:` — slightly more syntactically complex than a bare `expand()`. Contributors must understand the resolver pattern.
  - Debugging "why is rule X not running?" is harder when the DAG is partially dynamic. Use `snakemake --dry-run` which now prints `"The run involves checkpoint jobs, which will result in alteration of the DAG..."` as a hint.

- Neutral:
  - The fallback behaviour is preserved: if the parquet somehow produces zero species (pathological probe set, zero BLAST hits anywhere), `species_with_hits()` returns the full `SPECIES` list defensively. This shouldn't happen in practice but guards against an empty-DAG catastrophe.
  - The historical pre-B1 Snakefile versions remain bisectable — every commit in the B1 series (`557b9bf`, `95dec47`) leaves the pipeline runnable.

## Alternatives considered

- **Option A: Leave `SPECIES_POST` as a parse-time constant, document the quirk**. Rejected because it contradicts reproducibility by making the DAG filesystem-dependent. "Document the footgun" is rarely an acceptable answer when a structural fix is available.
- **Option B: Make every downstream rule expand over the full `SPECIES` list and let individual rule bodies early-exit for no-hit species**. Rejected because it wastes DAG edges and scheduler work, and because "some output files exist but are empty" is a fragile contract — Snakemake's rerun-triggering heuristics get confused.
- **Option C: Introduce a pre-parse shell hook that materialises `SPECIES_POST` as a YAML file, then read it at parse time**. Rejected because it adds an out-of-band execution step outside Snakemake's dependency graph, which defeats the point of Snakemake.
- **Option D: Write our own "dynamic input" mechanism via Snakemake's `params:` block and `run:` directive**. Rejected as reinventing the checkpoint mechanism that Snakemake already provides.

## Revisit trigger

- If Snakemake deprecates or significantly changes the `checkpoint` mechanism in a future release.
- If the pattern proves too cumbersome for new contributors — in which case we should consider extracting the `lambda wildcards:` pattern into a reusable helper function (partial application) rather than reverting.

## References

- Snakemake checkpoint docs: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution
- Landing commits: `557b9bf` (refactor: SPECIES_POST → checkpoint + species_with_hits resolver), `95dec47` (fix: three issues surfaced by toy-genome testing).
- Related: the `species_with_hits(wildcards)` function is defined at the top of `workflow/Snakefile` with a full docstring explaining the semantics.
