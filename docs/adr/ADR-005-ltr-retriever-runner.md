# ADR-005: Wrap LTR_retriever invocation in a Python runner script

- **Status**: Accepted
- **Date**: 2026-04-29
- **Deciders**: Jorge González García

## Context

Calling LTR_retriever from the Snakemake `ltr_retriever_setup` rule used to be an inline `shell:` block:

```
mkdir -p {params.out_dir} && \
cd {params.out_dir} && \
LTR_retriever -genome {input.genome_fa} -inharvest {input.scn} \
              -u {params.u} -miniden {params.iden} -threads {threads} {params.noanno} && \
for ext in pass.list.gff3 nmtf.pass.list LTRlib.fa; do \
    if [ -f {wildcards.genome}.fa.mod.$ext ]; then \
        mv {wildcards.genome}.fa.mod.$ext {wildcards.genome}.$ext ; \
    fi ; \
done
```

Several concerns sat tangled together: working-directory staging, binary invocation, stdout/stderr handling, the `.fa.mod.*` rename, and missing-output detection. The `[ -f X ] && mv X Y` pattern silently no-op'd if a file was missing — a real failure (e.g., a new LTR_retriever release using a different extension) would surface only as a downstream `FileNotFoundError` from `solo_ltr_integrator`. Logs went to Snakemake's default location, not a per-rule path, making remote debugging painful. And the rule body was difficult to test — the only way to exercise it was to actually run LTR_retriever, requiring the conda env, real genome data, and minutes of CPU.

The audit also introduced a new requirement: the prefilter now produces both `{genome}_retroviral.scn` and `{genome}_full.scn`, and the rule needs to pick which one feeds LTR_retriever based on `config.ltr_retriever.source_scn`. Adding that selection to the inline shell would have made it even more complex.

## Decision

Move all of LTR_retriever's invocation logic into `workflow/scripts/run_ltr_retriever.py` and have the Snakemake rule call it. The script decomposes into four building blocks:

- `resolve_source_scn(mode, retroviral, full)` — pure function picking which SCN feeds LTR_retriever based on `config.ltr_retriever.source_scn`. Raises `ValueError` on unknown mode (validator should have caught it; this is a belt-and-braces check at the runner boundary).
- `stage_workdir(workdir, genome_fa, scn, genome_name)` — creates the per-genome workdir under `LTR_RETRIEVER_DIR/{genome}/`, symlinks the genome FA and chosen SCN. Symlinks (not copies) keep gigabyte-scale genomes off the intermediate filesystem; LTR_retriever doesn't modify the input. Stale symlinks are replaced rather than reused.
- `run_binary(staged, params, log_file)` — `subprocess.run` of LTR_retriever with `cwd=workdir`. Tees stdout and stderr to the log file with `[stdout]` / `[stderr]` line prefixes. Returns the binary's exit code.
- `finalise_outputs(workdir, genome_name)` — renames `{genome}.fa.mod.{ext}` → `{genome}.{ext}` for the three expected extensions (`pass.list.gff3`, `nmtf.pass.list`, `LTRlib.fa`). If any expected output is missing after rename, raises `RuntimeError` with the workdir's directory listing.

The Snakemake rule body collapses to one line: `python scripts/run_ltr_retriever.py [args]`. Logs land at `LOG_DIR/ltr_retriever/{genome}.log`.

## Consequences

- **Positive:**
  - Each building block is unit-testable in isolation with `tmp_path` fixtures and a tiny bash shim binary — no real LTR_retriever required for most tests.
  - Failures are loud: missing outputs raise `RuntimeError` with a workdir listing rather than silently no-op'ing.
  - LTR_retriever's stdout + stderr are captured to a per-genome log file with line-level prefixes for easy grep.
  - Per-genome workdir staging keeps `.fa.mod.*` artefacts contained — `SPECIES_DB` is no longer polluted.
  - Symlink-based staging avoids duplicating gigabyte-scale genomes per rule invocation.
  - The `source_scn` enum routing lives in a single small function (`resolve_source_scn`), unit-tested.

- **Negative:**
  - One more script to maintain (~250 lines, stdlib-only).
  - Slight indirection: a contributor reading the Snakefile rule has to also open `run_ltr_retriever.py` to understand the full invocation.

- **Neutral:**
  - The script is stdlib-only — no new conda dependencies.
  - It takes the LTR_retriever binary path as an explicit `--ltr-retriever-binary` arg with a `shutil.which("LTR_retriever")` fallback, so tests can inject shim binaries.

## Alternatives considered

- **Keep the inline shell, harden it.** Add `set -e`, log redirection, explicit existence checks. → rejected: the rename loop is still a multi-line bash construct that resists mocking; testing would require a real LTR_retriever binary.
- **Use Snakemake's `script:` directive with a Python helper.** Closer to convention but requires the script to be Snakemake-aware (use the magic `snakemake` global), which makes it harder to test and harder to invoke standalone for debugging.

## Revisit trigger

- LTR_retriever moves to a stable Python API (e.g., a callable library entry point) — at that point we'd skip the subprocess + workdir staging and call the API directly.
- A new LTR_retriever release changes the output filename pattern away from `<genome>.fa.mod.*` — `EXPECTED_OUTPUT_EXTS` and the rename logic in `finalise_outputs` need a small update; the runner's overall structure is unchanged.

## References

- `workflow/scripts/run_ltr_retriever.py` — the runner.
- `workflow/Snakefile` `rule ltr_retriever_setup` — consumer.
- `tests/unit/test_run_ltr_retriever.py` — coverage of every building block + an end-to-end test with a fake bash shim binary.
- [ADR-003](ADR-003-ltr-retriever-pre-filter.md) — Coupling A pre-filter rationale.
- `.claude/memory/gotchas.md` #20, #21 — operational notes on LTR_retriever behaviour.
