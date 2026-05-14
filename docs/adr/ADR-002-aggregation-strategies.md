# ADR-002: Configurable metadata aggregation strategies

- **Status**: Accepted (amended 2026-05-14 — see [Amendment](#amendment-2026-05-14))
- **Date**: 2026-04-24
- **Deciders**: Jorge González García

## Context

`ranges_analysis.R` collapses overlapping tBLASTn hits into merged genomic ranges via `plyranges::reduce_ranges_directed`. Each merged range inherits metadata (virus, label, probe, species) from its contributing hits. The pre-refactor implementation hardcoded one strategy for every metadata field:

```r
reduce_ranges_directed(
  virus = paste(sort(unique(virus)), collapse = "; "),
  label = paste(sort(unique(label)), collapse = "; "),
  ...
)
```

This has three limitations:

1. **Information destroyed.** Multi-value metadata becomes a single `"; "`-joined string. Downstream consumers must re-split to compute counts, modes, or filter by individual labels.
2. **Non-idiomatic for the underlying formats.** GFF3 natively supports multi-value attributes (`Parent=A,B`). Parquet supports `list<string>` columns. Writing semicolon-joined strings to both ignores what those formats are designed for.
3. **Non-configurable.** Biology questions differ — sometimes the caller wants "the strongest single hit per merged range" (max bitscore), sometimes "every contributor" (full list), sometimes "only when contributors agree" (strict consensus). One fixed strategy can't serve all.

The same aggregation question re-appears in the LTR_retriever solo-LTR workstream: how should solo LTRs inherit probe labels from the ERVs that contributed to their consensus? The existing code had no hook for either case.

## Decision

Introduce a small, named vocabulary of aggregation strategies applied per metadata field via config:

- `list` — native multi-value (`CharacterList` in R, `list<string>` in parquet, comma-separated in GFF3). **Default.**
- `concatenate` — single string joined by `concat_separator`. Backward-compatible option.
- `best` — row with the highest `best_tiebreaker` wins. Tiebreaker is `bitscore` (default), `identity`, or `align_length`.
- `majority` — mode.
- `first` — alphabetical first unique value.
- `strict` — single value if unanimous, else `strict_marker`.

Applied independently to `virus`, `label`, `probe`, `species` via `parameters.aggregation` in `config.yaml`. A separate `parameters.solo_ltr_aggregation` block uses the same vocabulary for the LTR_retriever workstream (with a tiebreaker specific to consensus-member counts).

A single `aggregate_values()` helper in `ranges_analysis.R` dispatches on strategy; every `reduce_ranges_directed` call threads it through.

## Consequences

- **Positive:**
  - `list` as default preserves full contributor information. Downstream code can compute any summary (mode, first, best, …) without re-parsing strings.
  - GFF3 and parquet outputs exploit their native multi-value encoding — interop with external tooling (genome browsers, pandas, arrow-based readers) improves.
  - The LTR_retriever workstream and `ranges_analysis.R` share a single vocabulary; one mental model, one helper.
  - Users tune the strategy per biology question without editing R code.
  - The `best` strategy with `bitscore` tiebreaker gives a principled "dominant hit" representation — unlike the previous arithmetic-mean-over-bitscores it replaces.
- **Negative:**
  - Default output shape changes from semicolon-joined strings to multi-value. Downstream scripts that parse semicolon-separated strings (none known in this repo, but potentially in user's ad-hoc analysis) will need updating. Mitigated by `concatenate` strategy as an opt-in backward-compat path.
  - Six strategies is more surface area than one. Documented in `docs/configuration.md` and enforced by yamale schema regex.
- **Neutral:**
  - `reduce_ranges_directed` itself is unchanged — it accepts arbitrary expressions as aggregators, and the helper is called as a named argument.

## Alternatives considered

- **Keep the hardcoded concatenation.** Rejected — fails items 1 and 2 in Context. Downstream code complexity grows with every new consumer.
- **Make each aggregator a user-supplied function (functional style).** More flexible but requires users to write R. The named-strategy approach covers 95% of use cases without code.
- **Support only `list` (the good default) and drop the rest.** Tempting but removes a legitimate backward-compat escape hatch (`concatenate`) and useful single-value strategies (`best`, `strict`). The extra strategies are cheap to implement — one `switch` statement.
- **Use a more ambitious DSL (e.g. a mini-expression language in config).** Overkill for the current need; adds a grammar to maintain.

## Revisit trigger

- A user provides a concrete aggregation need not representable by the six strategies.
- Performance profiling shows the `list` default's `CharacterList` is a bottleneck for very large merged ranges (unlikely at typical bat-genome scale — tens of thousands of hits max per genome).
- The LTR_retriever workstream surfaces a solo-LTR-specific strategy that doesn't fit the shared vocabulary.

## References

- `workflow/scripts/ranges_analysis.R` — `aggregate_values()` helper and consumer sites.
- `data/config/schema.yaml` — `aggregation_schema` and `solo_ltr_aggregation_schema`.
- `docs/configuration.md` — end-user field-by-field reference.
- plyranges documentation — `reduce_ranges_directed` argument evaluation semantics.

## Amendment (2026-05-14)

The original decision shipped `list` as the default for `virus`/`label`. Two
operational findings prompted a default change to **`best`**:

1. **Plot integrity.** `list` (and `concatenate`) produce multi-value cells. The
   plot dataframe pipeline either explodes one locus into N rows (count
   inflation) or carries compound `"A; B; C"` category labels — either way the
   aggregate plots are not statistically meaningful. `best` yields one value per
   merged range, so plots count loci correctly.
2. **Reproducibility.** `best` resolves via `which.max(tiebreaker)`, which is
   order-dependent. With no input sort, exact ties were resolved
   non-deterministically across R / plyranges versions. Shipping `best` as the
   default *required* fixing this.

Changes:
- `parameters.aggregation.{virus,label}` default `list → best` in
  `data/config/config.yaml`.
- `reduce_first` / `reduce_global` (`range_analysis/reductions.R`) now sort
  their input via `sort_for_tiebreak()` before reduction — a fixed key chain
  (`bitscore` desc → `query_coverage` desc → `identity` desc → `evalue` asc →
  genomic position → `label` name). This makes `best` fully deterministic
  without touching the `aggregate_values()` dispatcher.
- The plot scripts (`plot2sort.R`, `stage_plot_generator.R`) emit a logged
  warning and an in-plot caption when `virus`/`label` use `list`/`concatenate`,
  so the inflation caveat is visible on the artifact.

`list`/`concatenate` remain fully supported opt-in strategies — the schema and
the six-strategy vocabulary are unchanged. Only the default moved.
