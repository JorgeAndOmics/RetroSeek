# Architectural Decision Records

Short, append-only records of significant architectural decisions. One file per decision, numbered sequentially.

## Format

Use [`ADR-000-template.md`](ADR-000-template.md) as the starting point. Keep each ADR concise (≈ 1 page): context, decision, consequences, alternatives considered.

## Index

| #     | Title                                                                 | Status   |
|-------|-----------------------------------------------------------------------|----------|
| 001   | [Single conda/mamba env instead of per-rule `--use-conda`](ADR-001-single-mamba-env.md) | Accepted |
| 002   | [Configurable metadata aggregation strategies](ADR-002-aggregation-strategies.md) | Accepted |

## When to write a new ADR

- Tool / library selected or replaced.
- Significant refactor of the pipeline structure or data layout.
- A reversible default is chosen that someone might want to revisit.
- A convention is locked in that affects future contributions.

## Status values

- **Proposed** — under discussion, not yet agreed.
- **Accepted** — decision made, implemented.
- **Deprecated** — superseded, replacement noted inline.
- **Superseded** — replaced by a later ADR (reference it).
