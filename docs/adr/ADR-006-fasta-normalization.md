# ADR-006: Canonicalise genome FASTA filenames to `.fa` via symlink

- **Status**: Accepted
- **Date**: 2026-04-29
- **Deciders**: Jorge González García

## Context

Genome FASTAs arrive at RetroSeek's `SPECIES_DB` directory with several plausible extensions: `.fa` (legacy), `.fna` (the common case from NCBI Datasets), `.fasta` (manual download convention), and `.ffn` (older GenBank exports). Every downstream rule — `blast_db_generator_setup`, `ltr_index_generator_setup`, `ltr_harvester_setup`, `ltr_retriever_setup` — hard-codes `{genome}.fa` as the input filename. The pipeline broke on hand-supplied non-`.fa` files unless the user manually renamed them first.

A previous workaround lived inside `blast_db_generator_setup` as an inline shell loop:

```
find . -maxdepth 1 -type f \( -iname '*.fasta' -o -iname '*.fna' \) \
| parallel 'mv {} "{.}.fa"'
```

This had three problems: (a) it lived in the wrong rule (BLAST DB generation) for a concern that affects every downstream consumer of the FASTA, (b) it physically renamed files (lossy — the original extension is gone), and (c) it ran inside a rule body that wasn't testable in isolation.

The audit also surfaced a related concern: `defaults.py`'s `SPECIES` auto-discovery scanned `SPECIES_DB` for `.fa` files only, so a fresh clone with only `.fna` genomes saw `SPECIES = []` and produced an empty DAG.

## Decision

Add a new `genome_fasta_normalizer_setup` rule that runs early in the DAG and produces a canonical `{genome}.fa` view of whatever variant exists. The rule is implemented in `workflow/scripts/genome_fasta_normalizer.py`:

- **Strategy**: symlink, not copy. Genomes can be gigabytes; LTR_retriever doesn't modify the input. A symlink is sufficient for every consumer.
- **Extension preference**: `.fa` > `.fna` > `.fasta` > `.ffn`. A real `.fa` (not a symlink) is left untouched. A correct symlink is left in place. A stale symlink is replaced atomically.
- **Ambiguity refusal**: if `.fa` is absent and TWO OR MORE of `.fna`/`.fasta`/`.ffn` coexist for the same genome, the script raises `RuntimeError` listing the candidates. The user must remove duplicates. Silent preference would mask the case where two variants reflect different genomes accidentally co-residing.
- **Cheap shape check**: peek at the first 16 bytes of the source; if the first non-whitespace byte isn't `>`, raise. Catches misnamed binary blobs (gzip, tarballs) before downstream tools spend CPU minutes on them.

Snakefile collision resolution: the normalizer and `genome_downloader_setup` both declare `{genome}.fa` as output. `ruleorder: genome_fasta_normalizer_setup > genome_downloader_setup` pins the normalizer to win. If a user has a local FASTA, no NCBI download fires; if not, the downloader produces `.fa` directly and the normalizer no-ops via its real-file path.

`defaults.py` `SPECIES` discovery widens to scan all four extensions so first-run on a fresh machine still sees genomes before the normalizer materialises the symlinks.

The previous inline `find / parallel mv` workaround is deleted from `blast_db_generator_setup` — the normalizer owns this concern, single source of truth.

## Consequences

- **Positive:**
  - Single source of truth for FASTA extension handling. No more workarounds scattered across rules.
  - Original files are preserved (symlink, not rename). Re-deriving from the source is trivial.
  - Hand-supplied genomes "just work" without manual renaming.
  - Failure is loud: ambiguous directories raise on first attempt; misnamed non-FASTA files are caught at normalization time, not at LTRharvest time.
  - Idempotent — re-runs are safe.

- **Negative:**
  - Don't edit `{genome}.fa` directly when it's a symlink. Edits go to the source `.fna`/`.fasta` instead, which is fine if the user knows; surprising if they don't. Documented in `.claude/memory/gotchas.md` #22.
  - Symlinks can confuse rsync / archive tools that aren't symlink-aware. Generally not an issue inside a single workstation; flag if it surfaces.

- **Neutral:**
  - The script is stdlib-only — no new conda dependencies.
  - The Snakemake rule has no declared input (it scans `SPECIES_DB` at runtime); Snakemake doesn't retrigger on `.fna` mtime changes, but neither did `genome_downloader_setup`, and the script's idempotency makes spurious reruns harmless.

## Alternatives considered

- **Copy, don't symlink.** → Rejected: doubles disk usage on gigabyte-scale genomes; copy time non-trivial.
- **Pick newest mtime when ambiguous.** → Rejected: brittle on systems where mtime is unreliable (`rsync`, archive extraction); silent preference masks user error.
- **Leave the inline rename in `blast_db_generator_setup`, ship a defensive double-write.** → Rejected: would still leave the rename concern owned by an unrelated rule, and DRY violation.
- **Have `genome_downloader_setup` always emit `.fna` and let the normalizer canonicalise.** → Rejected: unnecessary churn — the downloader already produces `.fa` directly via its `unzip -p ... > .fa` path.

## Revisit trigger

- A new common FASTA extension emerges (e.g., `.faa` for amino-acid sequences) — extend `EXT_PREFERENCE`.
- `SPECIES_DB` becomes shared across multiple pipelines that disagree on the canonical extension — at that point we'd version `EXT_PREFERENCE` per pipeline.
- Symlinks become a portability problem (Windows native filesystems pre-WSL, certain archival tools) — switch to copy-on-first-run.

## References

- `workflow/scripts/genome_fasta_normalizer.py` — the implementation.
- `workflow/Snakefile` — `genome_fasta_normalizer_setup` rule + `ruleorder`.
- `workflow/scripts/defaults.py` — widened `SPECIES` discovery.
- `tests/unit/test_genome_fasta_normalizer.py` — extension-preference + ambiguity-refusal + idempotency coverage.
- `.claude/memory/gotchas.md` #22 — operational notes.
