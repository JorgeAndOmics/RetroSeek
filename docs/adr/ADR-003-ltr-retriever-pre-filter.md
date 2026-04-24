# ADR-003: Retroviral-only pre-filter for LTR_retriever

- **Status**: Accepted
- **Date**: 2026-04-24
- **Deciders**: Jorge González García

## Context

RetroSeek needed to add **solo-LTR detection** — finding the single-LTR remnants of ancient retroviral integrations that LTRharvest cannot detect structurally because it only reports paired-LTR candidates. The chosen approach integrates LTR_retriever as a post-processor of LTRharvest output, which internally handles solo-LTR discovery by building consensus LTR sequences from intact ERVs and BLASTing them back against the genome.

However, LTR_retriever's Stage 1 (structural filtering of intact candidates) is retroviral-**agnostic**: it retains any structurally-sound LTR retrotransposon, including non-retroviral families like Copia, Gypsy, and BEL/Pao. If fed LTRharvest output directly, LTR_retriever would build consensus libraries covering all those families, and the solo LTRs it subsequently reports would be dominated by non-retroviral entries.

RetroSeek is specifically about **retroviral** integrations. Its `valid_ranges.gff3` (output of `ranges_analysis_setup`) is the domain-validated retroviral ERV track — the ERVs for which a user-specified probe has matched both structurally (via tBLASTn) and functionally (via Pfam domain match). This set is, by construction, retroviral.

Two ways to integrate this signal with LTR_retriever were considered:

1. **Pre-filter**: intersect LTRharvest SCN with `valid_ranges.gff3` before LTR_retriever sees it.
2. **Post-filter**: run LTR_retriever on unfiltered LTRharvest output, then drop non-retroviral families from its output.

## Decision

Use the **pre-filter** approach. A new script `workflow/scripts/ltr_retriever_prefilter.py` produces `data/ltr_scn/{genome}_retroviral.scn` by intersecting LTRharvest's `.scn` with `valid_ranges.gff3`, and LTR_retriever is invoked on this filtered SCN. (The underscore-separator filename is deliberate: a dot separator — `{genome}.retroviral.scn` — collides with LTRharvest's own `{genome}.scn` output pattern under Snakemake wildcard matching, because `{genome}` would greedily absorb `Antrozous_pallidus.retroviral`.)

Controlled by `config.ltr_retriever.restrict_to_retroviral` (default `true`). Setting `false` reverts to unfiltered behaviour for debugging or non-retroviral exploratory use.

## Consequences

**Positive:**

- **Tighter retroviral-specific consensus clustering.** LTR_retriever's Stage 2 family-building clusters only retroviral LTR sequences. Cluster boundaries reflect retroviral lineage structure rather than being contaminated by Copia/Gypsy consensuses that happen to share identity at the 80% threshold.
- **All downstream output is retroviral by construction.** No post-hoc cleanup — solo LTRs, family-to-ERV mappings, solo/intact ratios all describe retroviral biology only.
- **Cleaner solo-LTR labelling.** The probe-label propagation (Coupling B) operates on a family space that's already retroviral. Nearest-ERV fallback is more reliable because nearby valid ERVs are guaranteed to be retroviral too.
- **Performance win.** Smaller input → faster LTR_retriever runs. On large mammalian genomes, the candidate pool is typically reduced by 10-50×.

**Negative:**

- **Coupling on `ranges_analysis_setup` outputs.** Solo-LTR detection can't run until `valid_ranges.gff3` exists. Snakemake's DAG handles this correctly (it will run ranges_analysis first), but the dependency means solo-LTR detection is always a later-stage target.
- **Non-retroviral solos are invisible by default.** Users who want all LTR-retrotransposon solos must set `restrict_to_retroviral: false` and then filter results themselves. This is the right default for RetroSeek's purpose but would be the wrong default for a general TE-annotation tool.
- **Sensitivity bounded by probe quality.** If RetroSeek's probe set misses a retroviral lineage entirely, LTR_retriever won't find that lineage's solo LTRs either, because no intact ERVs from that lineage make it into `valid_ranges.gff3`. This amplifies the importance of comprehensive probe design.

**Neutral:**

- The pre-filter script is tiny (~250 lines, stdlib-only). Low maintenance burden.
- LTR_retriever's `-inharvest` accepts the filtered SCN without modification — no format translation needed.

## Alternatives considered

### Alternative 1 — Post-filter LTR_retriever's output

Run LTR_retriever on the unfiltered LTRharvest candidates, then drop any family (and its constituent solo LTRs) whose source ERVs don't overlap `valid_ranges.gff3`.

**Rejected because:** LTR_retriever's Stage 2 clustering is one-shot; families form based on the input set. If retroviral and non-retroviral LTRs happen to cluster together at the 80% identity threshold, a post-filter can't unmix them — dropping the family loses retroviral members too. Pre-filtering avoids the problem entirely.

### Alternative 2 — Run LTR_retriever with a retroviral-only reference library

Supply `-u retroviral_te_lib.fa` so LTR_retriever uses a curated external library instead of de novo clustering.

**Rejected because:** This requires maintaining a retroviral LTR library, which is exactly what RetroSeek's probe system already effectively does — at the amino-acid level for domain validation. Re-deriving at the nucleotide level for LTR matching duplicates effort and decouples from the probe-defined retroviral scope.

### Alternative 3 — Forgo LTR_retriever; do DIY BLAST-back

Earlier plan: extract flanking LTRs from RetroSeek's valid ERVs ourselves, build our own consensus via cd-hit, BLAST back against the genome, subtract overlaps with existing ERVs. Simpler in pipeline terms.

**Rejected because:** LTR_retriever's Stage 1 structural filtering (TSD, PBS/PPT, TRF screening) is not trivial to reimplement and meaningfully improves solo-LTR precision. Reinventing it would either produce noisier solos or require adding HMMER/TRF dependencies anyway. LTR_retriever packages the full pipeline, and the pre-filter coupling gets us retroviral specificity on top.

## Revisit trigger

- LTR_retriever's output format changes in a backward-incompatible way that breaks our integrator's parser. Mitigation: the integrator already has a nearest-ERV fallback that continues to produce (weaker) labels in that case.
- A user community emerges that needs non-retroviral LTR solos as a default — we'd flip `restrict_to_retroviral` to `false` and accept the broader output.
- Performance bottleneck on large probe sets — unlikely; the pre-filter is O(candidates × valid_ranges_on_chrom), both usually < 10k.

## References

- `workflow/scripts/ltr_retriever_prefilter.py` — the pre-filter script.
- `workflow/scripts/solo_ltr_integrator.py` — the label-propagation integrator (Coupling B, downstream of this pre-filter).
- `docs/solo_ltr.md` — full biological and algorithmic explainer for the whole solo-LTR workstream.
- `docs/configuration.md` — the `ltr_retriever` config section.
- LTR_retriever: Ou & Jiang 2018, [doi:10.1104/pp.17.01310](https://doi.org/10.1104/pp.17.01310). github.com/oushujun/LTR_retriever.
