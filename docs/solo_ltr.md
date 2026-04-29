# Solo-LTR detection in RetroSeek

This document explains how RetroSeek detects **solo LTRs** — the single-LTR remnants of ancient retroviral integrations that LTRharvest cannot find structurally — by integrating [LTR_retriever](https://github.com/oushujun/LTR_retriever) into the pipeline. It covers the biology, the mechanism, and the coupling with RetroSeek's existing probe-based detection.

If you just want to configure or run it, jump to [Configuration & usage](#configuration--usage). For the full biological and algorithmic detail, start from the top.

---

## Why solo LTRs matter

An intact endogenous retrovirus (ERV) integrates into the host genome as a provirus: two identical **long terminal repeats (LTRs)** flanking an internal coding region (`gag`, `pol`, `env`, etc.). Over evolutionary time — millions of years for ancient families — the two flanking LTRs, being perfect direct repeats, **recombine homologously**. The recombination event excises the entire internal region, leaving behind a **single LTR** at the original integration site.

This is not a degenerate case: in most mammalian genomes, solo LTRs outnumber intact ERVs by **one to two orders of magnitude**. They are silent genomic fossils of retroviral invasions that left no intact provirus behind. Each solo LTR marks exactly one ancestral integration event. Their count per family is a proxy for the **age and abundance of that lineage** — high solo/intact ratios indicate old, well-expanded invasions.

From a RetroSeek perspective, solo LTRs are valuable because:

- They **extend the ERV inventory** — finding integrations our probe-based tBLASTn search would miss entirely, because there's no retroviral protein left at the locus.
- They **enable dating proxies** — the solo/intact ratio per family is a simple, interpretable age signal.
- They **complete hotspot maps** — a hotspot defined only by intact ERVs systematically underrepresents high-activity lineages that have recombined heavily.

## Why LTRharvest alone cannot find them

LTRharvest detects paired-LTR structures: two sufficiently-similar direct repeats flanking an inner region of plausible length. It scans the genome with a suffix-array index and reports candidates with two LTRs. A solo LTR has no pair — only one LTR, surrounded by non-retroviral genomic context. LTRharvest silently ignores it.

No parameter adjustment fixes this. LTRharvest's search primitive **requires** the pair. Finding solo LTRs requires a different approach: start with the paired LTRs LTRharvest *did* find, build a consensus of them, and then search the genome for **single instances** of that consensus that don't overlap any intact ERV. This is what LTR_retriever's pipeline does.

---

## How LTR_retriever works

LTR_retriever ([Ou & Jiang 2018](https://doi.org/10.1104/pp.17.01310)) is a post-processor that consumes LTRharvest output (and optionally LTR_FINDER output) and produces a filtered, classified, and enriched ERV inventory — including solo LTRs. It runs in three conceptual stages.

### Stage 1 — Structural filtering of intact candidates

LTRharvest's reported paired-LTR structures have a high false-positive rate: many are tandem repeats, segmental duplications, or non-LTR elements that happen to resemble the pattern. LTR_retriever filters them with biological criteria:

- **Target Site Duplication (TSD).** Real retroviral integrations leave a 4–6 bp host sequence duplication flanking the LTR pair. LTR_retriever verifies TSD presence and typical length.
- **LTR similarity in biological range.** Too-similar LTRs (>99%) suggest very recent events or segmental duplication; too-divergent (<70%) often aren't related at all.
- **Internal retroviral features.** Primer Binding Site (PBS), PolyPurine Tract (PPT), and minimum internal length are checked via HMMER.
- **Tandem-repeat screening.** TRF (Tandem Repeat Finder) rejects candidates whose internal regions are tandem-repeat-dominant.
- **Nesting detection.** Resolves the case where an ERV integrated inside another ERV.

Surviving candidates are written to `*.pass.list` (plus `*.pass.list.gff3` in GFF3 format). These are the **intact ERVs** used in later stages.

> Note: Stage 1 filters for "real LTR retrotransposons," which includes non-retroviral LTR-retrotransposons like Copia, Gypsy, and BEL/Pao. It does not specifically filter for retroviruses. That's where RetroSeek's coupling kicks in — see [Coupling A](#coupling-a--retroviral-only-pre-filter) below.

### Stage 2 — Family classification via de novo clustering

LTR_retriever then clusters the intact ERVs' flanking-LTR sequences into **families**. A family is a **sequence-similarity cluster**, not a taxonomic classification — and this distinction matters.

Mechanically:

1. Extract both flanking LTRs from every intact ERV as FASTA.
2. Run all-versus-all BLAST on these LTR sequences.
3. Use the BLAST similarity graph (edges above a configurable identity threshold, typically 80%) to build connected-component clusters.
4. Each cluster is a family. For each family, build a consensus sequence from its member LTRs via multiple alignment.

Families are numbered generically: `family1`, `family2`, ... The numbers are **opaque** — they say "these LTR sequences cluster together by similarity" and nothing more. Two families may correspond to one retroviral genus in taxonomy; one family may split a genus that's diverged heavily; families are **within-genome, de novo** and carry no biological interpretation by themselves.

Output of this stage: `*.LTRlib.fa` — the consensus library, where each FASTA entry is one family's consensus sequence. Header format varies across LTR_retriever releases, but typically includes the family ID and a list of **source ERV IDs** that contributed to the consensus.

Optional: if you supply `-u custom_te_lib.fa`, LTR_retriever can cross-annotate families against a known TE library (Dfam, Repbase, etc.), assigning Copia/Gypsy/… superfamily labels. RetroSeek **doesn't use this** (`noanno: true` in config) because RetroSeek already provides meaningful biological labelling through its probe-based system — see [Coupling B](#coupling-b--multi-label-probe-propagation) below.

### Stage 3 — Solo-LTR discovery via BLAST-back

This is the payoff stage. LTR_retriever takes the consensus library from Stage 2 and **BLASTs it against the whole genome** (using `blastn`). For every genome hit:

- If the hit **overlaps an intact ERV** from Stage 1, it's already accounted for as an LTR flank of that ERV — skip.
- If the hit is **full-length** (close to the consensus's length) and **standalone** (doesn't overlap any intact ERV), it's a **solo LTR**.
- If the hit is **partial** (significantly shorter than the consensus), it's a **truncated LTR fragment**.

Solo LTRs and truncated LTRs are written together to `*.nmtf.pass.list` ("non-matching-full pass list"). Each entry records the genomic coordinates, family assignment, length, identity to the consensus, and strand. RetroSeek's integrator reads this file.

The key conceptual point: **solo LTRs come from the genome matching the consensus, not from filtering LTRharvest false positives**. The family consensus is built from intact ERVs. A solo LTR is a real genomic position that used to have an intact provirus but has since lost its internal region via homologous recombination — leaving a sequence that still matches the family's LTR consensus.

---

## How RetroSeek couples with LTR_retriever

RetroSeek brings two things LTR_retriever doesn't: (1) a probe-based classification that separates retroviral from non-retroviral LTR-retrotransposons, and (2) meaningful biological labels (GAG, POL, ENV, PRO, etc.) anchored to the user's probe set. We use those at two points in the LTR_retriever flow.

### Coupling A — retroviral-only pre-filter

**Problem.** LTR_retriever's Stage 1 keeps all structurally-sound LTR retrotransposons, including Copia/Gypsy/BEL/Pao. Their LTRs cluster into families in Stage 2, those families' consensuses are BLASTed back in Stage 3, and the solo LTRs reported include non-retroviral solos. Most of them, for typical eukaryotic genomes.

**Solution.** Before handing LTRharvest output to LTR_retriever, we **intersect the SCN file with RetroSeek's `valid_ranges.gff3`** — the domain-validated retroviral ERV track. Only LTRharvest candidates that overlap a retroviral-confirmed valid range survive into LTR_retriever's input. As a result:

- LTR_retriever's Stage 1 sees a candidate pool restricted to retroviral-confirmed positions.
- Stage 2's consensus library is built from retroviral LTR sequences only.
- Stage 3's BLAST-back hits match retroviral consensuses.
- Solo LTRs reported are guaranteed retroviral by construction.

This is implemented in `workflow/scripts/ltr_retriever_prefilter.py`. From a single read pass over the input `.scn` + `.des` (LTRharvest's sequence-index descriptor, mapping seq-nr back to chromosome names) + `valid_ranges.gff3`, the script emits **both** SCN files:

- `data/ltr_scn/{genome}_retroviral.scn` — rows whose `[s(ret), e(ret)]` overlaps a `valid_ranges` interval on the same chromosome. Closed-interval, 0-indexed semantics. This is the file LTR_retriever currently consumes.
- `data/ltr_scn/{genome}_full.scn` — every well-formed source row, byte-equivalent to the LTRharvest stdout (modulo malformed rows). This is the unfiltered passthrough, kept on disk for inspection and as the backing file for the upcoming `source_scn: full` mode.

Both files always materialise; which one drives LTR_retriever is becoming a runtime config decision (see `config.ltr_retriever.restrict_to_retroviral` today; an enum `source_scn: retroviral | full` lands in a follow-up). Keeping the full SCN on disk makes ad-hoc inspection trivial and supports side-by-side comparisons.

The filename underscore separator (`_retroviral.scn`, `_full.scn`, not dotted forms) avoids a Snakemake wildcard-resolution collision with LTRharvest's own `{genome}.scn`.

### Coupling B — multi-label probe propagation

**Problem.** LTR_retriever's family IDs (`family1`, `family2`, …) are sequence-similarity clusters, opaque without biological context. We want solo LTRs labelled with the probe families RetroSeek knows about (GAG, POL, ENV, etc.), because:

- Those labels drive hotspot analysis, plot categorisation, and probe-pair detection.
- Users think in probe-family terms, not in LTR_retriever's internal cluster IDs.

**Solution.** `workflow/scripts/solo_ltr_integrator.py` propagates probe labels from RetroSeek's `valid_ranges.gff3` onto LTR_retriever's solo LTRs using a **hybrid two-tier approach**.

#### Primary path — consensus-family mapping

For each solo LTR:

1. Look up its family in `LTRlib.fa` (the headers list source-ERV IDs contributing to the consensus).
2. Match those source-ERV IDs against `valid_ranges.gff3` — the RetroSeek ERVs that seeded this family's consensus.
3. Collect probe labels from every matching valid ERV.
4. Assign the union of those labels to the solo LTR.
5. Record `label_source=family` and list the contributing ERVs in the GFF3 attributes.

This is biologically principled: the sequence similarity that clustered the consensus drives the label inheritance. A solo LTR matched to `family1` inherits labels from the specific intact ERVs whose LTRs are in family1, regardless of where those ERVs sit on the chromosome.

#### Fallback path — nearest-ERV

The primary path may fail for several reasons:

- The solo LTR's family is unresolved (e.g., LTR_retriever couldn't confidently classify it).
- `LTRlib.fa` header format doesn't surface source-ERV IDs parseably (version-dependent).
- Source-ERV IDs from LTR_retriever (e.g., `LTR_retrotransposon5`, which references LTRharvest's candidate numbering) don't match RetroSeek's ID nomenclature in `valid_ranges.gff3`.
- A family exists but none of its source ERVs are in `valid_ranges.gff3` (e.g., filtered out at RetroSeek's domain-validation stage).

When the primary path yields zero labels, the integrator falls back:

1. Find the **nearest valid ERV on the same chromosome** within `config.ltr_retriever.nearest_erv_max_distance` bp (default 10 kb).
2. Inherit its probe labels.
3. Record `label_source=nearest_erv` and list the single nearest ERV as the contributor.

Nearest-ERV is a weaker signal — retroviruses don't always integrate in tight clusters — but it handles the edge cases the primary path doesn't. For a genome where the pre-filter has already constrained families to be retroviral, nearest-ERV is usually correct: solo LTRs tend to cluster near intact ERVs of the same lineage because invasion events have spatial structure.

#### Records left in the output

Every solo LTR in the final `{genome}.gff3` has these attributes:

```
ID=soloLTR_<n>
family=<LTR_retriever family ID or "unresolved">
probe_labels=<comma-separated labels, or "none">
contributing_ervs=<comma-separated ERV IDs, or "none">
label_source=family | nearest_erv | none
```

`label_source` lets downstream analyses filter by confidence. For instance, a strict analysis might keep only `label_source=family` solos; a permissive analysis keeps both.

---

## Solo/intact ratio as an age proxy

With solo LTRs and intact ERVs both labelled by probe family, RetroSeek computes a per-family ratio:

```
solo_to_intact_ratio = solo_count / intact_count
```

High ratios mean the lineage has been present long enough for many of its integrations to have undergone LTR–LTR recombination. Low ratios suggest a recent invasion, where intact proviruses still dominate the inventory. It's one of the classical ways to order retroviral lineages by age without requiring a full molecular-clock analysis.

Two counting modes are emitted side-by-side in the ratio CSV:

- **`label_mode=exclusive`** — counts a solo/intact LTR only if its labels form a single-element set. A record labelled just `POL` counts in POL's row. A record labelled `POL,GAG` doesn't count in exclusive rows. This view is the "pure signal" — it excludes multi-labelled records whose family attribution is ambiguous.
- **`label_mode=shared`** — counts a solo/intact LTR in every family it claims. A record labelled `POL,GAG` counts in both POL's and GAG's rows. This view is the "inclusive signal" — no record is dropped, at the cost of inflating counts for multi-labelled elements.

Both modes are written so downstream analysis can pick either. `label_mode` is a column in the CSV; filter on it to choose.

Per-genome CSVs live at `results/tables/solo_intact_ratio/{genome}.csv`. The aggregate across species is `results/tables/solo_intact_ratio/all_species.csv`, built by an inline `pandas.concat` in the Snakemake aggregate rule.

---

## Configuration & usage

### Prerequisites

LTR_retriever must be in the conda env. After pulling the latest `data/config/environment.yml`:

```bash
mamba env update -f data/config/environment.yml
# or for a first-time install
mamba env create -f data/config/environment.yml
```

### Running solo-LTR detection

With a config pointing at genomes + probes + `execution.entrez_email` set:

```bash
./RetroSeek --solo-ltr-detection --cores all
```

This triggers the full dependency chain: genome download → BLAST DB → suffix arrays → LTRharvest (with SCN output) → probe extraction → tBLASTn → ranges analysis → pre-filter → LTR_retriever → solo-LTR integrator → ratio aggregate. Any stage already completed is skipped by Snakemake's freshness tracking.

Expected outputs per genome:

| Path | Content |
|---|---|
| `data/ltr_scn/{genome}.scn` | LTRharvest screen-format SCN (all candidates). |
| `data/ltr_scn/{genome}_retroviral.scn` | SCN filtered to retroviral-confirmed candidates (Coupling A). |
| `data/ltr_scn/{genome}_full.scn` | SCN passthrough — byte-equal to LTRharvest stdout, kept for inspection. |
| `results/tracks/ltr_retriever/{genome}/{genome}.pass.list.gff3` | LTR_retriever-filtered intact ERVs. |
| `results/tracks/ltr_retriever/{genome}/{genome}.nmtf.pass.list` | Solo + truncated LTRs. |
| `results/tracks/ltr_retriever/{genome}/{genome}.LTRlib.fa` | Consensus library. |
| `results/tracks/solo_ltr/{genome}.gff3` | **Annotated solo LTRs with probe_labels.** |
| `results/tables/solo_intact_ratio/{genome}.csv` | Per-family solo/intact counts + ratio. |
| `results/tables/solo_intact_ratio/all_species.csv` | Aggregated across all genomes. |

### Tunable parameters

All under `config.ltr_retriever` in `data/config/config.yaml`. See [`docs/configuration.md`](configuration.md#ltr_retriever) for the full field reference. Brief summary:

| Key | Default | Effect |
|---|---|---|
| `substitution_rate` | `1.3e-8` | Mammalian bp-substitutions/site/year — used by LTR_retriever for age estimation. Use `7e-9` for plants. |
| `min_ltr_similarity` | `91` | Percent identity floor for LTR pairs. LTR_retriever's `-miniden`. |
| `threads_per_genome` | `4` | CPU threads passed to each LTR_retriever invocation. |
| `noanno` | `true` | Skip LTR_retriever's internal TE library annotation (RetroSeek has its own). |
| `restrict_to_retroviral` | `true` | **Coupling A toggle.** Pre-filter SCN by `valid_ranges.gff3`. |
| `nearest_erv_max_distance` | `10000` | bp window for Coupling B's nearest-ERV fallback. |

Probe-label aggregation across contributors follows the same strategy vocabulary as `parameters.aggregation` (see ADR-002). Configure via `parameters.solo_ltr_aggregation`.

---

## Caveats and known limitations

### LTR_retriever output-format variance

Between LTR_retriever releases, the header format of `LTRlib.fa` and the column layout of `nmtf.pass.list` have changed slightly. Our integrator parser is defensive — it uses regex-heuristic extraction for family IDs and source-ERV lists, and falls back to nearest-ERV if the primary parse yields no labels. If a new LTR_retriever version produces a format our parser can't extract, solo LTRs will still be produced (via fallback) but with `label_source=nearest_erv` rather than `family` for every record. If you see this in `solo_ltr/{genome}.gff3`, either update the parser or pin LTR_retriever to a known-good version.

### Synthetic toy genomes typically produce zero solos

The toy genomes at `/mnt/v/databases/toy-genomes/` (generated by `tests/fixtures/build_toy_genomes.py`) plant tightly-paired LTRs with identical sequences — both flanks of every planted ERV are perfect matches of each other. LTRharvest finds all of them as intact, LTR_retriever keeps them all as intact, and Stage 3's BLAST-back finds no solos because there are no unpaired LTR-like positions. This is expected. Real bat genomes will produce solo LTRs; toy genomes exist for rule-chain smoke-testing, not biological realism.

### Pre-filter excludes non-retroviral LTR-retrotransposons

By design — Coupling A's entire purpose. If you want Copia/Gypsy/other LTR-retrotransposon solos, set `config.ltr_retriever.restrict_to_retroviral: false`. LTR_retriever will then see the unfiltered LTRharvest output, and the solo LTRs include every LTR-retrotransposon class its consensus library captures. You'll likely need to adjust the probe-label propagation strategy since those solos won't have RetroSeek-known probe labels.

### Solo/intact ratios are approximations

The solo/intact ratio is a **crude** proxy for lineage age. Confounders:

- **Unequal discovery sensitivity.** Solo LTRs degrade faster than intact ERVs (single-copy sequences accumulate mutations without gene-conversion repair), so very ancient lineages may have solo LTRs diverged past the detection threshold.
- **Integration preference.** Some lineages integrate into heterochromatin, where repair is less efficient; others prefer euchromatin.
- **Non-homologous deletion.** Some intact ERVs are lost via deletion rather than LTR–LTR recombination, reducing both counts.

For rigorous lineage dating, use a proper molecular-clock analysis on the LTR sequences themselves (LTR_retriever's age estimates, available via the `substitution_rate` config key, are one such mechanism).

---

## Further reading

- **LTR_retriever paper:** Ou & Jiang 2018, *Plant Physiology* 176:1410–1422. [doi:10.1104/pp.17.01310](https://doi.org/10.1104/pp.17.01310). Explains the three-stage pipeline and the family-clustering algorithm in detail.
- **LTRharvest paper:** Ellinghaus et al. 2008, *BMC Bioinformatics* 9:18. [doi:10.1186/1471-2105-9-18](https://doi.org/10.1186/1471-2105-9-18). Covers the suffix-array-based paired-LTR detection.
- **ADR-003:** [`docs/adr/ADR-003-ltr-retriever-pre-filter.md`](adr/ADR-003-ltr-retriever-pre-filter.md) — decision record for Coupling A.
- **Configuration reference:** [`docs/configuration.md`](configuration.md) — the `ltr_retriever` section.
- **Architecture overview:** [`docs/architecture.md`](architecture.md) — RetroSeek's full rule graph including the LTR_retriever branch.
