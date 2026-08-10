# ADR-007: Per-locus ERV taxonomic classification (genus calls)

- **Status**: Accepted (genus-rooting superseded in part by [ADR-008](ADR-008-rank-agnostic-classification.md))
- **Date**: 2026-06-18
- **Deciders**: Jorge González García

## Context

RetroSeek detects ERV loci with high confidence (probe tBLASTn + LTR-domain validation) but **assigned virus/genus poorly**. The legacy assignment transferred the *best-bitscore probe label* onto each locus. Retroviral `pol`/RT is conserved across genera, so a single ERV is homologous to many reference viruses at once; picking the single strongest probe hit over that cross-genus cloud is an uncalibrated argmax (and the alternative - concatenating the whole genus set - is uninterpretable). Neither *classifies*. The probe set is also a detection tool, not a taxonomy: it cannot resolve genera it has no probe for, and it conflates "what found this locus" with "what this locus is".

We needed a taxonomic assignment that is **calibrated, reproducible, probe-/gene-agnostic, mosaic-aware**, and works on mammalian genomes - without weakening detection.

## Decision

Add a **per-locus classification stage** that reclassifies each valid LTR-element locus from its *own* marker sequence against an **independent, genus-comprehensive reference**, per gene, with a method cascade:

1. **Loci** are the LTR_retrotransposon elements: valid-tier features are grouped by their `Parent=` (emitted natively by `range_analysis/validation.R`), then gene-partitioned.
2. **Per gene**, each marker region is blastx'd against the reference, then classified by:
   - **phylogenetic placement** (MAFFT -> EPA-ng -> gappa) onto a per-gene reference tree, for the configured `placement_genes` - but only when placement *resolves a genus* (else fall back, avoiding over-backoff on divergent HERVs);
   - **weighted-LCA** (MEGAN/detectEVE top-percent paradigm over the reference taxonomy) otherwise;
   - **presence** for genes that are diagnostic of a single genus (auto-detected from the reference, e.g. REX/TAX -> *Deltaretrovirus*).
3. **Per locus**, the per-gene calls are combined by marker reliability (the `parameters.main_probes` order, placement preferred), yielding `genus_call` + `rank` + `confidence` + `method`, plus `is_mosaic` / `mosaic_composition` when member genes disagree.

The resulting per-locus table **is the genus-founded ERV assembly**: each LTR-element locus is one provirus carrying its consensus genus, mosaic composition, and structural metrics (`completeness`, `canonical_order`, `n_main_genes`). The legacy probe `virus`/`label` are retained as **detection provenance** (`probe_label_set`), not the taxonomic truth.

**Reference = build-from-NCBI rules** (`taxonomy_reference` + `taxonomy_reference_trees`), not a committed blob: RefSeq proteins per genus via Entrez (curated fixed genus list, balanced per (genus, gene)), an NCBI-derived `taxonomy.tsv`, and the one curated piece - `erv_class.tsv` (Class I/II/III; Jern/Blomberg - not an NCBI rank) - tracked under `data/config/`. The build writes a **provenance manifest** (genera, counts, content hash, Biopython version).

**Defaults**: `placement_genes: [POL]` (POL's tree is reliable); GAG is shipped but opt-in (its reference alignment is low-identity ~19.7%). Search is `blastx` (reuses BLAST+ already in the env - no new search dependency).

### Alternatives considered

- **Genome-wide DIAMOND scan (Mode B)** - abandoned: redundant with the high-quality valid tier, and added a dependency. Everything is Mode A (valid-tier).
- **Commit a frozen reference blob** - rejected: the build-from-NCBI rule + manifest is auditable and self-updating, and the Entrez non-determinism is contained (build-once cache; downstream trees are `--seed`-deterministic from `parameters.seed`).
- **Placement as an override** - rejected: it sent ~89% of human loci to subfamily. Placement now only *refines* when it resolves a genus.

## Consequences

**Positive**

- **Calibrated, literature-consistent genus calls.** Trial on 5 model genomes: bats + human **Class I (gamma) dominant**, mouse **Class II (beta) dominant** - the human<->mouse inversion reproduced; **96.7%** leave-one-out genus accuracy.
- **Reproducible.** Trees seeded from `parameters.seed`; the classification layer regenerates byte-identically; a manifest pins the reference snapshot.
- **Probe-/gene-agnostic.** Gene reliability, mosaic set, and diagnostic genes are all data-derived (`main_probes` + reference), never hard-coded - a REX-only run resolves REX loci.
- **Mosaic-aware**, with IGV-ready GFF3/BED tracks and a plot panel concordant with the tables by construction.

**Negative / costs**

- Adds external tools (`mafft`, `iqtree`, `raxml-ng`, `epa-ng`, `gappa`) to the env. Per-locus marker regions are cut with Biostrings (Bioconductor, `extract_region_fasta.R`), not bedtools - range/sequence work stays in Bioconductor.
- A one-time network reference build (`make reference` / `RetroSeek --build-reference`) is required before `--classify`.
- GAG placement is low-confidence (kept opt-in); ENV stays an LCA-only marker.

## References

- `docs/configuration.md` (`classification` section), `docs/taxonomy_classification/` (design + trial report), `docs/architecture.md`.
- Maksakova 2006 PLoS Genet; Belshaw 2005 MBE; Hayward 2013 PNAS - composition anchors.
