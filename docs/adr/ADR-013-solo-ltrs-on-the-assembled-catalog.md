# ADR-013 - Solo LTRs on the assembled ERV catalog

- **Status**: Accepted
- **Date**: 2026-08-11
- **Deciders**: Jorge Gonzalez Garcia
- **Refines**: [ADR-003](ADR-003-ltr-retriever-pre-filter.md) (retroviral pre-filter),
  [ADR-005](ADR-005-ltr-retriever-runner.md) (runner)
- **Builds on**: [ADR-010](ADR-010-orphan-clustering-and-authoritative-catalog.md),
  [ADR-012](ADR-012-hotspot-on-the-assembled-catalog.md)

## Context

Solo LTRs are the single-LTR remnants left when a provirus's two flanking LTRs
recombine homologously and excise the internal region. In mammalian genomes they
outnumber intact ERVs by one to two orders of magnitude, and each marks one
ancestral integration. They were the last stage ADR-010 deferred.

The stage was fully scaffolded - rules, CLI flag, conda dependency, 262 lines of
`docs/solo_ltr.md`, two ADRs - and had **never run**. `results/tracks/solo_ltr/`,
`results/tracks/ltr_retriever/` and `results/tables/solo_intact_ratio/` were all
empty. Investigating why surfaced four independent defects, three of which the
fourth was hiding.

1. **No input, and repairing it would cascade.** The prefilter needed
   `data/ltr_scn/{genome}.scn` plus the suffix array's `.des`. Antrozous has no
   SCN and every index file is zero bytes after a space-reclaim event.
   Regenerating the index would also have invalidated LTRharvest and LTRdigest
   for **every** genome, because `ltr_harvester_setup` declares
   `input: rules.ltr_index_generator.input` - an `expand()` over all species.
   One missing stdout redirect implied a full-pipeline re-run.
2. **The wrong file.** The integrator read `nmtf.pass.list`, glossed in
   `docs/solo_ltr.md` as "non-matching-full pass list" holding solo and truncated
   LTRs. It is nothing of the sort: LTR_retriever's own banner calls it
   `(Non-TGCA LTR-RTs)` and its summary prints "Total intact non-TGCA LTR-RTs
   found". These are **intact** elements whose termini lack the canonical TGCA
   motif. The stage would have reported intact elements as solo LTRs.
3. **The config disabled the evidence.** Real solo detection runs off
   `bin/solo_finder.pl`, which reads the whole-genome RepeatMasker table
   `{genome}.out`. That file only exists when annotation runs, and
   `ltr_retriever.noanno` was `true`.
4. **Stale vocabulary.** Labels were propagated as `probe_labels` from
   `valid_ranges.gff3` - the pre-ADR-007/008 concept `taxon_call` replaced -
   and the primary family-to-ERV path could never fire, because LTR_retriever
   emits `LTR_retrotransposonN` while `valid_ranges.gff3` carries RetroSeek IDs.
   Every solo would have fallen through to the weak nearest-ERV path, as the
   document's own caveats conceded.

Defect 1 masked 2, 3 and 4: with no output, a confidently-wrong table never
appeared.

## Decision

### Rebuild the SCN from the GFF3 instead of the index

`scn_from_ltrharvest_gff3.py` reconstructs `{genome}.scn` from
`{genome}.gff3`, which carries every SCN field: the element span, its two
`long_terminal_repeat` children, `ltr_similarity=` and `seq_number=`. Verified
against Desmodus rotundus, the one model genome holding both artefacts: all
**9,893 data rows byte-identical**, two-space separators and two-decimal
similarity included. Antrozous rebuilds to 26,499 rows in 0.85 s.

The rule takes the **GFF3 only** as input, never
`rules.ltr_index_generator.input`, so no cascade. `seq_number=` also supplies the
`seq-nr -> chromosome` map, retiring the prefilter's `.des` input - which
mattered, since Antrozous's `.des` is zero bytes.

Similarity is carried as **text, not float**: `90.90` through a float renders as
`90.9` and silently breaks byte-equality.

### Fix the prefilter's coordinate frame

The byte-equal reconstruction proves SCN `s(ret)` equals the GFF3 start exactly.
The prefilter assumed the SCN was 0-based and shifted GFF3 starts by `- 1`,
widening every valid interval by one base. Both frames are now treated as
1-based closed. The corrected filter reproduces the catalog counts exactly:
**406/406** for Desmodus and **905/905** for Antrozous ltr-flanked loci.

### Take solos from the annotation, not from `nmtf.pass.list`

`run_ltr_retriever.py` no longer passes `-noanno` and chains LTR_retriever's own
helpers: `find_LTR.pl -lib {genome}.LTRlib.fa` then
`solo_finder.pl -i {genome}.out -info {genome}.LTR.info`, emitting
`chrom, start, end, locus, library_id, coverage`. Reusing the published criteria
(coverage 0.8-1.2 of the library LTR, minimum 80 bp, 300 bp clear of internal
regions) rather than reimplementing them. `{genome}.out` is a declared output, so
a missing annotation fails loudly instead of yielding zero solos.

### Inherit taxonomy by coordinate, not by family ID

LTR_retriever names each library sequence after the genomic span of the intact
element that seeded it (`>{chr}:{start}..{end}#LTR/{fam}`, `annotate_lib.pl`), so
the library ID is a **coordinate**, not an opaque `family1`. The integrator
overlaps that span against `{genome}.loci.csv` and inherits `taxon_call`,
`rank`, `segment` and `erv_class` from the largest-overlap locus
(`label_source=library`). This follows sequence homology - the solo matched
*that* element's library entry - and turns the primary path from never-fires into
the normal case, deleting the family-ID mapping problem rather than debugging it.

Overlap rather than exact equality, because LTR_retriever adjusts element
boundaries and a locus sits *inside* its element. Only `source=ltr-flanked` loci
donate: an orphan has no LTR element and cannot have seeded the library. The
nearest-locus fallback survives for solos whose library ID lacks coordinates, and
is recorded distinctly so it can be filtered out.

Ratio grouping follows ADR-012: `group_by` = `segment` (default) | `taxon_call` |
`none`, replacing `probe_family`. Groups with intact loci but no solos are
emitted with `solo_count = 0` rather than dropped.

### Solos as the catalog's third tier

`reconcile_catalog` generalises from "ltr-flanked versus everything else" to an
explicit precedence: **`ltr-flanked` > `solo-ltr` > `orphan`**. A solo overlapping
a real provirus is that element's flank, not a solo; a solo outranks an orphan
because homology against a curated library beats proximity-inferred gene hits.
Precedence is evaluated against survivors only, so a solo beaten by a provirus
does not go on to displace an orphan.

Solos join `catalog.csv` but deliberately **not** the plotting frame: all 22
figures key on gene-content columns a solo has lost by definition, and folding
solos in would silently change every one of them.

The coupling is opt-in via `classification.include_solo_ltr` (default `false`).
Demanding solos unconditionally would make every `--classify` run wait on a
whole-genome RepeatMasker pass.

## Consequences

- **Positive**: the stage runs at all; solos carry the authoritative taxonomy at
  any rank; the catalog gains its most numerous tier; a 20 GB index rebuild and a
  five-genome cascade are avoided; a 1 bp prefilter bug is gone; `noanno` and
  `parameters.solo_ltr_aggregation` retire with the vocabulary they served.
- **Negative**: whole-genome RepeatMasker is now on the critical path and is by
  far the slowest step. Mitigated by the 24-34x prefilter reduction (26,499 ->
  905 candidates for Antrozous) and by `threads_per_genome` feeding
  RepeatMasker's `-pa`.
- **Coupling**: `classification.include_solo_ltr: true` makes `--classify`
  depend on the solo-LTR stage.
- **Honest reading of a zero result**: on a real mammalian genome, solo LTRs
  normally outnumber intact proviruses. Zero solos means the search is wrong, not
  that the biology is quiet, and the integrator warns explicitly when a genome
  with LTR-flanked loci yields none.

## Alternatives considered

**Rebuild Antrozous's suffix array.** Rejected: unnecessary once the SCN proved
reconstructible, and it would have cascaded LTRharvest and LTRdigest across all
five genomes for a file that is a pure re-rendering.

**Drop LTR_retriever for a native blastn back-search.** RetroSeek already has
the curated LTR arms, the BLAST databases and `extract_region_fasta.R`, so this
was viable and cheaper. Rejected by the maintainer in favour of keeping the
published tool and its validated solo criteria. Revisit if RepeatMasker runtime
proves untenable at 100-genome scale.

**Keep `noanno: true` and derive solos from the library ourselves.** Rejected for
the same reason: it would have made the solo call ours rather than the published
one, for a stage whose entire value is precision.

## Revisit trigger

- RepeatMasker runtime blocks the 100-genome study.
- LTR_retriever changes its library naming away from `{chr}:{start}..{end}`,
  breaking the coordinate join. The nearest-locus fallback keeps producing
  (weaker) labels, and a jump in `label_source=nearest_locus` is the signal.
- Solos are wanted in the plots, not only the catalog.

## References

- `workflow/scripts/solo_ltr/scn_from_ltrharvest_gff3.py` - SCN reconstruction.
- `workflow/scripts/solo_ltr/ltr_retriever_prefilter.py` - Coupling A.
- `workflow/scripts/solo_ltr/run_ltr_retriever.py` - runner + solo_finder chain.
- `workflow/scripts/solo_ltr/solo_ltr_integrator.py` - taxonomy inheritance.
- `workflow/scripts/taxonomy/taxonomy_plot_generator.R` - `reconcile_catalog`.
- `docs/solo_ltr.md` - biology and mechanism.
- Ou & Jiang 2018, *Plant Physiology* 176:1410-1422. doi:10.1104/pp.17.01310.
