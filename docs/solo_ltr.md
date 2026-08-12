# Solo-LTR detection in RetroSeek

This document explains how RetroSeek detects **solo LTRs** - the single-LTR remnants of ancient retroviral integrations that LTRharvest cannot find structurally - by integrating [LTR_retriever](https://github.com/oushujun/LTR_retriever) into the pipeline. It covers the biology, the mechanism, and the coupling with RetroSeek's existing probe-based detection and taxonomic classification.

If you just want to configure or run it, jump to [Configuration and usage](#configuration-and-usage). For the full biological and algorithmic detail, start from the top.

---

## Why solo LTRs matter

An intact endogenous retrovirus (ERV) integrates into the host genome as a provirus: two identical **long terminal repeats (LTRs)** flanking an internal coding region (`gag`, `pol`, `env`, etc.). Over evolutionary time - millions of years for ancient families - the two flanking LTRs, being perfect direct repeats, **recombine homologously**. The recombination event excises the entire internal region, leaving behind a **single LTR** at the original integration site.

This is not a degenerate case: in most mammalian genomes, solo LTRs outnumber intact ERVs by **one to two orders of magnitude**. They are silent genomic fossils of retroviral invasions that left no intact provirus behind. Each solo LTR marks exactly one ancestral integration event. Their count per family is a proxy for the **age and abundance of that lineage** - high solo/intact ratios indicate old, well-expanded invasions.

From a RetroSeek perspective, solo LTRs are valuable because:

- They **extend the ERV inventory** - finding integrations our probe-based tBLASTn search would miss entirely, because there is no retroviral protein left at the locus.
- They **enable dating proxies** - the solo/intact ratio per lineage is a simple, interpretable age signal.
- They **complete hotspot maps** - a hotspot defined only by intact ERVs systematically underrepresents high-activity lineages that have recombined heavily.

## Why LTRharvest alone cannot find them

LTRharvest detects paired-LTR structures: two sufficiently-similar direct repeats flanking an inner region of plausible length. It scans the genome with a suffix-array index and reports candidates with two LTRs. A solo LTR has no pair - only one LTR, surrounded by non-retroviral genomic context. LTRharvest silently ignores it.

No parameter adjustment fixes this. LTRharvest's search primitive **requires** the pair. Finding solo LTRs requires a different approach: start with the paired LTRs LTRharvest *did* find, build a library from them, annotate the whole genome with that library, and then keep the matches that stand alone. That is what LTR_retriever does.

---

## How LTR_retriever finds solo LTRs

LTR_retriever ([Ou and Jiang 2018](https://doi.org/10.1104/pp.17.01310)) post-processes LTRharvest output into a filtered, classified ERV inventory, and can then annotate the whole genome with the library it built.

### Stage 1 - Structural filtering of intact candidates

LTRharvest's reported paired-LTR structures have a high false-positive rate: many are tandem repeats, segmental duplications, or non-LTR elements that resemble the pattern. LTR_retriever filters them with biological criteria: **target site duplication** (a 4-6 bp host duplication flanking a real integration), **LTR similarity in a biological range**, **internal retroviral features** (PBS, PPT, minimum internal length, via HMMER), **tandem-repeat screening** with TRF, and **nesting detection** for an ERV integrated inside another.

Survivors are written to `{genome}.pass.list` (plus `.pass.list.gff3`). These are the **intact ERVs**.

> Stage 1 filters for "real LTR retrotransposons", which includes non-retroviral families like Copia, Gypsy and BEL/Pao. It does not filter for retroviruses - that is what RetroSeek's pre-filter adds, see [Coupling A](#coupling-a---retroviral-only-pre-filter).

### Stage 2 - The LTR library

LTR_retriever clusters the intact elements' sequences into a non-redundant library, `{genome}.LTRlib.fa`. Crucially for RetroSeek, each library sequence is **named after the genomic span of the intact element that seeded it**:

```
>Chr1:106472..118130#LTR/unknown
```

The name is a coordinate, not an opaque `family1` identifier. That is what makes the taxonomy coupling below a plain interval join.

### Stage 3 - Whole-genome annotation, then solo calling

LTR_retriever runs **RepeatMasker** over the genome using that library, producing `{genome}.out` - a table of every genomic region matching a library sequence. Two of its own helpers then turn that into a solo list:

```
bin/find_LTR.pl   -lib {genome}.LTRlib.fa              > {genome}.LTR.info
bin/solo_finder.pl -i {genome}.out -info {genome}.LTR.info > {genome}.solo_list
```

`find_LTR.pl` records where the LTR regions sit inside each library sequence, so `solo_finder.pl` can tell an LTR match from an internal-region match. A match is called solo when it:

- lies in an **LTR region** of the library sequence, not the internal region;
- covers **0.8 to 1.2** of that library LTR's length (a partial match is a truncated fragment, an over-long one is spurious);
- is at least **80 bp** long with a Smith-Waterman score above 300;
- sits at least **300 bp clear** of any internal-region annotation - if the internal region is right there, this LTR is a flank of an intact element, not a solo;
- differs in divergence by more than 4% from other LTRs matching the same library entry nearby, which separates a genuine lone LTR from the two flanks of one element.

Output columns:

```
chrom    start    end    chrom:start..end    library_id    coverage
```

### Where solo LTRs are NOT

`{genome}.nmtf.pass.list` is **not** the solo list, despite the suggestive name. "nmtf" means **non-motif**: these are *intact* LTR-RTs whose termini lack the canonical TGCA motif. LTR_retriever's own result banner labels the file `(Non-TGCA LTR-RTs)` and its summary prints "Total intact non-TGCA LTR-RTs found".

An earlier revision of this workstream read "nmtf" as "non-matching-full" and wired the integrator to that file. Because the stage never ran (see [SCN reconstruction](#scn-reconstruction) below), the error never produced output - but it would have reported intact elements as solo LTRs, with real coordinates and plausible counts. This is recorded in [ADR-013](adr/ADR-013-solo-ltrs-on-the-assembled-catalog.md).

---

## SCN reconstruction

LTR_retriever consumes LTRharvest's screen-format `.scn`. That file is a **stdout redirect**, which makes it the most losable artefact in the pipeline: for *Antrozous pallidus* it was gone, along with every byte of the suffix-array index.

Rebuilding the index was not an option. `ltr_harvester_setup` declares `input: rules.ltr_index_generator.input`, which is an `expand()` over **every** genome, so regenerating one index would invalidate LTRharvest and LTRdigest for all five species - a full-pipeline re-run costing days, to recover one redirect.

It is also unnecessary. The SCN contains nothing the GFF3 lacks:

| SCN column | GFF3 source |
|---|---|
| `s(ret)` / `e(ret)` | `LTR_retrotransposon` start / end |
| `l(ret)` | computed, `end - start + 1` |
| `s(lLTR)` / `e(lLTR)` | first `long_terminal_repeat` child |
| `s(rLTR)` / `e(rLTR)` | second `long_terminal_repeat` child |
| `sim(LTRs)` | `ltr_similarity=` attribute |
| `seq-nr` | `seq_number=` attribute |

`workflow/scripts/solo_ltr/scn_from_ltrharvest_gff3.py` does the reconstruction. Verified against Desmodus rotundus, the one model genome holding both artefacts: all **9,893 data rows byte-identical**, including two-space field separators and two-decimal similarity. Antrozous rebuilds to 26,499 rows in under a second.

Similarity is carried through as **text**, never parsed to float: `90.90` round-tripped through a float renders as `90.9` and silently breaks the match.

The `seq_number=` attribute sits on the same GFF3 row as the chromosome name, so this module also supplies the `seq-nr -> chromosome` mapping the pre-filter used to read from the index's `.des` file - which for Antrozous is zero bytes.

---

## How RetroSeek couples with LTR_retriever

RetroSeek brings two things LTR_retriever does not: a probe-based validation that separates retroviral from non-retroviral LTR-retrotransposons, and a per-locus **taxonomic call** (ADR-007/008). Both are used.

### Coupling A - retroviral-only pre-filter

**Problem.** LTR_retriever's Stage 1 keeps all structurally-sound LTR retrotransposons, including Copia, Gypsy and BEL/Pao. Their sequences enter the library, the library annotates the genome, and the solo LTRs reported are dominated by non-retroviral entries - most of them, in a typical eukaryotic genome.

**Solution.** Before handing LTRharvest output to LTR_retriever, **intersect the SCN with RetroSeek's `valid_ranges.gff3`** - the domain-validated retroviral ERV track. Only candidates overlapping a retroviral-confirmed range survive, so the library is built from retroviral sequence only and every solo it finds is retroviral by construction.

`workflow/scripts/solo_ltr/ltr_retriever_prefilter.py` emits both SCN files from one read pass:

- `data/ltr_scn/{genome}_retroviral.scn` - rows overlapping a valid range. The default LTR_retriever input.
- `data/ltr_scn/{genome}_full.scn` - every well-formed row, byte-equivalent to the source.

Both always materialise; `ltr_retriever.source_scn` picks which one feeds the runner. Keeping the full SCN on disk makes side-by-side comparison trivial.

Measured reduction on the model genomes, which also cross-checks the filter - the survivor count matches the classified LTR-flanked locus count exactly:

| genome | LTRharvest candidates | retroviral | reduction |
|---|---|---|---|
| *Antrozous pallidus* | 26,499 | 905 | 29.3x |
| *Desmodus rotundus* | 9,893 | 406 | 24.4x |
| *Molossus molossus* | 36,119 | 1,069 | 33.8x |

Both the SCN and GFF3 are **1-based closed**, so coordinates compare directly. An earlier revision believed the SCN was 0-based and shifted GFF3 starts by `- 1`, widening every valid interval by one base; the byte-equal reconstruction above disproved that.

### Coupling B - taxonomy inheritance

**Problem.** A solo LTR carries no genes, so RetroSeek's probe-based classification cannot reach it directly. It still needs a lineage, because that is what drives grouping, plots and the catalog.

**Solution.** `workflow/scripts/solo_ltr/solo_ltr_integrator.py` inherits the call from the intact element the solo's sequence came from.

#### Primary path - the library element

Each solo's `library_id` is the genomic span of its seeding element (`Chr1:106472..118130#LTR/unknown`). Parse it, overlap it against `{genome}.loci.csv`, and inherit `taxon_call`, `rank`, `segment` and `erv_class` from the classified locus with the largest overlap. `label_source=library`.

This follows **sequence homology, not proximity**: a solo on chromosome 4 whose sequence matched an element on chromosome 1 inherits that element's lineage, which is exactly the biological claim - they descend from the same invasion.

Overlap rather than exact coordinate equality, because LTR_retriever adjusts element boundaries during filtering, and a RetroSeek locus sits *inside* its LTR element rather than sharing its edges. Only `source=ltr-flanked` loci donate: an orphan has no LTR element, so it cannot have contributed sequence to an LTR library.

#### Fallback - nearest classified locus

If the library name carries no coordinates, or its span overlaps no classified locus, fall back to the nearest classified locus on the same chromosome within `nearest_locus_max_distance` bp. `label_source=nearest_locus`.

Proximity is a weaker signal than homology - retroviruses do not always integrate in tight clusters - so it is recorded distinctly. A strict analysis keeps only `label_source=library`; a permissive one keeps both. A sudden rise in `nearest_locus` is the signal that LTR_retriever has changed its library naming.

Solos matching neither get `label_source=none` and empty taxonomy. They are **never dropped**: the count is a real observation even when the lineage is not resolvable.

---

## Solo LTRs in the catalog

Solo LTRs join `catalog.csv` as a third tier beside `ltr-flanked` and `orphan`, with `source=solo-ltr` and `structure_class=solo_ltr`. `reconcile_catalog` keeps the catalog non-overlapping by tier precedence:

```
ltr-flanked  >  solo-ltr  >  orphan
```

A "solo" overlapping a real provirus is that element's flank, not a solo, so the LTR-confirmed call wins. A solo outranks an orphan because a homology match against a curated library is stronger evidence than a proximity-inferred cluster of gene hits. Precedence is evaluated against survivors only, so a solo beaten by a provirus does not go on to displace an orphan.

Solos enter the catalog but **not** the classification plots: all 22 figures key on gene-content columns a solo has lost by definition, and folding solos in would silently change every one of them.

The tier is opt-in via `classification.include_solo_ltr`, because demanding it makes `--classify` wait on the whole-genome RepeatMasker pass.

---

## Solo/intact ratio as an age proxy

With solos and intact proviruses both carrying a taxonomic call, RetroSeek computes a per-group ratio:

```
solo_to_intact_ratio = solo_count / intact_count
```

where the denominator is the count of **LTR-flanked catalog loci** in that group. Grouping follows the ADR-012 vocabulary - `segment` (default, rolled up to `classification.segment_rank`), `taxon_call`, or `none`.

High ratios mean the lineage has been present long enough for many of its integrations to have undergone LTR-LTR recombination. Low ratios suggest a recent invasion where intact proviruses still dominate. Groups with intact loci but no solos are emitted with `solo_count = 0` rather than dropped - a lineage whose proviruses have not recombined away is a finding, not a missing row.

Per-genome CSVs live at `results/tables/solo_intact_ratio/{genome}.csv`, the aggregate at `all_species.csv`.

---

## Configuration and usage

### Prerequisites

LTR_retriever must be in the conda env (it pulls Perl, HMMER, TRF, cd-hit and RepeatMasker as transitive dependencies):

```bash
mamba env update -f data/config/environment.yml
```

### Running solo-LTR detection

```bash
./RetroSeek --solo-ltr-detection --configfile /abs/path/to/config.local.yaml
```

This triggers the chain: SCN reconstruction -> pre-filter -> LTR_retriever (including whole-genome RepeatMasker) -> solo_finder -> integrator -> ratio aggregate. Snakemake skips anything already complete.

**Runtime.** The whole-genome RepeatMasker pass dominates and is the reason this stage is slow on mammalian assemblies. `ltr_retriever.threads_per_genome` is passed straight through as RepeatMasker's `-pa`, so raise it. Always dry-run first (`-n`) and read the `Job stats` table.

Expected outputs per genome:

| Path | Content |
|---|---|
| `data/ltr_scn/{genome}_from_gff3.scn` | SCN reconstructed from the LTRharvest GFF3. |
| `data/ltr_scn/{genome}_retroviral.scn` | SCN filtered to retroviral-confirmed candidates (Coupling A). |
| `data/ltr_scn/{genome}_full.scn` | Unfiltered passthrough, kept for inspection. |
| `results/tracks/ltr_retriever/{genome}/{genome}.pass.list` | Intact LTR-RTs (the ratio denominator's source). |
| `results/tracks/ltr_retriever/{genome}/{genome}.out` | RepeatMasker whole-genome annotation. |
| `results/tracks/ltr_retriever/{genome}/{genome}.LTRlib.fa` | The LTR library. |
| `results/tracks/ltr_retriever/{genome}/{genome}.solo_list` | Raw solo calls from `solo_finder.pl`. |
| `results/tracks/solo_ltr/{genome}.gff3` | **Solo LTRs annotated with `taxon_call`.** |
| `results/tables/solo_ltr/{genome}.solo_ltr.csv` | Per-solo table in catalog vocabulary. |
| `results/tables/solo_intact_ratio/{genome}.csv` | Per-group solo/intact counts + ratio. |
| `results/tables/solo_intact_ratio/all_species.csv` | Aggregated across genomes. |

### Tunable parameters

All under `config.ltr_retriever`; see [`docs/configuration.md`](configuration.md#ltr_retriever) for the full reference.

| Key | Default | Effect |
|---|---|---|
| `substitution_rate` | `1.3e-8` | Mammalian bp-substitutions/site/year, for LTR_retriever's age estimates. `7e-9` for plants. |
| `min_ltr_similarity` | `91` | Percent identity floor for LTR pairs (`-miniden`). |
| `threads_per_genome` | `4` | Threads for LTR_retriever, and RepeatMasker's `-pa`. The main runtime lever. |
| `source_scn` | `retroviral` | **Coupling A toggle.** `retroviral` or `full`. |
| `group_by` | `segment` | Grouping for the ratio table: `segment`, `taxon_call` or `none`. |
| `nearest_locus_max_distance` | `10000` | bp window for Coupling B's nearest-locus fallback. |

Plus `classification.include_solo_ltr` (default `false`) to fold solos into `catalog.csv`.

There is no `noanno` key. Whole-genome annotation produces the RepeatMasker table solo detection reads, so this stage always runs it.

---

## Caveats and known limitations

### Zero solos means something is wrong

On a real mammalian genome, solo LTRs normally outnumber intact proviruses by 1-2 orders of magnitude. A genome with LTR-flanked loci but no solos indicates a broken search, not a quiet genome, and the integrator logs a warning saying so. Check the RepeatMasker section of the LTR_retriever log first.

Synthetic toy genomes are the honest exception: `tests/fixtures/build_toy_genomes.py` plants tightly-paired identical LTRs, so every planted element is found intact and nothing stands alone. Toy genomes exist for rule-chain smoke-testing, not biological realism.

### Sensitivity is bounded by probe quality

Coupling A means the library is built only from probe-validated retroviral elements. If RetroSeek's probe set misses a retroviral lineage entirely, no intact element from that lineage reaches `valid_ranges.gff3`, nothing from it enters the library, and its solo LTRs are invisible. This amplifies the importance of comprehensive probe design.

### Non-retroviral LTR-retrotransposons are excluded

By design - Coupling A's entire purpose. Set `source_scn: full` for Copia/Gypsy/other solos, and expect most of them to have no RetroSeek taxonomy.

### Solo/intact ratios are approximations

The ratio is a **crude** age proxy. Confounders:

- **Unequal discovery sensitivity.** Solo LTRs degrade faster than intact ERVs (single-copy sequence accumulates mutations without gene-conversion repair), so the oldest lineages may have solos diverged past detection.
- **Integration preference.** Lineages differ in their chromatin targeting, and repair efficiency differs with it.
- **Non-homologous deletion.** Some proviruses are lost by deletion rather than LTR-LTR recombination, reducing both counts.
- **Denominator scope.** Intact counts come from catalogued LTR-flanked loci, which require retroviral *gene* evidence, so the denominator is itself bounded by probe coverage.

For rigorous dating, use a molecular-clock analysis on the LTR sequences themselves.

---

## Further reading

- **LTR_retriever paper:** Ou and Jiang 2018, *Plant Physiology* 176:1410-1422. [doi:10.1104/pp.17.01310](https://doi.org/10.1104/pp.17.01310).
- **LTRharvest paper:** Ellinghaus et al. 2008, *BMC Bioinformatics* 9:18. [doi:10.1186/1471-2105-9-18](https://doi.org/10.1186/1471-2105-9-18).
- **ADR-003:** [`docs/adr/ADR-003-ltr-retriever-pre-filter.md`](adr/ADR-003-ltr-retriever-pre-filter.md) - the pre-filter decision.
- **ADR-013:** [`docs/adr/ADR-013-solo-ltrs-on-the-assembled-catalog.md`](adr/ADR-013-solo-ltrs-on-the-assembled-catalog.md) - SCN reconstruction, the solo-source correction, taxonomy inheritance, and the catalog tier.
- **Configuration reference:** [`docs/configuration.md`](configuration.md).
- **Architecture overview:** [`docs/architecture.md`](architecture.md).
