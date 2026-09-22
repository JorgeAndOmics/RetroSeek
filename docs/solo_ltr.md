# Solo LTRs

## What a solo LTR is, and why they matter

A provirus inserts with two identical LTRs flanking its coding region. Those two
LTRs are a homologous pair sitting a few kilobases apart on the same chromosome, so
they can recombine with each other. When they do, everything between them is
excised and a **single LTR** is left at the integration site.

The provirus is gone. What remains marks the spot.

In mammalian genomes solo LTRs outnumber intact proviruses by one to two orders of
magnitude, and each one records an ancestral integration. A survey that counts only
intact proviruses therefore misses most of the history it is trying to reconstruct.
It also misses it *non-randomly*: recombination has had longer to act on older
insertions, so the integrations that are missing are disproportionately the ancient
ones.

RetroSeek's probe search cannot find them. It looks for retroviral protein, and a
solo LTR has none: the coding sequence is exactly what was excised.

## Why detection is a subtraction, not a classification

There is no sequence feature that says "solo". A solo LTR is simply an LTR, and
what makes it a solo is its *context*: no partner, no internal region.

So the method does not ask "is this a solo?". It finds every copy of a
known-retroviral LTR and then removes the copies that are demonstrably something
else. Every LTR in the genome began as one of a pair flanking a provirus, which
leaves exactly three possible fates, and that is what makes the subtraction
complete:

| fate | second LTR | internal region | is it a solo? |
|---|---|---|---|
| flank of an intact element | present, still pairable | present | no, LTRharvest already catalogued it |
| monoLTR at an orphan | lost or too diverged to pair | **still there** | no, the provirus is damaged rather than excised |
| **solo LTR** | recombined away with | **excised** | **yes** |

The middle row is the one worth dwelling on, because it is where this method beats
the alternative. A lone LTR sitting beside surviving *pol* is not a solo: that
provirus was mutilated, not cleanly excised by recombination. Distinguishing the
two requires a map of retroviral coding sequence that the LTR search missed, which
is precisely what RetroSeek's **orphan tier** is. A tool that only subtracts its own
intact elements has no way to represent that class at all.

## The method

### Bait

The LTR arms of elements that host a catalogued ERV locus, from
`tracks/flanking_ltr/{genome}.gff3`, restricted to elements whose ID appears in the
`parent` column of `{genome}.loci.csv`.

This is what makes the method classifier-free. The probe search already found
retroviral protein inside those elements, so their LTRs are the LTRs of a
retrovirus **by construction**. Every blastn hit is therefore a copy of a
known-retroviral LTR, and nothing downstream has to judge whether a sequence looks
retroviral.

It is the same move the pipeline makes one level up, where the probe CSV holds
retroviral proteins used as tBLASTn bait. Here the bait is nucleotide, and one
level down.

Arms shorter than `solo_ltr.min_bait_length` are dropped. See "Why the length floor
is not optional" below.

### Search

`blastn` against the per-genome nucleotide database the pipeline already builds
(`blast_db_generator`). No database is created for this stage and no new dependency
is introduced.

### Acceptance

A hit is accepted when it satisfies all of:

| criterion | config key | source |
|---|---|---|
| alignment covers 0.8 to 1.2 of the bait arm | `min_coverage`, `max_coverage` | Ou and Jiang 2018, as published |
| alignment is at least 80 bp | `min_alignment_length` | Ou and Jiang 2018, as published |
| bait and hit are both near full length | `min_bait_length`, `min_hit_length` | ours, and required |
| identity to the bait is at least 95% | `min_identity` | ours, and required |

### Subtraction

Surviving hits are merged into candidate loci per sequence, then classified, first
match winning:

1. overlaps an `LTR_retrotransposon` from the LTRdigest track, so it is a flanking
   arm of a catalogued element;
2. lies within `solo_ltr.orphan_pad` of an orphan locus, so it is a monoLTR beside
   surviving coding sequence;
3. otherwise, it is a **solo LTR**.

Element overlap is checked first on purpose. An intact element that happens to have
an orphan nearby is still a catalogued element, and calling it a monoLTR would
remove it from the solo/intact ratio's denominator and inflate the ratio.

### Taxonomy

Each solo inherits `taxon_call`, `rank`, `segment` and `erv_class` from the
classified locus of the element whose arm caught it. That is **homology**, not
proximity: the solo's DNA matched that element's LTR, wherever the two sit on the
chromosome. A nearest-locus fallback exists for the rare case where the seeding
element resolves to no classified locus, and it is recorded as a distinct
`label_source` so it can be filtered.

## Why the length floor is not optional

With the published criteria alone, the method returns **199,816 solo candidates for
Desmodus rotundus**, a ratio of 492 solos per intact element. That is an order of
magnitude above anything in the literature, so it is wrong.

The cause is a mismatch of assumptions. Ou and Jiang's coverage rule is a fraction
of a family **consensus**, which is full length by construction. Our bait is
individual arms running 102 to 999 bp, and 80% of a 102 bp arm is 82 bp. An 82 bp
LTR-derived match, in a genome several percent LTR by mass, is not evidence of
anything.

Requiring bait and hit to be at least 300 bp at 95% identity brings the same genome
to **11,197 solos at 27.6:1**, which is inside the published range and within about
18% of what LTR_retriever reports for the same genome by a completely different
route.

## What the method gives up

**Ancient solos.** The identity requirement is an age filter: it keeps solos whose
family still has a close modern relative. A curated family consensus sits closer to
the ancestral sequence than any surviving copy and would reach further back. This
is the one capability no native implementation matches.

**Target-site duplication does not rescue it.** Recombination happens between the
LTRs and leaves the TSD intact, so a flanking direct repeat should be positive
evidence. Measured against a matched random null it is enriched 2.23x, confirming
the candidates are not noise, but the absolute rate is 1.7% against 0.8%: a 4 to
6 bp exact repeat does not survive tens of millions of years.

**Family names are our own.** No `MLT1A` or `LTR12C`, because there is no curated
library involved.

## Outputs

All paths are under the configured results and data roots. `{genome}` is one of the model 5. Every
file is a declared Snakemake output, so a missing or stale one is rebuilt.

**Figures** (`results/plots/classification/solo_ltr/`)

| file | what it is |
|---|---|
| `{genome}.solo_ltr.pdf` | twelve pages: funnel, identity by fate, length vs identity, distance to orphan, candidates per sequence, solos per seeding element, divergence as time, the drawn tree, tree enrichment against its null, the family census, the solo-only tree, and the largest family subtrees. Pages 8 to 12 only when `tree.enable` |
| `all_species.solo_ltr.pdf` | two pages: solo:intact per genome against the published range, and the three fates per genome |

**Tables for people** (`results/tables/solo_ltr/`)

| file | what it is |
|---|---|
| `{genome}.solo_ltr.csv` | one row per solo, in `catalog.csv`'s column vocabulary: coordinates, `taxon_call`, `rank`, `segment`, `erv_class`, `library_id` (the seeding element) and `label_source` |
| `{genome}.candidates.csv` | **every** candidate, all three fates, with best identity, hit count, seeding bait and element, and distance to the nearest orphan. The evidence behind every call |
| `{genome}.funnel.csv` | counts at every stage, where each rejected hit failed, and the solo:intact ratio |
| `{genome}.ratio.csv` | solo:intact per taxonomic segment |
| `{genome}.tree_summary.csv` | the tree's positive control, clustering against its null, tip census, seed |
| `{genome}.tree_adjacency.csv` | which fate sits beside which on the tree, as enrichment |
| `{genome}.tree_families.csv` | one row per LTR family cut from the tree: its kind (no intact member, with one, or no solos), what it holds, its diameter, and whether it is drawn |
| `{genome}.tree_tips.csv`, `.tree_segments.csv` | the full tree's drawing coordinates, each tip with its class and family |
| `{genome}.solo_tree_tips.csv`, `.solo_tree_segments.csv` | the solo-only tree's drawing coordinates |
| `{genome}.family_tree_tips.csv`, `.family_tree_segments.csv` | the drawn family subtrees' coordinates |
| `solo_report.csv` | one row per genome: solos, monoLTRs, intact flanks, intact loci, ratio |

**Genome-browser track** (`results/tracks/solo_ltr/`)

| file | what it is |
|---|---|
| `{genome}.gff3` | one `solo_LTR` feature per solo, taxonomy in the attributes, for IGV |

**Trees** (`results/trees/solo_ltr/`)

| file | what it is |
|---|---|
| `{genome}.treefile` | the tree, Newick, for FigTree or iTOL |
| `{genome}.solos.treefile` | the tree pruned to its solos, Newick |
| `{genome}.tips.bed`, `.tips.fna`, `.tips.afa` | which tips were sampled, their sequences, and the alignment |
| `{genome}.iqtree`, `.log` | IQ-TREE's report and log: model parameters, likelihood, the seed |
| `{genome}.mldist`, `.bionj`, `.ckp.gz` | IQ-TREE byproducts (distance matrix, starting tree, checkpoint) |
| `{genome}.uniqueseq.phy` | written by IQ-TREE only when some tips have identical sequences |

**Working files** (`data/`, not meant for reading)

| file | what it is |
|---|---|
| `solo_bait/{genome}.bait.bed`, `.bait.fna` | the bait: LTR arms of ERV-bearing elements |
| `solo_blast/{genome}.hits.tsv.gz` | raw blastn hits, gzipped. The only large output: 3.4 GB for the model 5, 1.4 GB of it Homo. Kept so thresholds can be re-swept without re-running blastn |
| `tables/solo_ltr/{genome}.solo_list.tsv` | the six-column hand-off from detector to annotator |
| `tables/solo_ltr/{genome}.*.parquet` | parquet twins of the CSV tables, for the pipeline |
| `tables/solo_ltr/{genome}.manifest.yaml` | provenance: input md5s, every threshold used, the funnel counts |

**Logs** (`logs/solo_bait_builder/`, `solo_blaster/`, `solo_finder/`, `solo_annotator/`, `solo_tree/`, `solo_tree_views/`): one `{genome}.log` each.

The candidate table keeps all three fates deliberately. Evidence is recorded, not
gated (the ADR-015 principle): distance to the nearest orphan is a column, because
in a repeat-dense region a genuine solo can sit near an unrelated orphan by chance,
and a reader should be able to see that rather than have it decided for them.

## The evidence tree

### What it is for, and what it is not for

The tree answers one question the detector cannot: are the three fates of a lone
LTR real biological classes, or artefacts of where the thresholds were drawn? If
solos were simply mis-called flanking arms, they would scatter among flanking arms
on an LTR phylogeny. If the fates are real, they cluster by LTR family, some
families being far richer in solos than others.

It cannot show LTR families that survive only as solos, and nothing here claims
that. Every solo was caught by one bait arm at 95% identity or more, its **seed**,
and a seed is an intact ERV-bearing element: no detected solo is without an intact
relative in its genome.

**It detects nothing and decides nothing.** A tree cannot say a sequence is a solo,
because being a solo is about genomic context (no partner, no internal region), not
sequence. No output of this stage is computed from the tree: fates come from the
subtraction, and taxonomy comes from the element whose LTR arm caught the solo. The
tree is evidence about the method, read by a person.

### How it is built

One tree per genome, built by rule `solo_tree_setup` in five steps:

1. **Choose the tips** (`solo_tips.py`). Up to `tree.n_solo_tips` (200) solos and
   `tree.n_mono_tips` (200) monoLTRs-at-orphans are sampled from the candidate
   table. Then ERV-bearing elements: the seed of every sampled solo, plus a random
   fill up to `tree.n_element_tips` (300), with **both** LTR arms of each element
   kept. Sampling is seeded (`tree.seed`). Each tip's class is written into its
   name as a prefix (`FLANK__`, `SOLO__`, `MONO__`), and a solo's name also carries
   its seed.
2. **Extract the sequences** with `taxonomy/extract_region_fasta.R`, the same
   extractor the bait uses. They come out as they lie on the + strand.
3. **Align** with `mafft --auto --adjustdirectionaccurately`, which turns every
   sequence to the orientation of the others first (about half of them: an LTR can
   sit on either strand). MAFFT marks a turned sequence with an `_R_` prefix, which
   is stripped so tip names stay parseable.
4. **Infer** with `iqtree -m GTR+G -fast -seed` (`tree.model`, `tree.fast`,
   `tree.seed`): maximum likelihood, two search iterations, no bootstrap. **Single
   threaded on purpose**: with the seed fixed, two 8-thread runs still wrote
   different trees, while 1-thread runs are byte-identical. The genomes build in
   parallel instead, so the whole tree stage stays reproducible end to end.
5. **Measure and lay out.** `tree_stats.py` computes the statistics below;
   `solo_tree_layout.py` turns the tree into drawing coordinates using the same
   `taxonomy/tree_layout.py` code the host and taxon trees use (ADR-011's
   coordinate bridge: Python lays the tree out, R draws it, and no R tree package
   is needed).

**Why elements are sampled rather than every arm taken.** The clustering statistic
is a comparison against class abundance, and it saturates when one class dominates.
An earlier version took every bait arm, which made the Mus musculus tree 96.5%
flanking arms: the permutation null rose to 0.94 and the enrichment collapsed to
1.03x, a number that measured tip composition rather than biology.

**Why the orientation step matters.** Until 2026-09-22 the alignment ran without
it. About half the solos then entered the alignment reverse-complemented, and the
tree put them 0.2 to 3.6 substitutions per site from the arm they match at 95% or
more; orientation predicted the misplacement exactly (Desmodus: all 99 reversed
solos far from their seed, all 101 others beside it). The statistics, families and
conclusions drawn from those trees were withdrawn. The seed control below exists so
that a failure of this kind shows up in the summary.

### What it measures

**The seed control.** A solo sits within about 0.05 substitutions per site of its
seed's arm, so on a correct tree nearly every solo is within 0.1 of it. This is the
sharpest check the tree gets.

**The arm control.** An element's two arms were identical the day it inserted, so
they tend to come out as sister tips. The rate is blurred wherever a family burst
left many near-identical copies (an arm's sister is then as likely another
element's arm), and it cannot see strand errors, since both arms of an element
share a strand.

**Same-class sisters against a permutation null.** For every tip: does its sister
group contain at least one tip of its own class? The observed fraction is compared
with the same fraction after the class labels are shuffled `tree.permutations`
times over the fixed topology. Without that null the number means nothing, because
any structured tree shows some clustering.

Results on the model 5, rebuilt 2026-09-23:

| genome | seed control | arm control | same-class observed | null | enrichment |
|---|---|---|---|---|---|
| Antrozous pallidus | 99% | 61% | 79% | 58% | 1.36x |
| Desmodus rotundus | 100% | 64% | 82% | 56% | 1.48x |
| Homo sapiens | 97% | 62% | 83% | 56% | 1.50x |
| Molossus molossus | 99% | 62% | 79% | 56% | 1.39x |
| Mus musculus | 92% | 60% | 77% | 57% | 1.34x |

The fates cluster well above the null in every genome: solos, monoLTRs and flanking
arms are not spread evenly across LTR families. The class-by-class adjacency table
(`tree_adjacency.csv`) counts every leaf of a tip's sister group, so a bushy tree
lets a few large sister groups dominate it (Mus); read it with that in mind.

### Families, and the two views derived from them

The tree is cut into LTR families (`tree_families.py`): a family is a maximal clade
whose largest tip-to-tip distance is at most `tree.family_max_distance` (0.2
substitutions per site). 0.2 is the transposable-element convention, the 80-80-80
rule's 80% identity, and it was where family counts stopped depending on the cut (measured before the orientation fix, to be re-checked).

Each family with solos is either **with an intact member** (at least one sampled
flanking arm) or **with no intact member**. Because every solo's seed is on the
tree, a family with no intact member can only arise where the tree separates a
solo from its seed: after the rebuild that is 0 solos in Antrozous, 1 in Desmodus
and Molossus, 15 in Homo and 16 in Mus, matching the solos that miss the seed
control. The label is kept as a diagnostic, not as a finding.

| genome | families with solos | with no intact member |
|---|---|---|
| Antrozous pallidus | 13 | 0 |
| Desmodus rotundus | 16 | 1 (1 solo) |
| Homo sapiens | 29 | 13 (15 solos) |
| Molossus molossus | 15 | 1 (1 solo) |
| Mus musculus | 43 | 13 (16 solos) |

Two views are derived from the same tree, with no new inference:

- **The solo-only tree**: the evidence tree pruned to its solos, coloured by family.
  Pruning keeps every relationship the full tree inferred among them. Also written
  as Newick (`{genome}.solos.treefile`) for a tree viewer.
- **Family subtrees**: the `tree.family_panels_per_kind` (3) largest families of
  each kind, each drawn as its own small tree, side by side.

### Honest limits

- **Exploratory, not publication-grade.** `-fast` and no bootstrap, so individual
  branches carry no support values. The statistics above summarise the whole tree
  and are robust to that; a reading of any single clade is not.
- **Unrooted.** IQ-TREE infers an unrooted tree, so the root in the drawing is
  arbitrary. Clade membership is meaningful; left-to-right depth near the root is
  not.
- **A sample.** 200 of up to 67,539 solos. Enough for the statistics, not a census.
- **Family, not genus.** LTRs are short and fast-evolving, so an LTR tree resolves
  families but not genera. Genus comes from protein domains; a solo inherits genus
  through its seeding element.
- **Families depend on the cut.** A family is defined by `family_max_distance`, a
  choice rather than a measurement. Its stability from 0.2 upwards was measured on
  the pre-2026-09-22 trees and has not yet been re-measured on the corrected ones.

## Running it

```bash
./RetroSeek --solo-ltr-detector --configfile /abs/path/config.local.yaml --cores 8
```

Every threshold is a `solo_ltr.*` field in `config.yaml`, documented in
[`docs/configuration.md`](configuration.md#solo_ltr) and validated by
`schema.yaml`. The scripts take them as required arguments and carry no defaults of
their own, so the config file is the single source of truth.

## History

This replaces an LTR_retriever-based route (ADR-003, ADR-005, ADR-013). That route
was retired because its Gypsy versus Retrovirus label cannot be used on mammals:
the classifier is TEsorter's, TEsorter has no *env* model, and mammalian ERVs sit
inside the Gypsy lineage on the gag-pol axis it uses. The consequence was that
95.4% of Desmodus solos derived from families it called Gypsy, while every one of
them inherited a retroviral genus call.

See [ADR-017](adr/ADR-017-native-solo-ltr-detection.md) for the decision, the
calibration against a length-matched random-window null, and the alternatives that
were considered and rejected.

## References

- Ou and Jiang 2018, *Plant Physiology* 176:1410. doi:10.1104/pp.17.01310
- Bailly-Bechet, Haudry and Lerat 2014, *Mobile DNA* 5:13. doi:10.1186/1759-8753-5-13
- Novak et al. 2024, *NAR Genomics and Bioinformatics* 6:113
