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

| path | content |
|---|---|
| `tracks/solo_ltr_native/{genome}.gff3` | the solo track, taxonomy in the attributes |
| `tables/solo_ltr_native/{genome}.solo_ltr.csv` | one row per solo, in the catalog's column vocabulary |
| `tables/solo_ltr_native/{genome}.candidates.csv` | **all three fates** with their evidence, not just solos |
| `tables/solo_ltr_native/{genome}.funnel.csv` | counts at every stage, plus the solo/intact ratio |
| `tables/solo_ltr_native/{genome}.ratio.csv` | solo/intact per taxonomic segment |
| `tables/solo_ltr_native/{genome}.tree_summary.csv` | the tree's controls and clustering statistics |
| `trees/solo_ltr/{genome}.treefile` | the LTR evidence tree, Newick |
| `plots/classification/solo_ltr/` | the figure panel |

The candidate table keeps all three fates deliberately. Evidence is recorded, not
gated (the ADR-015 principle): distance to the nearest orphan is a column, because
in a repeat-dense region a genuine solo can sit near an unrelated orphan by chance,
and a reader should be able to see that rather than have it decided for them.

## The evidence tree

The stage builds an LTR nucleotide phylogeny over all three fates: every bait arm,
plus a seeded sample of solos and monoLTRs.

**It detects nothing.** A tree cannot tell you a sequence is a solo, because being a
solo is about genomic context. What it can say is whether the three fates are real
classes, and it carries a positive control that comes for free: an element's two
arms were identical the day it inserted, so they must be sister tips. On Desmodus,
284 of 406 elements (70.0%) recover them as sisters, which says the tree carries
real signal.

The substantive result is that solos sit beside other solos 2.16x more often than
class abundance predicts, and beside flanking arms 0.55x as often. Some LTR families
in this genome exist predominantly or entirely as solos, with no intact
representative for LTRharvest to find, which is evidence the method reaches
genuinely new material rather than rediscovering what was already catalogued.

The tree is **exploratory**: `-fast`, no bootstrap, and a sample of the solos. No
classification decision depends on it, and none should.

## Running it

```bash
./RetroSeek --solo-ltr-native --configfile /abs/path/config.local.yaml --cores 8
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
