# RetroSeek - Fixing Taxonomic Assignment of ERV Loci

**What this document is.** A theory write-up of why RetroSeek currently struggles to say *which*
virus/genus an ERV locus belongs to, and a proposed approach to fix it reliably on mammal genomes.
It is written three times, for three audiences:

1. **Part 1 - Plain-language version** (no biology or computing background assumed).
2. **Part 2 - Intermediate version** (some genomics/stats vocabulary, but kept light).
3. **Part 3 - Full technical methods** (the actual proposed algorithm and data structures).

All three describe the *same* problem and the *same* solution, at increasing depth. Read whichever
matches you; they stand alone.

> Status: design write-up only. No code or pipeline behaviour is changed by this document.

---
---

# PART 1 - The plain-language version

## What RetroSeek is trying to do

Buried inside the DNA of every mammal are the fossil remains of ancient virus infections. Millions
of years ago, certain viruses (called **retroviruses** - HIV is a modern example) inserted copies of
themselves into the DNA of egg and sperm cells. Those copies got passed down to every descendant,
generation after generation, and they are still sitting in the genome today, mostly broken and
silent. We call these fossils **ERVs** (endogenous retroviruses).

RetroSeek is a tool that reads a mammal's genome and hunts for these viral fossils. For each one it
finds, we'd like to answer two questions:

- **Where is it?** (its location in the genome)
- **What is it?** (which kind of ancient virus it came from)

## The problem in one sentence

**RetroSeek is very good at the "where", but unreliable at the "what".**

It can point to a spot and say "there is definitely a viral fossil here" with high confidence. But
when you ask "*which* virus?", the answer it gives is often little better than a guess - and that is
the flaw this document is about.

## Why the "what" goes wrong - an analogy

Imagine you're trying to identify an unknown, weathered old coin. You have a reference book with
photos of coins from many countries. You lay your coin next to the book and look for matches.

Here's the catch: coins from *many different countries* all share the same basic features - a round
shape, a face in the middle, a ring of text around the edge. So your worn coin "matches" dozens of
reference photos at once, from lots of different countries. You're staring at a **pile of partial
matches**, not one clean answer.

Now, the way RetroSeek currently decides is roughly: *"Of all these partial matches, which single one
looks sharpest? Call it that country."* But "sharpest match" is a terrible way to choose here,
because the sharpness depends on lighting, photo quality, and how worn each reference photo is - not
on which country the coin is really from. You can easily pick the wrong country because one reference
photo just happened to be clearer.

That is almost exactly RetroSeek's situation:

- The "reference book" is a panel of known viral protein sequences (the **probes**).
- Viral fossils from *different* virus groups genuinely look alike, because the key viral genes
  (especially the one called *pol*) barely changed across virus groups over evolution. They share the
  "round shape and a face" - the deep family resemblance.
- So a single fossil in the genome lights up matches to **many** different reference viruses at once.
  This is the "**pile of overlapping matches**" you've noticed.
- RetroSeek then collapses that pile down to one answer by picking the single **highest-scoring**
  match. And that winner is close to arbitrary, for the same reason "sharpest photo" was arbitrary
  for the coin. The location is solid; the identity is a coin toss dressed up as a decision.

## The idea behind the fix

Stop choosing "the sharpest match from the pile." Instead, **identify the fossil on its own terms.**

Concretely: instead of copying the label off whichever reference happened to score highest, we take
the fossil's **own DNA**, read out the actual gene it contains, and then ask a much smarter question:

> "Given everything we know about how all these virus groups are related to each other - like a giant
> family tree of viruses - where does *this specific fossil* belong on that tree?"

Placing something on a family tree is far more reliable than picking a single best match, because the
tree uses *all* the evidence at once and understands that many groups share common features. It can
say "this clearly sits in the cat-virus branch," or - just as importantly - it can say "this sits
somewhere near the base where three branches meet, so I can only narrow it to that broader group, not
one twig." That honesty is the whole point.

Three things make this trustworthy:

1. **We read the fossil's real sequence**, not a borrowed label. The fossil speaks for itself.
2. **We get a confidence number** with every answer. A shaky placement is reported as shaky - or as
   "unclassified" - instead of being dressed up as a confident wrong answer.
3. **We answer gene-by-gene.** Viral fossils are often *patchwork* - one part of the fossil came from
   one virus group and another part from a different one (viruses swap parts; it's normal). So
   instead of forcing one label onto the whole thing, we report each part's identity and flag the
   patchwork. This is the honest picture of what's actually in the genome, and it directly explains
   why you were seeing "overlapping things collapse into overlapping things."

## What changes for someone using RetroSeek

- Instead of one possibly-wrong virus name per fossil, you get an identity **with a confidence
  score**, and a clear **"not sure"** option when the evidence is weak.
- You can see when a fossil is a **patchwork** of multiple virus groups.
- The tool stays flexible: you can still choose which viral genes to look for. The new method doesn't
  depend on any one particular gene - it works with whatever you give it.

That's the whole story. The rest of this document says the same thing with progressively more
precision.

---
---

# PART 2 - The intermediate version

## The pipeline as it stands

RetroSeek detects ERV (endogenous retrovirus) integrations by **tBLASTn**: it takes a panel of known
retroviral **protein probes** (defined by the user in a CSV - each probe is one protein accession
tagged with a genus/`Label`, a virus `Name`, and a gene `Probe` such as POL/GAG/ENV) and searches
them against each genome. Every hit is a candidate piece of a viral fossil.

Hits are then merged into **loci** using genomic-range tooling (overlapping hits on the same strand
are reduced into a single range). Loci are validated against the genome's own LTR/retrotransposon
annotations (LTRharvest + LTRdigest Pfam domains), and multi-gene loci are chained into "ERV-like"
candidates. A separate statistical step (a negative-binomial GLM) finds **hotspots** where ERVs are
over-represented.

**The location machinery is sound.** The detection-and-merge logic produces robust, reproducible
locus boundaries supported by independent evidence (BLAST hits *plus* LTR structure *plus* domain
annotation).

## Where taxonomy assignment breaks

The taxonomic label of a locus is currently produced as a **by-product of merging**, not by a
dedicated classifier. When overlapping hits collapse into one locus, each hit carries its probe's
`virus`/`Label`, and an **aggregation strategy** decides the surviving label. The default strategy is
`best`: keep the label of the single hit with the **highest bitscore** (with a deterministic
tie-break chain).

This fails for a real biological reason:

- Retroviral *pol* (and its reverse-transcriptase core) is **conserved across genera**. A single
  genomic ERV is genuinely homologous to many reference viruses spanning several genera, so it
  receives a **cloud of cross-genus hits** at one locus. This is expected, not an artefact.
- **Raw tBLASTn bitscores are not comparable across probes.** Bitscore scales with alignment length
  and with each probe's particular divergence from the locus. So "highest bitscore" mostly reflects
  *which probe is longest / least diverged*, not *which genus the locus truly belongs to*.
- Therefore `best` resolves a real biological ambiguity with an **uncalibrated argmax**. The
  alternatives in the current design (`strict` -> "ambiguous"; `concatenate`/`list` -> keep them all)
  don't classify either - they either give up or defer the problem to plots.

Net effect: **location is decided by converging independent evidence; identity is decided by a noisy
tie-break.** That asymmetry is the "fatal theoretical flaw."

## The proposed approach, conceptually

Replace **label transfer** with **sequence classification of the locus itself**, built on three
principles.

### Principle 1 - Classify the locus's own sequence, per marker

The probe hit already tells us, for each gene present at a locus, *which gene* it is and *where* it
sits (coordinates, strand, frame). Use that to extract the **locus's own genomic subsequence** under
each hit and translate it. That translated marker - not the probe - is the query we classify. This
keeps the method **probe-agnostic**: it works for whatever markers the user declares (POL, GAG, ENV,
or anything custom), and hard-codes none of them.

### Principle 2 - Place the marker on a reference phylogeny, with confidence

For each gene marker we maintain a **reference framework**: a set of genus-labelled reference
sequences, a fixed multiple-sequence alignment, and a maximum-likelihood reference tree. The locus
marker is aligned into that reference and **phylogenetically placed** onto the tree. Placement is the
field-standard way retroviruses are classified to genus (the reverse-transcriptase phylogeny
approach). Crucially, placement yields a **confidence** (a likelihood-weight ratio over candidate
branches), which we turn into:

- **a confident genus** when one branch dominates,
- **a broader group (an LCA / higher rank)** when the signal is spread across sibling branches,
- **"unclassified retrovirus"** when nothing places well.

This is the principled replacement for `best`/`strict`: it uses all reference evidence at once and it
*knows when it doesn't know*.

### Principle 3 - Report per gene, and flag mosaics

ERVs are frequently **recombinant mosaics** (e.g. *gag* from one lineage, *pol*/*env* from another).
So classify each marker independently and report a genus + confidence **per gene**, with a
locus-level summary that **flags disagreement** rather than hiding it. This is what actually explains
"overlapping ERV-like sequences collapsing into overlapping ERV-like sequences" - several genera's
signals legitimately coexist, and the honest output is to say so.

## Where the reference framework comes from

To stay probe-agnostic *and* reproducible, the recommended sourcing is **hybrid**:

- A **curated, pinned, versioned seed** of reference panels for the common retroviral genera
  (so the default path is high-quality and reproducible across runs and machines), plus
- **automatic expansion/building** for any custom probe the user introduces (so nothing is locked to
  a fixed gene set).

The two extremes are also viable: a *fully curated* package (highest accuracy, most maintenance) or
*fully auto-derived* from the user's probe accessions + fetched relatives (zero curation, noisier,
needs pinning for reproducibility).

## What you'd gain

- A genus call **with a confidence score** and a real **"unclassified"** outcome.
- **Mosaic detection** instead of a single forced label.
- Determinism/reproducibility preserved (pinned references; deterministic placement settings).
- The existing strengths - robust loci, hotspot stats, LTR coupling - untouched; this slots in as a
  **classification stage** consuming the loci that detection already produces.

---
---

# PART 3 - The full technical version (proposed methods)

## 3.1 Problem statement, precisely

Let a genome `G` yield a set of candidate ERV loci `{L_i}` from tBLASTn of a user-defined probe panel
`P`. Each probe `p  in  P` is a tuple `(accession, gene(p), genus(p), virus(p))`. Detection produces,
per locus, a set of HSPs `H_i = { h }`, each `h` carrying `(probe(h), seqname, start, end, strand,
frame, bitscore, evalue, pct_identity, aln_len)`.

Current taxonomy assignment computes, after range reduction,

```
genus(L_i) = genus( argmax_{h  in  H_i} tiebreak_rank(h) ),  tiebreak_rank ~ (bitscore, qcov, identity, -evalue, pos, label)
```

i.e. **label transfer under an argmax of raw bitscore**. The estimator is statistically inconsistent
for genus identity because:

1. **Marker conservation** -> `H_i` spans multiple `genus(p)` with overlapping score distributions
   (high cross-genus homology of *pol*/RT in particular).
2. **Bitscore non-comparability** -> for a fixed locus, `bitscore(h)` is monotone in `aln_len(probe(h))`
   and in `-divergence(probe(h), L_i)`; it is **not** monotone in `P(genus(p) = true genus | data)`.
   Comparing bitscores across probes of differing length and divergence is uncalibrated.
3. **No null/abstention** -> the estimator always emits a genus; there is no mechanism for LCA-backoff
   or abstention under ambiguity (`strict` degenerates to a non-informative "ambiguous" string;
   `concatenate`/`list` defer resolution to visualization).

The detection estimator (locus *location*), by contrast, integrates independent evidence (BLAST  &  LTR
structure  &  Pfam domain), which is why location is reliable while identity is not.

## 3.2 Design goals & invariants

- **G1 Probe-agnostic.** No gene is privileged. The classifier operates per declared marker
  `gene(p)`; if the user omits POL it must still function on whatever markers exist.
- **G2 Self-evidencing.** Classify the locus's *own* translated sequence, not a borrowed probe label.
- **G3 Calibrated confidence + abstention.** Every call carries a confidence; ambiguity backs off to
  a higher rank (LCA) or to `unclassified`.
- **G4 Mosaic-aware.** Per-gene calls + an explicit mosaic flag at the locus level.
- **G5 Deterministic & reproducible.** Pinned reference artifacts + fixed model/seed settings ->
  byte-stable outputs (consistent with the pipeline's existing determinism contract).
- **G6 Additive.** Slots in as a new stage consuming existing loci; detection, hotspot, and LTR
  coupling are unchanged.

## 3.3 Pipeline placement

```
... detection (tBLASTn) -> range reduction -> validation (valid_ranges) -+
                                                                       v
                                            [ NEW: per-marker classification stage ]
                                                                       |
        per-locus, per-gene: extract -> translate -> align -> place -> resolve -> confidence
                                                                       v
                                  classified loci (genus/gene + LWR + mosaic flag)
                                                                       |
                              +----------------------------------------+
                              v
        downstream: ERV-like assembly, plot dataframes, GFF3 tracks (now carry calibrated taxonomy)
```

The classifier replaces the *taxonomic* role of the `best`/`strict` aggregation strategies. Range
*geometry* (reduction, overlaps, n_hits, hotspot input) is unaffected; only the taxon fields change
provenance.

## 3.4 Stage A - marker extraction

For each locus `L_i` and each distinct gene `g` present in `H_i`:

1. Select the representative HSP(s) for `(L_i, g)`. Default: the HSP maximising query coverage of the
   probe for that gene (coverage, not bitscore, to favour completeness of the marker region), subject
   to a minimum `aln_len` (reuse `parameters.probe_min_length`).
2. Extract the **genomic** interval under the HSP from the genome FASTA, with a small flank, on the
   HSP strand. Translate in the HSP frame (tBLASTn reports frame), yielding the locus's own protein
   marker `q_{i,g}`. Stop-codon-aware: split at internal stops, keep the longest ORF spanning the HSP,
   or retain a gappy translation flagged as `pseudo` (degraded ERVs are expected - degradation itself
   is signal, not a failure).

Output: query markers `{ q_{i,g} }`, one per (locus, gene). No gene is assumed - the set of `g` is
whatever the probe panel produced.

> Reuse note: HSP coordinates/strand/frame already exist in the parsed BLAST objects
> (`seq_utils.blaster_parser` / `obj2dict`); genome FASTA access already exists for the BLAST DB and
> hotspot N-masking. Extraction is new glue over existing data, not new infrastructure.

## 3.5 Stage B - the reference framework (per gene)

For each gene `g`, a versioned **reference package** `R_g = (S_g, A_g, T_g, tau_g, M_g)`:

- `S_g`: genus-labelled reference protein sequences spanning retroviral genera + major ERV clades for
  gene `g`.
- `A_g`: a fixed reference MSA of `S_g` (the alignment the queries are aligned *into*).
- `T_g`: a maximum-likelihood reference tree on `A_g` with fitted model parameters (the placement
  scaffold), pinned.
- `tau_g`: taxonomy map, tip -> {genus, subfamily, family, ...}, enabling LCA backoff at any rank.
- `M_g`: optional per-genus profile HMM(s) for the fast gate (see 3.7) and for marker boundary
  detection.

**Sourcing (recommended: hybrid).**
- **Curated seed:** a pinned, vetted panel for the common genera, shipped/cached as a release
  artifact -> reproducible default, controls GenBank label noise and tree topology.
- **Auto-expansion:** for any gene/probe lacking a curated package, build `R_g` automatically:
  expand probe accessions with homologs (Entrez neighbours or a bundled retroviral protein DB),
  inherit genus labels from `genus(p)` and from fetched relatives' taxonomy, then
  `align -> infer tree -> fit model`, and **pin the result** so subsequent runs reuse it.
- Reproducibility is enforced by hashing `R_g` into the run manifest; placement settings are fixed.

This satisfies G1 (any user gene gets a package, curated or built) and G5 (pinned, hashed).

## 3.6 Stage C - alignment + phylogenetic placement

For each query `q_{i,g}`:

1. **Align into `A_g`** without disturbing it (e.g. `hmmalign` against the gene HMM, or
   `mafft --add`), producing a query row in reference coordinates.
2. **Place** the aligned query onto `T_g` by maximum likelihood (e.g. **EPA-ng** or **pplacer**),
   under `T_g`'s fitted model. Output: a set of candidate edges `{e}` each with a
   **likelihood-weight ratio** `LWR(e)`, plus pendant/distal branch lengths.

Placement, not best-hit, because it (a) conditions on the full reference topology (so shared
cross-genus conservation is modelled, not penalised), (b) yields a probability-like weight per edge,
and (c) localises divergent/degraded markers to the correct clade even with low identity.

## 3.7 Stage C' - optional HMM gate (cost control)

To bound cost on genomes with many loci, optionally gate with profile HMMs first
(`hmmsearch q_{i,g}` vs `M_g`): use the top-scoring genus *clade* to (i) confirm the marker boundary
for extraction and (ii) restrict placement to the relevant subtree of `T_g`. HMM bit scores are
length-normalised and E-value-calibrated, so the gate is sound for *coarse* clade selection even
though it remains best-hit in spirit; the **fine** genus call is still made by placement (G3). This
is the "hybrid" engine; pure placement (skip the gate) and pure-HMM (skip placement) are degenerate
configurations of the same stage.

## 3.8 Stage D - resolution to a taxon + confidence

Given placement weights `{(e, LWR(e))}` and taxonomy `tau_g`:

1. Let `e* = argmax_e LWR(e)`.
2. **Confident genus:** if `LWR(e*) >= theta_high` and `e*` maps to a single genus -> assign `genus(e*)`,
   confidence `= LWR(e*)`.
3. **LCA backoff:** else take the minimal edge set `E*` with cumulative LWR >= `theta_mass`; assign the
   **lowest common ancestor** rank of `tau_g(E*)` (e.g. subfamily or "Gammaretrovirus-like"),
   confidence `= sum_{e in E*} LWR(e)`.
4. **Abstention:** if even the LCA is above family level, or total placed mass `< theta_min` -> assign
   `unclassified_retrovirus`.

Thresholds `(theta_high, theta_mass, theta_min)` are config-exposed with documented defaults; resolution is
deterministic given placement output (G5). This is the formal replacement for the
`aggregation.virus/label` strategies: `best` -> `theta_high` rule; `strict`/"ambiguous" -> LCA backoff with
a *meaningful* rank instead of an opaque string.

## 3.9 Stage E - locus summary & mosaic detection

Per locus `L_i`, collect per-gene calls `{ (g, taxon_{i,g}, conf_{i,g}) }`:

- **Locus taxon:** the per-gene call with the highest confidence (configurable: or a designated
  "anchor" gene if the user marks one, but no gene is anchored by default - G1).
- **Mosaic flag:** TRUE if >=2 genes yield confident calls (`conf >= theta_high`) to **different** genera
  (or different LCA clades). Record the per-gene genera so the mosaic composition is inspectable.
- **Degradation flag:** carried from Stage A (`pseudo` translations / internal stops), since heavily
  degraded markers should temper downstream interpretation.

This makes recombination a first-class, reported property (G4) rather than the cause of silent
mislabelling.

## 3.10 Output schema (additive)

Per-locus fields gain (names illustrative; finalise against `plot_dataframe.R` and GFF3 writers):

```
genus_call            # locus-level resolved taxon (may be a higher rank or 'unclassified')
genus_confidence      # LWR-derived confidence of the locus call
genus_rank            # rank at which the call resolved (genus | subfamily | family | unclassified)
per_gene_taxon        # map gene -> resolved taxon
per_gene_confidence   # map gene -> confidence
is_mosaic             # bool
mosaic_composition    # e.g. "POL:Gammaretrovirus; ENV:Betaretrovirus"
is_degraded           # bool (pseudogenised/stop-interrupted marker)
ref_package_version   # hash/version of R_g used (provenance, reproducibility)
```

Legacy `virus`/`label`/`probe` fields are retained for continuity but are explicitly **provenance of
detection** (which probe found the locus), decoupled from `genus_call` (what the locus *is*). This
decoupling is the conceptual heart of the fix.

## 3.11 Validation strategy (how we'd know it works)

- **Known-truth recovery:** on genomes with curated ERV annotations (e.g. well-characterised
  human/mouse ERV families), measure genus-level precision/recall vs the current `best` baseline.
- **Cross-genus confusability test:** synthetic loci built from a known genus's marker, perturbed to
  N% divergence; verify (a) correct genus while resolvable, (b) graceful LCA backoff as divergence
  rises, (c) abstention beyond the reference's resolving power - never a confident wrong call.
- **Mosaic recovery:** chimeric synthetic loci (gag from genus X, pol from genus Y); verify per-gene
  calls + mosaic flag.
- **Determinism:** identical inputs + pinned `R_g` -> byte-identical classification outputs (extend
  the existing reproducibility checks).
- **Probe-agnosticism:** run with non-POL-only panels (e.g. ENV-only, or a custom marker) and confirm
  the stage produces calibrated calls without code changes.

## 3.12 Risks, costs, open questions

- **Reference quality is the ceiling (G1/G5 tension).** Auto-derived packages inherit probe blind
  spots and GenBank label noise; curated seeds mitigate but need maintenance as retroviral taxonomy
  evolves. The hybrid sourcing is the pragmatic compromise; the curated-vs-auto balance is the main
  decision still open.
- **Compute.** Per-locus align+place is heavier than label transfer. The HMM gate (3.7), per-gene
  caching of identical markers, and parallelism bound this; placement on pre-aligned queries against a
  fixed tree is cheap relative to detection.
- **Degraded/short markers.** Very fragmentary loci may not place; abstention (3.8) is the designed,
  honest outcome - not a bug.
- **Tooling dependency.** Adds `mafft`/`hmmalign` + `epa-ng`/`pplacer` (+ `gappa` for placement
  parsing) to `environment.yml`; HMMER is already present. To be recorded as a dependency decision.
- **Anchor-gene policy.** Default locus summary picks the highest-confidence gene; whether to let users
  designate an anchor gene (without hard-coding one) is an open config-surface question.

## 3.13 One-paragraph summary

The current pipeline decides *where* an ERV is by converging independent evidence (robust) but decides
*what* it is by an uncalibrated argmax of raw bitscores over a biologically-expected cross-genus hit
cloud (unreliable). The fix is to stop transferring probe labels and instead classify each locus's own
translated marker by **phylogenetic placement** into a per-gene, genus-labelled, versioned reference,
**probe-agnostically**, producing **per-gene genus calls with calibrated confidence**, **LCA backoff
and abstention** under ambiguity, and an explicit **mosaic** flag - turning "overlapping things that
collapse to overlapping things" from a failure mode into an accurate, honest description of the locus.
