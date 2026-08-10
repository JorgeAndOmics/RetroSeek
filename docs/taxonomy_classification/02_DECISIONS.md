# Decisions Log - Taxonomy Classification

Append-only. Each entry: what was decided, why, and what would reverse it.

---

## D1 - Engine: weighted-LCA over translated-marker homology (not phylogenetic placement, yet)

**Decision.** Classify each locus by searching its *own* translated sequence against a genus-labelled
retroviral protein DB (`diamond blastx`), then assigning the **lowest common ancestor** of all
strong hits in the retroviral taxonomy, with confidence from the score-mass supporting that node.

**Why.**
- Placement tools (`mafft`, `epa-ng`, `gappa`, `iqtree`) are absent from every conda env, and adding
  them trips the CLAUDE.md "new dependency" stop condition. DIAMOND (in `detectEVE` env), HMMER and
  BLAST (in `retroseek` env) are already present.
- Weighted-LCA is the literature-standard taxonomic-binning method (MEGAN; the detectEVE EVE pipeline
  builds on it). It is **not** best-hit: it integrates all strong hits and backs off to a higher rank
  under cross-genus ambiguity - exactly the Section 3.8 confident/LCA/abstain behaviour, via taxonomy rather
  than a tree.
- Simplest runnable approach; matches the "prefer the simplest solution" project rule.

**What would reverse it.** If LCA cannot resolve genus on cases where placement provably could (e.g.
ERV clades absent from NCBI taxonomy), escalate to Phase-2 phylogenetic placement - but only after
explicit approval to add the four tools to `environment.yml`.

---

## D2 - Reference data sourced via Entrez, pinned as a committed artifact

**Decision.** Build the genus-labelled retroviral protein reference by fetching curated accessions
from NCBI (Entrez; network confirmed reachable), reusing the `probe_extractor.py` fetch pattern.
Commit the resolved FASTA + taxonomy map so runs are reproducible and offline-stable.

**Why.** Network is available now but not guaranteed; pinning the resolved sequences removes the
runtime NCBI dependency and gives byte-stable classification (G5). These are public reference
sequences - committing them is consistent with the "example data only" rule (not the genome catalog,
machine paths, or study probe CSV).

**What would reverse it.** If a maintained offline retroviral protein DB with taxonomy is located
under `/mnt/v/databases`, prefer pinning a dereplicated subset of that instead of bespoke Entrez pulls.

---

## D3 - Validation anchored on human + mouse; bats corroborate

**Decision.** Treat `Homo_sapiens` and `Mus_musculus` as the primary literature-truth validation
(their ERVomes are exhaustively characterised), and the 3 model bats as generalisation/corroboration
(DrERV->Gammaretrovirus in `Desmodus_rotundus` is a concrete anchor present in the probe CSV).

**Why.** Reliable "consistent with literature" requires a ground truth; human/mouse provide it.

**What would reverse it.** Nothing expected; if human/mouse detection proves infeasible here, validate
the engine on curated known reference ERV sequences (Dfam/RepBase exemplars) instead of de-novo loci.

---

## D4 - Data-access boundary honoured

`testing-genomes` (5 model genomes) is sanctioned. The real 102/103-bat study data and `*.local.yaml`
real catalog are not read. Directory listings of public reference areas are used only to locate
resources, never to read study sequence contents.

---

## D5 - Non-lossy reduction: separate geometry from labeling; defer taxonomy

**Decision.** Replace the legacy "merge + best-bitscore aggregate" with a two-part step:
1. **Geometry (lossless):** merge overlapping same-strand ranges into a locus extent, keeping
   a back-pointer (`revmap`) to *every* contributing hit. (`GenomicRanges::reduce(with.revmap=TRUE)`
   / `plyranges::reduce_ranges_directed` - already used in `reductions.R`.)
2. **Evidence retention (no aggregation):** attach the full multiset of member hits
   `(gene, genus, bitscore, ...)` to each locus - grouped **by gene**. Taxonomy is resolved
   downstream by `weighted_lca` **per gene**, with a locus summary + **mosaic flag** when
   genes disagree.

**Why.** The legacy `aggregation: best`/`strict` decides genus *at merge time* by bitscore
rank, destroying the mosaic before it can be examined. Decoupling geometry from labeling and
deferring taxonomy keeps all evidence and makes recombination first-class.

**Status.** Implemented in `taxonomy_classify_genome.py` (per-gene LCA + `is_mosaic`). For the
integrated pipeline this principle should move into `reductions.R` (retain revmap/evidence
instead of best-aggregating the `virus`/`label` fields). Genome-scan (Mode B) loci are mostly
single-gene fragments so report few mosaics; richer multi-gene loci (valid tier) will show more.

## D6 - "Valid hits for all 5" is not achievable; even field = uniform Mode B (candidate tier)

**Decision.** Provide an even field by running the *same* genome-wide DIAMOND classification on
all 5 genomes (identical reference + `--sensitive` + reduction), not true valid-tier hits.

**Why.** A "valid" hit requires `ltrharvest`+`ltrdigest` (LTR-overlap + Pfam-domain) output,
which does not exist for human/mouse/Molossus-new and was excluded from regeneration. The only
even, available-information field is the homology/candidate tier produced uniformly. The two
bats' domain-validated Mode-A results remain as a higher-quality cross-check, not the even field.

## D7 - Taxonomy hierarchy is data-derived from NCBI (not hard-coded)

**Decision.** Build the parent/rank hierarchy from NCBI Taxonomy for whatever genera the
reference contains (`taxonomy_build_hierarchy.py` -> `reference/taxonomy.tsv`); `taxonomy_lca.load_taxonomy()`
loads it and the runners call it. The literal maps in `taxonomy_lca.py` are now only a
test/fallback default. ERV class (I/II/III) stays a small curated map (a Jern/Blomberg biological
grouping, not an NCBI rank - nothing to derive it from).

**Why.** Removes manual maintenance, makes the hierarchy probe-agnostic at the taxonomy level,
and fixes silent staleness (e.g. the spuma genus rename). Verified behaviour-preserving: Desmodus
Mode-A identical (859 loci, Gamma 460 / Beta 302) before vs after.
