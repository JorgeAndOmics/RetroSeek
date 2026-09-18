# ADR-015: Domain evidence from one symmetric scan, curated by accession

- **Status**: Accepted
- **Date**: 2026-09-17
- **Deciders**: Jorge González García
- **Refines**: [ADR-009](ADR-009-anchored-domain-tiering-and-structure-class.md) (domain tiering), [ADR-010](ADR-010-orphan-clustering-and-authoritative-catalog.md) (orphan tier)

## Context

Two defects, found together.

**Orphans had no domain evidence at all.** Domains reached the catalog only through
`gt ltrdigest -hmms`, which by construction searches inside LTRharvest elements.
Orphan loci sit outside every element, so `taxonomy_classify_loci.py` defaulted
them to `non_domain`. Measured on the catalog: **30,123 of 30,415 `non_domain`
rows (99.04%) had never been assessed**, yet were indistinguishable from the 292
that genuinely meant "looked, found nothing". A `groupby('domain_tier')` reported
71.9% of loci as domain-free, which was not true.

**The LTR-flanked tier's domain labels were wrong.** ADR-009's `domain_tier` was
decided by substring-matching the Pfam domain **name** against regexes in
`config.domains` (`ranges/granges_build.R`, unanchored `grepl`,
`ignore.case = TRUE`). Measured over the 233,874 domain instances in the model 5:

- **43.8% of the retroviral-diagnostic signal was invisible.** `rve` (integrase
  core, 7,111 instances), `RVP` (retroviral protease, 5,070), `IN_DBD_C` (5,139),
  `MLVIN_C` (2,313) and `GP41` (1,721) match no pattern, so
  `extract_domains_with_probes` discarded them. An element carrying only `rve`
  and `RVP` was labelled `domain_unlisted`.
- **43,739 instances were falsely probe-labelled.** The POL pattern `ase` matched
  `Transposase_22` (an L1 ORF1p domain) 29,081 times, against `RVT_1`'s 29,101,
  so roughly 35% of everything labelled POL was LINE machinery. The GAG pattern
  `zf` pulled in 6,425 host zinc fingers.
- The `domains` keys (`POL, GAG, ENV, PR15`) and `parameters.main_probes`
  (`POL, GAG, ENV`) had drifted apart, and the tie-break sorted on probe-key
  length, so `PR15` outranked every real probe.

The regex list began life in ADR-009's predecessor as a **filter** that deleted
LTR-anchored hits without a matching domain. ADR-009 stopped it deleting and made
it label instead, but the mechanism survived into a role it was never designed for.

A third problem sits underneath both: LTRdigest searches the whole element
(6,919 bp mean) while an orphan locus is a hit span (681 bp mean). Equal method
over unequal search space is still biased, so the target is method symmetry, not
merely "add domains to orphans".

## Decision

**Scan both tiers over the catalogued locus span with one procedure, and
interpret the result through one curated table.**

1. **One scan, both tiers.** `hmmsearch --cut_ga`, six-frame translation, over the
   locus span, for LTR-flanked and orphan loci alike. The locus span is the only
   region both tiers possess, so it is the only basis on which a cross-tier
   comparison is a measurement rather than an artifact.
2. **`--cut_ga`, Pfam's curated per-family gathering thresholds.** Bit scores do
   not depend on how the search was framed, so GA results are comparable in a way
   E-values are not. At GA only **39 distinct families** fire across the orphan
   marker sets, against 6,175 reported by LTRdigest at its loose E-value setting,
   which is what makes a ~150-model curated subset very nearly lossless and turns
   a 7.8-hour full-Pfam scan into minutes.
3. **Curation on accessions, not name regexes.** `data/config/pfam_domain_classes.tsv`
   assigns each family to one of `retroviral_diagnostic`, `retroelement_shared`,
   `non_ltr`, `dna_transposon`, `other`. Accession is the key because Pfam
   guarantees accession stability and explicitly reserves renaming.
4. **`domain_tier` keeps its ADR-009 name, values and ordering; only its mechanism
   changes.** `domain_selected` now means "carries a family curated as
   retroviral_diagnostic or retroelement_shared".
5. **Three new columns**: `domain_evidence` (strongest class present, else
   `none`), `domain_names` (distinct families found), `domain_source` (`scan` or
   `not_scanned`). `domain_source` is what makes "no domains" distinguishable
   from "never looked".
6. **Keyed on the SET of distinct families, never on hit counts.** LTRdigest
   chains fragments of one model into a single feature and `hmmsearch` does not,
   so any count-based statistic would differ between them for reasons unrelated
   to biology.
7. **Evidence never gates admission.** No locus is filtered, and `confidence`
   (which means taxon-call confidence) is untouched.
8. **`domain_hit_class` is deleted**, with `hit_domain_mode` and
   `.positional_hit_class`. It was produced, exported to the track, parsed back
   by the classifier, and then never read: not in `LOCI_COLUMNS`, no plot, no
   filter. It only existed because the regex mapped each domain name to a probe.

## LTRdigest is kept, and is not rerun

`gt ltrdigest` is invoked byte-identically. Its `protein_match` and `RR_tract`
features keep feeding the structural completeness panel and the circle plots.
What changes is only that the **classification** signal now comes from a scan
that records accessions and bit scores.

`-pdomcutoff` is deliberately left unset. Source inspection shows LTRdigest sums
`log(E-value)` per strand across all `protein_match` children, picks the winner,
and **deletes** the losing strand's features, so changing the cutoff would
silently flip element strands and alter domain content together. That is far
outside the scope of a domain-evidence feature.

Why scan rather than trust LTRdigest's calls: it ships no domain database
(`-hmms` is opt-in and domain search is off without it), the 2009 paper
(Steinbiss et al., NAR 37:7002, doi:10.1093/nar/gkp759) reports no benchmark of
domain accuracy, it does not reconstruct ORFs, it reports only the optimal chain
(`-allchains no`), **27.6% of the rows it emits exceed its own documented
`-pdomevalcutoff` of 1E-6**, bit scores are parsed internally but never written,
and accessions are never parsed at all.

This split is the mainstream architecture rather than a deviation from it:
LTR_retriever consumes LTRharvest output and delegates domain annotation to
TEsorter's own `hmmscan` against REXdb or GyDB.

## Consequences

- **`domain_tier` values change**, so catalogs from before this ADR are not
  comparable on that column. Movement is expected in both directions: gains from
  the 21,735 previously-invisible diagnostic instances, losses from the 43,739
  falsely probe-labelled ones.
- **Orphan loci gain domain evidence for the first time**, and every locus now
  carries honest provenance in `domain_source`.
- **Two annotations now exist over the same LTR-flanked loci**, LTRdigest's and
  the scan's. Their concordance is reported and is the evidence on which a later
  decision to retire `-hmms` (worth about 20 hours per genome) can rest.
- **Stage-level columns are renamed to say what they mean**: `n_probe_domains` ->
  `n_selected_domains`, `domain_probes` -> `domain_classes`, and the probe x
  domain heatmap now plots the domain's curated class rather than a probe guess,
  so a POL hit sitting on an L1 domain is visible as exactly that.
- **A family absent from the curated table is invisible to the scan.** At GA this
  is close to moot (39 families fire, 33 already listed), but it is not
  self-monitoring: a periodic full-Pfam audit is the documented maintenance step.
- **`config.domains` and `parameters.hit_domain_mode` are removed**, which also
  retires the PR15/PRO namespace mismatch.
