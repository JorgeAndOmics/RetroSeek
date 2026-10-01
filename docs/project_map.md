# RetroSeek on one page

What the pipeline does, in what order, what it writes, and which number answers
which question. For the reasons behind each part, follow the ADR links.

## In one paragraph

RetroSeek looks for endogenous retroviruses (ERVs) in genomes. It searches each genome
with protein probes from known retroviruses (tBLASTn), finds LTR retrotransposon
structures independently (LTRharvest, LTRdigest), and joins the two: a probe hit inside
an LTR element is an **LTR-flanked** locus, a probe hit with no LTR element around it is
an **orphan**. Every locus is then classified to a viral genus and compared with the
nearest known virus. Separately, the LTRs of those elements are used as bait to find
**solo LTRs**, the single LTRs a provirus leaves behind when it recombines away.

## The stages, in order

Run with `./RetroSeek <flags> --configfile <config>`; `--downstream` runs every
Analysis and Figures stage. Setup and Discovery are the heavy, run-once steps.

| Phase | Flag | What it does | Main output |
|---|---|---|---|
| Setup | `--download-genomes` | genome FASTAs from NCBI | the genomes |
| Setup | `--download-hmm` | pinned Pfam release | `Pfam-A.hmm` |
| Setup | `--download-dfam` | pinned Dfam curated models (optional) | `Dfam-curated_only-1.hmm` |
| Setup | `--build-reference` | reference virus proteins, taxonomy, gene trees | the classification reference |
| Setup | `--probe-extractor` | probe sequences from NCBI | the probe set |
| Indexing | `--blast-dbs`, `--suffix-arrays` | search indexes per genome | BLAST databases, suffix arrays |
| Discovery | `--ltr-candidates` | LTRharvest: LTR element candidates | element track |
| Discovery | `--ltr-domains` | LTRdigest: domains and signals inside each element | annotated element track |
| Discovery | `--blast` | every probe against every genome | the hit table |
| Analysis | `--ranges-analysis` | joins hits and elements: LTR-flanked loci and orphans | `tables/ranges_analysis/` |
| Analysis | `--domain-scan` | curated Pfam domains in each locus | `tables/domains/` |
| Analysis | `--classify` | genus call and nearest virus per locus; the catalog | `tables/taxonomy_classification/` |
| Analysis | `--segment` | the catalog split by taxon | `tables/taxonomy_classification/segments/` |
| Analysis | `--solo-ltr-detector` | solo LTRs, LTR families, Dfam labels | `tables/solo_ltr/` |
| Analysis | `--hotspot-detection` | genomic windows rich in integrations | `tables/hotspots/` |
| Analysis | `--placement-trees` | placement figures on the reference gene trees | `plots/placement/` |
| Figures | `--generate-global-plots` | homology, integration and structure PDFs, and `highlights.pdf`, a short summary of the catalog | `plots/*.pdf` |

## Which file answers which question

| Question | Where to look | Read this |
|---|---|---|
| What did we find, where, and which genus? | `catalog.csv` | `taxon_call`, `segment`, `confidence_tag`, `source` (LTR-flanked or orphan) |
| Is a locus a complete provirus? | `catalog.csv` | `structure_class`, `completeness`, `canonical_order` |
| Which known virus is it closest to? | `catalog.csv` | `nearest_virus` **with** `nearest_virus_identity`: most loci sit at 40 to 50%, a distant relative, not that virus (figure: "How far from a known virus" in `taxonomy.pdf`) |
| Which loci might be new? | `loss_analysis/novel_candidates/` | loci with no reference hit at all |
| How many solo LTRs? | `{genome}.solo_by_class.csv` | the **ERV LTR** row, not the raw total (see below) |
| Which LTR family is an element or solo in? | `{genome}.ltr_families.csv`, `{genome}.solo_ltr.csv` | `ltr_family`; its Dfam name in `ltr_family_dfam.csv` |
| Where do integrations cluster? | `tables/hotspots/` | `{genome}.hotspots.csv` |
| The main results at a glance? | `plots/highlights.pdf` | loci per host, lineage mix, lineage matrix, focal lineages |

The figures for each answer are in `results/plots/`: one PDF per stage, each opening
with a key page that lists its pages and colours.

## The trap: three solo ratios

Three numbers look alike and are not:

1. **Solo candidates per intact locus** (`solo_report.csv`, the ratio page): every
   solo candidate, whatever its bait was. Most of these are copies of LINE or SINE
   repeats that reached the bait by mistake, so this number is inflated.
2. **Solos per intact element of one family** (`{genome}.ltr_family_ratio.csv`): per
   LTR family, for reading one family at a time.
3. **Solos per bait element of one repeat group** (`{genome}.solo_by_class.csv`): the
   ERV LTR row is the retroviral solo count. Use this one for biology.

## The settings that matter most

| Key | What it changes |
|---|---|
| `species:` | which genomes are studied |
| `parameters.main_probes` | which genes count as main (completeness) |
| `classification.gene_priority` | which gene's call wins when genes disagree |
| `parameters.gene_order` | the expected 5' to 3' gene order (`canonical_order`) |
| `classification.segment_rank` | the rank the catalog is split by (genus by default) |
| `solo_ltr.*` thresholds | what counts as a solo LTR |
| `solo_ltr.families.identity` | how alike two LTRs must be to share a family (0.8) |
| `solo_ltr.families.dfam` | whether families get Dfam names (needs `--download-dfam`) |
| `plots.focal_lineages` | the lineages `highlights.pdf` shows on their own page |

Every key is described in `configuration.md`.

## Where things are

- `results/plots/`: the PDFs, one level deep; per-genome PDFs in a folder per stage.
- `results/tables/`: the CSV tables, one folder per stage.
- `results/tracks/`: GFF3 and BED tracks for a genome browser.
- `results/trees/`: evidence and placement trees.

## Deeper reading

`architecture.md` for the whole design; the ADRs in `adr/` for each decision, most
relevant: ADR-007 and ADR-008 (classification), ADR-017 (solo LTRs), ADR-022 (the
three gene lists), ADR-023 (LTR families and Dfam), ADR-024 (nearest virus).
