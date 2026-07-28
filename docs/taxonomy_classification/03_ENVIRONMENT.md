# Environment & Resource Inventory

Ground truth of what is actually runnable here. Captured 2026-06-15. Update when it changes.

## Tools

| Tool | Location | Notes |
|---|---|---|
| `tblastn`/`blastn`/`makeblastdb` | `/usr/bin` (system) | present |
| `hmmsearch`/`hmmalign`/`hmmbuild` | `retroseek` env | present (via LTR_retriever→HMMER dep) |
| `datasets` (NCBI) | `retroseek` env | present |
| `diamond` | `detectEVE` env | present — used for the weighted-LCA engine |
| `mafft`, `epa-ng`, `gappa`, `iqtree`/`FastTree` | **absent everywhere** | needed only for Phase-2 placement; would require new deps |
| `snakemake`, `R`/`Rscript`, `python` | present | base + envs |

Conda envs present: `retroseek` (the pipeline env), `detectEVE`, `phylo-bat103`, others.

## Network

- NCBI eutils **reachable** (`curl -sI https://eutils.ncbi.nlm.nih.gov` → OK).
- `conda search` on bioconda is very slow (timed out at 20s) — installs feasible but slow; prefer
  `mamba`. Pin fetched data to avoid depending on this at runtime.

## Data

- **5 model genomes** in `/mnt/v/databases/testing-genomes/`: `Homo_sapiens`, `Mus_musculus`,
  `Desmodus_rotundus`, `Antrozous_pallidus`, `Molossus_molossus`. FASTAs + `.fai` + BLAST DBs +
  suffixerator indices present (detection partially staged).
- **Local `results/tracks/valid/`**: `valid_ranges` GFF3 present for the 3 model bats (and many other
  bats from prior runs); **no human/mouse** yet.
- `/mnt/v/databases/` also holds `bat103`, `globus-bat103*`, `local`, `toy-genomes`, `bat-experimental`
  — not used (boundary / not needed).
- No Pfam HMM cache locally (pipeline downloads on demand); no obvious bundled viral protein DB found.

## Implications

- Build the classifier on **existing bat `valid_ranges`** immediately (no re-running the ~40h/genome
  `ltrdigest`).
- Human/mouse need detection OR a known-reference-locus validation shortcut (see 02_DECISIONS.md D3).
- Engine = DIAMOND + weighted-LCA (no new deps). Placement deferred to Phase 2.

## Canonical dev-data location (use this; never /tmp)

All classifier dev data lives in **`data/taxonomy_dev/`** (gitignored, stable across sessions):
`reference/` (pinned reference + DIAMOND db), `scans/<Genome>.hits.tsv` (raw DIAMOND hits),
`loci/<Genome>.{modeA,modeB}.csv` (results), `valid_tracks/` (Mode-A inputs). Source genomes
stay at `/mnt/v/databases/testing-genomes/*.fa`. Regeneration commands: `data/taxonomy_dev/README.md`.
