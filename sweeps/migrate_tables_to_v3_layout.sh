#!/usr/bin/env bash
# =============================================================================
# migrate_tables_to_v3_layout.sh — one-time table-layout migration
# =============================================================================
# RetroSeek's table I/O was reorganised: every table group now lives in its own
# subdirectory, with the pipeline-internal Parquet copy under
# `data/tables/<name>/` and the user-facing CSV copy under
# `results/tables/<name>/`.
#
# An existing results tree from the *old* layout has table files at the old
# paths. Snakemake would otherwise see the new-path outputs as "missing" and
# try to re-run the upstream rules — including `blast_pkl2parquet`, which
# re-reads the cached tBLASTn pickles and fails if the environment's Biopython
# has dropped `Bio.Blast.Record`.
#
# This script relocates the pickle-blocked table groups (`probe_dict`,
# `full_genome_blast`) to their v3 homes so those upstream rules stay satisfied,
# and tidies `segmented_species`. It is idempotent — safe to re-run.
#
# It does NOT touch `plot_dataframes/`, `stage_dataframe/` or
# `tables/manifest/`: those were superseded by the Phase-2 rename
# (`ranges_analysis_setup` now regenerates that data as the `ranges_analysis`
# table group + `results/manifest/`). The script only *flags* them for removal.
#
# Run with the RetroSeek conda env active (needs PyYAML), from the repo root:
#   bash sweeps/migrate_tables_to_v3_layout.sh data/config/config.yaml
# The config's `root.data_root_folder` / `root.results_root_folder` must be
# absolute paths (the norm for server / local configs).
# =============================================================================
set -euo pipefail

CONFIG="${1:?usage: migrate_tables_to_v3_layout.sh <config.yaml>}"
[[ -f "$CONFIG" ]] || { echo "config not found: $CONFIG" >&2; exit 1; }

# Resolve the data/results roots from the config (same idiom as setup_and_run.sh).
eval "$(CONFIG_PATH="$CONFIG" python -c '
import os, yaml
cfg = yaml.safe_load(open(os.environ["CONFIG_PATH"]))
print("DATA_ROOT="    + cfg["root"]["data_root_folder"])
print("RESULTS_ROOT=" + cfg["root"]["results_root_folder"])
')"

echo "RetroSeek table-layout migration (v3)"
echo "  data root    : $DATA_ROOT"
echo "  results root : $RESULTS_ROOT"

# move <src> <dst> — relocate a file if present; idempotent (no-op once moved).
move() {
  local src="$1" dst="$2"
  if [[ -e "$src" ]]; then
    mkdir -p "$(dirname "$dst")"
    mv "$src" "$dst"
    echo "  moved  $src"
    echo "      -> $dst"
  elif [[ -e "$dst" ]]; then
    echo "  ok     already at $dst"
  else
    echo "  skip   $src (not present)"
  fi
}

echo
echo "== Required: pickle-blocked groups (keeps blast_pkl2parquet/probe_extractor satisfied) =="
move "$RESULTS_ROOT/tables/probe_dict.parquet"        "$DATA_ROOT/tables/probe_dict/probe_dict.parquet"
move "$RESULTS_ROOT/tables/probe_dict.csv"            "$RESULTS_ROOT/tables/probe_dict/probe_dict.csv"
move "$RESULTS_ROOT/tables/full_genome_blast.parquet" "$DATA_ROOT/tables/full_genome_blast/full_genome_blast.parquet"
move "$RESULTS_ROOT/tables/full_genome_blast.csv"     "$RESULTS_ROOT/tables/full_genome_blast/full_genome_blast.csv"

echo
echo "== Tidy: segmented_species parquet copies move under data/tables/ =="
echo "   (species_segmenter re-runs cheaply and regenerates both formats anyway;"
echo "    this just keeps results/tables/ clean in the meantime.)"
if [[ -d "$RESULTS_ROOT/tables/segmented_species" ]]; then
  mkdir -p "$DATA_ROOT/tables/segmented_species"
  for f in "$RESULTS_ROOT/tables/segmented_species/"*.parquet; do
    [[ -e "$f" ]] || continue
    mv "$f" "$DATA_ROOT/tables/segmented_species/$(basename "$f")"
    echo "  moved  $(basename "$f")"
  done
fi

echo
echo "== Superseded by the Phase-2 rename — no longer read or written =="
echo "   These hold only regenerable derived data; ranges_analysis_setup now"
echo "   produces results/tables/ranges_analysis/ + results/manifest/ instead."
echo "   Safe to delete once you are satisfied:"
for d in "$RESULTS_ROOT/tables/plot_dataframes" \
         "$RESULTS_ROOT/tables/stage_dataframe" \
         "$RESULTS_ROOT/tables/manifest"; do
  [[ -d "$d" ]] && echo "     rm -rf $d"
done

echo
echo "Migration complete. probe_dict + full_genome_blast are now at their v3"
echo "paths, so blast_pkl2parquet / probe_extractor stay satisfied and the"
echo "tBLASTn pickles are never re-read. species_segmenter + ranges_analysis"
echo "re-run (no pickles involved) and emit the v3 layout."
