#!/usr/bin/env bash
# =============================================================================
# RetroSeek parameter-sweep bootstrap.
#
# Fill in the CONFIG block below, then run:
#   bash sweeps/setup_and_run.sh
#
# What it does, in order:
#   1. Verifies the conda env exists.
#   2. Verifies data/config/config.local.yaml is sane (path keys + species).
#   3. Verifies upstream artifacts (FASTA, ltrdigest, segmented_species,
#      probe_dict) are on disk for every species in config.local.yaml —
#      missing any of these would silently trigger hours of upstream
#      recomputation.
#   4. Verifies snakemake is >= 8.0 (needed for `--rerun-triggers mtime`).
#   5. Patches CORES in sweep.py to your server's value.
#   6. Dry-runs the first profile to confirm the DAG resolves correctly
#      (exactly 7 jobs per profile, no upstream rules in the plan).
#   7. Asks for confirmation, then launches the full sweep.
#   8. Runs the aggregator after the sweep finishes.
#
# Exits non-zero on any verification failure. Designed so you can run it
# repeatedly: re-running after a partial completion picks up only the
# missing profiles (Snakemake's mtime check + this driver's per-run-dir
# layout handle that automatically).
# =============================================================================

set -euo pipefail

# ==================== CONFIG (edit these) ====================
# Name of the conda environment with snakemake + R retroseek deps.
CONDA_ENV="retroseek"

# Number of cores snakemake should grab per profile invocation.
# Leave ~25% headroom for OS / I/O. On a 32-core box, 24 is reasonable.
CORES=8

# Set NONINTERACTIVE=yes to skip the "Proceed?" prompt (useful for tmux/nohup).
NONINTERACTIVE="no"

# Optional: restrict the sweep to a subset of profiles by run_id. Empty = all.
# Example: PROFILE_SUBSET=(01_baseline 05_id85)
PROFILE_SUBSET=()
# =============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BASE_CONFIG="$REPO_ROOT/data/config/config.local.yaml"
PROFILES_CSV="$SCRIPT_DIR/profiles.csv"
SWEEP_PY="$SCRIPT_DIR/sweep.py"
AGGREGATE_R="$SCRIPT_DIR/aggregate.R"

red()    { printf '\033[31m%s\033[0m\n' "$*"; }
green()  { printf '\033[32m%s\033[0m\n' "$*"; }
yellow() { printf '\033[33m%s\033[0m\n' "$*"; }
banner() { printf '\n\033[1;36m[%s/8] %s\033[0m\n' "$1" "$2"; }
die()    { red "ERROR: $*"; exit 1; }

# -----------------------------------------------------------------------------
banner 1 "Verifying conda env '$CONDA_ENV'..."
if ! conda env list 2>/dev/null | awk '{print $1}' | grep -qx "$CONDA_ENV"; then
  die "conda env '$CONDA_ENV' not found. Available envs:
$(conda env list 2>/dev/null | awk 'NF>1 && !/^#/ {print "  - " $1}')"
fi
green "  OK"

# -----------------------------------------------------------------------------
banner 2 "Verifying base config '$BASE_CONFIG'..."
[[ -f "$BASE_CONFIG" ]] || die "config.local.yaml not found at $BASE_CONFIG"
[[ -f "$PROFILES_CSV" ]] || die "profiles.csv not found at $PROFILES_CSV"
[[ -f "$SWEEP_PY" ]] || die "sweep.py not found at $SWEEP_PY"
[[ -f "$AGGREGATE_R" ]] || die "aggregate.R not found at $AGGREGATE_R"

# Pull path roots + species list out of YAML via python (safer than grep).
# `conda run` doesn't reliably forward stdin through a heredoc, so we pass
# the config path via env var and use `python -c` with single-quoted code.
CONFIG_SUMMARY=$(BASE_CONFIG_PATH="$BASE_CONFIG" conda run -n "$CONDA_ENV" python -c '
import os, sys, yaml
cfg = yaml.safe_load(open(os.environ["BASE_CONFIG_PATH"]))
spp = list(cfg.get("species", {}).keys())
if not spp:
    sys.exit("config.local.yaml has no species entries")
print("DB_ROOT="      + cfg["root"]["db_root_folder"])
print("RESULTS_ROOT=" + cfg["root"]["results_root_folder"])
print("DATA_ROOT="    + cfg["root"]["data_root_folder"])
print("SPECIES="     + ",".join(spp))
') || die "failed to parse $BASE_CONFIG"
eval "$CONFIG_SUMMARY"
IFS=',' read -ra SPECIES_ARR <<<"$SPECIES"

echo "  db_root_folder      : $DB_ROOT"
echo "  results_root_folder : $RESULTS_ROOT"
echo "  active species (${#SPECIES_ARR[@]}): ${SPECIES_ARR[*]}"
green "  OK"

# -----------------------------------------------------------------------------
banner 3 "Verifying upstream artifacts for all ${#SPECIES_ARR[@]} species..."
MISSING=()
for g in "${SPECIES_ARR[@]}"; do
  for p in "$DB_ROOT/$g.fa" \
           "$RESULTS_ROOT/tracks/ltrdigest/$g.gff3" \
           "$DATA_ROOT/tables/segmented_species/$g.parquet"; do
    [[ -e "$p" ]] || MISSING+=("$p")
  done
done
[[ -e "$DATA_ROOT/tables/probe_dict/probe_dict.parquet" ]] \
  || MISSING+=("$DATA_ROOT/tables/probe_dict/probe_dict.parquet")

if (( ${#MISSING[@]} > 0 )); then
  red "  MISSING ${#MISSING[@]} required upstream artifact(s):"
  printf '    %s\n' "${MISSING[@]}"
  red "  Run the upstream rules (ltr_digester, species_segmenter, probe_extractor) first."
  exit 1
fi
green "  OK — every species has fasta + ltrdigest + segmented_species; probe_dict present"

# -----------------------------------------------------------------------------
banner 4 "Verifying snakemake version >= 8.0..."
SM_VER=$(conda run -n "$CONDA_ENV" snakemake --version 2>/dev/null | head -1)
SM_MAJOR=${SM_VER%%.*}
if [[ -z "$SM_MAJOR" || "$SM_MAJOR" -lt 8 ]]; then
  die "snakemake version $SM_VER is too old. Need >= 8.0 for --rerun-triggers mtime."
fi
green "  OK — snakemake $SM_VER"

# -----------------------------------------------------------------------------
banner 5 "Patching CORES=$CORES in sweep.py..."
sed -i.bak -E "s/^(CORES[[:space:]]*=[[:space:]]*)[0-9]+/\1$CORES/" "$SWEEP_PY"
grep -E "^CORES[[:space:]]*=" "$SWEEP_PY" | sed 's/^/  /'
green "  OK"

# -----------------------------------------------------------------------------
banner 6 "Dry-run on first profile..."
FIRST_PROFILE=$(awk -F, 'NR==2 {print $1}' "$PROFILES_CSV")
[[ -n "$FIRST_PROFILE" ]] || die "couldn't read first profile from $PROFILES_CSV"
echo "  Using run_id '$FIRST_PROFILE' for dry-run"

TEMP_CFG="$SCRIPT_DIR/runs/$FIRST_PROFILE/config.yaml"
mkdir -p "$(dirname "$TEMP_CFG")"
SWEEP_SCRIPT_DIR="$SCRIPT_DIR" \
SWEEP_PROFILES_CSV="$PROFILES_CSV" \
SWEEP_FIRST_PROFILE="$FIRST_PROFILE" \
SWEEP_TEMP_CFG="$TEMP_CFG" \
conda run -n "$CONDA_ENV" python -c '
import os, sys, csv, yaml
sys.path.insert(0, os.environ["SWEEP_SCRIPT_DIR"])
import sweep
target = os.environ["SWEEP_FIRST_PROFILE"]
with open(os.environ["SWEEP_PROFILES_CSV"], newline="") as f:
    row = next(r for r in csv.DictReader(f) if r["run_id"] == target)
cfg = sweep.render_config(row)
with open(os.environ["SWEEP_TEMP_CFG"], "w") as f:
    yaml.safe_dump(cfg, f, sort_keys=False)
' || die "failed to render dry-run config"

DRY_OUT=$(cd "$REPO_ROOT/workflow" && \
  RETROSEEK_CONFIG="$TEMP_CFG" conda run -n "$CONDA_ENV" snakemake plot_generator \
    --configfile "$TEMP_CFG" \
    --forcerun ranges_analysis_setup plot_generator_setup \
    -c "$CORES" --rerun-triggers mtime --dry-run 2>&1) || true

# Look only inside the Job stats table (not informational prose around it).
# Snakemake 9 prints upstream rule names in "Rules with missing metadata:"
# notes even when those rules aren't actually scheduled — the prior check
# false-positived on that.
JOB_STATS=$(echo "$DRY_OUT" | awk '/^Job stats:/,/^total/' || true)
BAD_RULES=$(echo "$JOB_STATS" | grep -E "^\s*(ltr_harvester_setup|ltr_digester_setup|blast_pkl2parquet|species_segmenter_setup|full_genome_blaster_setup|probe_extractor|pfam_hmm_downloader)\s+[1-9]" || true)
if [[ -n "$BAD_RULES" ]]; then
  red "  Dry-run would re-run upstream rule(s) with non-zero counts:"
  echo "$BAD_RULES" | sed 's/^/    /'
  red "  This means upstream artifacts aren't where snakemake expects them."
  red "  Re-check the paths printed in step 2 — config.local.yaml's results_root_folder must match where ltrdigest/segmented_species actually live."
  exit 1
fi

N_RANGES=$(echo "$DRY_OUT" | grep -oP "ranges_analysis_setup\s+\K[0-9]+" | head -1)
N_PLOT=$(echo "$DRY_OUT"   | grep -oP "plot_generator_setup\s+\K[0-9]+" | head -1)
N_RANGES=${N_RANGES:-0}
N_PLOT=${N_PLOT:-0}
EXPECTED_RANGES=${#SPECIES_ARR[@]}
if [[ "$N_RANGES" != "$EXPECTED_RANGES" || "$N_PLOT" != "1" ]]; then
  red "  Dry-run job counts wrong: ranges_analysis_setup=$N_RANGES (expected $EXPECTED_RANGES), plot_generator_setup=$N_PLOT (expected 1)"
  echo "$DRY_OUT" | tail -20
  exit 1
fi
green "  OK — dry-run produced $N_RANGES ranges_analysis_setup + $N_PLOT plot_generator_setup, no upstream rules"

# -----------------------------------------------------------------------------
banner 7 "Confirmation"
if (( ${#PROFILE_SUBSET[@]} > 0 )); then
  N_PROFILES=${#PROFILE_SUBSET[@]}
  echo "  Subset mode: $N_PROFILES profile(s) — ${PROFILE_SUBSET[*]}"
else
  N_PROFILES=$(($(wc -l < "$PROFILES_CSV") - 1))
  echo "  Full sweep: $N_PROFILES profile(s) from $PROFILES_CSV"
fi
echo "  Genomes per profile : ${#SPECIES_ARR[@]}"
echo "  Cores               : $CORES"
echo "  Rough ETA           : ~30–60 min/profile (genome-size dependent) ≈ $((N_PROFILES * 45)) min total"

if [[ "$NONINTERACTIVE" != "yes" ]]; then
  echo ""
  read -r -p "  Proceed with sweep? [y/N] " ans
  case "${ans,,}" in
    y|yes) ;;
    *) yellow "  Aborted before sweep launch."; exit 0 ;;
  esac
fi

# -----------------------------------------------------------------------------
banner 8 "Running sweep + aggregator"
cd "$REPO_ROOT"
START_TS=$(date +%s)
LOG="$SCRIPT_DIR/sweep.log"

set +e
if (( ${#PROFILE_SUBSET[@]} > 0 )); then
  conda run --no-capture-output -n "$CONDA_ENV" python -u "$SWEEP_PY" "${PROFILE_SUBSET[@]}" 2>&1 | tee "$LOG"
else
  conda run --no-capture-output -n "$CONDA_ENV" python -u "$SWEEP_PY" 2>&1 | tee "$LOG"
fi
SWEEP_RC=$?
set -e

ELAPSED=$(( $(date +%s) - START_TS ))
echo ""
echo "  Sweep finished in $((ELAPSED / 60))m $((ELAPSED % 60))s (rc=$SWEEP_RC)"

if (( SWEEP_RC != 0 )); then
  yellow "  Sweep had failures. Running aggregator on whatever did succeed..."
fi

conda run -n "$CONDA_ENV" Rscript "$AGGREGATE_R" "$SCRIPT_DIR"

green ""
green "Done. Outputs:"
green "  sweeps/runs/<run_id>/        — per-profile output trees (tracks, plots, tables)"
green "  sweeps/aggregate/             — three CSVs spanning every run × genome × virus × label × category"
green "  sweeps/sweep.log              — full run log for debugging"
exit "$SWEEP_RC"
