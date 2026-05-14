#!/usr/bin/env python3
"""
RetroSeek parameter sweep driver.

Runs `ranges_analysis_setup` + `plot_generator_setup` once per row in
profiles.csv, with all upstream rules (LTRharvest / LTRdigest / BLAST)
left untouched. Each row's outputs are moved into sweeps/runs/<run_id>/
so the production results dir is reusable for the next iteration.

Usage:
    # Activate the retroseek conda env first.
    python sweeps/sweep.py                  # run every row
    python sweeps/sweep.py 01_baseline 04_id80   # run a subset by run_id
"""

import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

REPO         = Path(__file__).resolve().parent.parent
SWEEPS       = REPO / "sweeps"
PROFILES_CSV = SWEEPS / "profiles.csv"
BASE_CONFIG  = REPO / "data" / "config" / "config.local.yaml"

CORES         = 8
MERGE_OPTION  = "virus"   # locked per spec

# Subdirectories under `results/` to move out after each run.
# `tracks/ltrharvest`, `tracks/ltrdigest`, `tracks/ltr_retriever`,
# `tracks/solo_ltr`, `tracks/hotspots`, `tables/segmented_species`,
# `tables/probe_dict`, `tables/full_genome_blast` are intentionally
# left in place — they are upstream and time-heavy to regenerate.
MOVE_RESULTS_SUBDIRS = [
    "tracks/original",
    "tracks/candidates",
    "tracks/valid",
    "tracks/flanking_ltr",
    "tables/overlap_matrix",
    "tables/ranges_analysis",   # user-facing CSV copies of the per-genome tables
    "manifest",                 # run manifests (moved out of tables/)
]

# Subdirectories under `data/` to move out after each run — the pipeline-internal
# parquet tables the plot generators consume. The next profile would otherwise
# overwrite them.
MOVE_DATA_SUBDIRS = [
    "tables/ranges_analysis",
    "tables/overlap_matrix",
]

# Snakemake's `--forcerun` removes the output dirs of forced rules before
# rebuilding. The R script does not create its own parent dirs, so we
# re-create them on every iteration to keep the rule-shell happy.
ENSURE_RESULTS_DIRS = MOVE_RESULTS_SUBDIRS + ["plots"]
ENSURE_DATA_DIRS = MOVE_DATA_SUBDIRS


def deep_set(cfg: dict, dotted_key: str, raw_value: str) -> None:
    parts = dotted_key.split(".")
    cur = cfg
    for p in parts[:-1]:
        cur = cur.setdefault(p, {})
    # Coerce numeric-looking strings to int / float.
    try:
        v: object = int(raw_value)
    except ValueError:
        try:
            v = float(raw_value)
        except ValueError:
            v = raw_value
    cur[parts[-1]] = v


def render_config(row: dict) -> dict:
    with open(BASE_CONFIG) as f:
        cfg = yaml.safe_load(f)
    for k, raw in row.items():
        if k == "run_id" or raw == "":
            continue
        deep_set(cfg, k, raw)
    cfg["parameters"]["merge_option"] = MERGE_OPTION
    return cfg


def move_outputs(results_root: Path, data_root: Path, run_dir: Path) -> None:
    # results/ subdirs land under run_dir/results/...; data/ subdirs under
    # run_dir/data/... — both roots carry a `tables/ranges_analysis`, so the
    # prefix keeps the CSV and parquet copies from colliding in the run dir.
    for root, subdirs, prefix in (
        (results_root, MOVE_RESULTS_SUBDIRS, "results"),
        (data_root, MOVE_DATA_SUBDIRS, "data"),
    ):
        for sub in subdirs:
            src = root / sub
            dst = run_dir / prefix / sub
            if not src.exists():
                continue
            dst.parent.mkdir(parents=True, exist_ok=True)
            if dst.exists():
                shutil.rmtree(dst)
            shutil.move(str(src), str(dst))

    # Top-level PNGs from plots/ only — leave circle_plots/ and
    # hotspot_pdfs/ subdirs in place since those rules aren't being run.
    plots_src = results_root / "plots"
    plots_dst = run_dir / "results" / "plots"
    plots_dst.mkdir(parents=True, exist_ok=True)
    if plots_src.exists():
        for p in plots_src.glob("*.png"):
            shutil.move(str(p), str(plots_dst / p.name))


def run_one(row: dict) -> bool:
    run_id  = row["run_id"]
    run_dir = SWEEPS / "runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    cfg = render_config(row)
    cfg_path = run_dir / "config.yaml"
    with open(cfg_path, "w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)

    results_root = Path(cfg["root"]["results_root_folder"])
    data_root = Path(cfg["root"]["data_root_folder"])
    for sub in ENSURE_RESULTS_DIRS:
        (results_root / sub).mkdir(parents=True, exist_ok=True)
    for sub in ENSURE_DATA_DIRS:
        (data_root / sub).mkdir(parents=True, exist_ok=True)

    print(f"\n[{run_id}] config → {cfg_path}")
    print(f"[{run_id}] snakemake invocation (RETROSEEK_CONFIG={cfg_path}):")
    # defaults.py reads its config at import time from RETROSEEK_CONFIG;
    # --configfile alone only affects config[...] lookups inside rules,
    # leaving PATH_DICT / SPECIES bound to the repo-default config.
    cmd = [
        "snakemake", "plot_generator",
        "--configfile", str(cfg_path),
        "--forcerun", "ranges_analysis_setup", "plot_generator_setup",
        "-c", str(CORES),
        "--keep-going",
        "--rerun-triggers", "mtime",
    ]
    print("    " + " ".join(cmd))
    env = {**os.environ, "RETROSEEK_CONFIG": str(cfg_path)}
    # Snakefile shell rules use repo-relative paths like `scripts/ranges_analysis.R`,
    # which resolve correctly only when snakemake's cwd is `workflow/` — matches
    # how the canonical `RetroSeek` wrapper at repo root invokes snakemake.
    r = subprocess.run(cmd, cwd=REPO / "workflow", env=env)
    if r.returncode != 0:
        print(f"[{run_id}] FAILED (exit {r.returncode}); leaving production results in place")
        return False

    move_outputs(results_root, data_root, run_dir)
    print(f"[{run_id}] outputs moved → {run_dir}")
    return True


def main() -> int:
    if not PROFILES_CSV.exists():
        print(f"missing {PROFILES_CSV}", file=sys.stderr)
        return 2
    with open(PROFILES_CSV, newline="") as f:
        rows = list(csv.DictReader(f))

    wanted = set(sys.argv[1:]) if len(sys.argv) > 1 else None
    if wanted is not None:
        rows = [r for r in rows if r["run_id"] in wanted]
        if not rows:
            print(f"no profiles matched {sorted(wanted)}", file=sys.stderr)
            return 2

    print(f"sweep: {len(rows)} profile(s) — {[r['run_id'] for r in rows]}")
    n_ok = 0
    for row in rows:
        if run_one(row):
            n_ok += 1
    print(f"\nsweep complete: {n_ok}/{len(rows)} ok")
    return 0 if n_ok == len(rows) else 1


if __name__ == "__main__":
    sys.exit(main())
