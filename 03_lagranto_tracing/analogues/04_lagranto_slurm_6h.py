#!/usr/bin/env python3
"""LAGRANTO 15-day backward trajectories for one ranked MJ analogue - 6-hourly.

Faithful port of lagranto_top20_slurm.py. The trajectory setup is unchanged
(run_lagranto_analogue.sh: same start-file box/profile, 360 h caltra, trace).
This driver already fans out the 6-hourly start hours [0,6,12,18]; paths come
from config_analogues.sh (scratch under the isolated RERUN tree) and the runner
is the co-located 04_run_lagranto_analogue.sh.

PREREQUISITE (one-time, manual): set  set interval = "6"  in
  $LAGRANTO/goodies/datelist.sh   (was "3")
so LAGRANTO steps on the 6-hourly P/S files.
"""

import argparse
import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd
from joblib import Parallel, delayed


# Paths come from config_analogues.sh (env); defaults = the RERUN tree.
_DEF_FINAL = "${LAGRANTO_RERUN}/analogues_windfix"
_DEF_SCRATCH = "${SCRATCH}/lagranto/RERUN_2026_MAY_6h_windfix"
PFILE_BASE = Path(os.environ.get("ANALOGUE_FINAL_BASE", _DEF_FINAL))
# P/S files are READ from here (default: the existing old tree, see config_analogues.sh);
# trajectories are WRITTEN under PFILE_BASE
PS_BASE = Path(os.environ.get("ANALOGUE_PS_BASE", _DEF_FINAL))
SCRATCH_BASE = Path(os.environ.get("ANALOGUE_SCRATCH_BASE", _DEF_SCRATCH)) / "lagranto_slurm"

SCRIPT_DIR = Path(__file__).resolve().parent
RUNNER = SCRIPT_DIR / "04_run_lagranto_analogue.sh"
# TCLIM traced +/-10 hPa (gradient along the trajectory, stage 05); shared with the event
TRACEVARS_FILE = os.environ.get("TRACEVARS_FILE", str(SCRIPT_DIR.parent / "tracevars_tc10hpa.shellscript"))
HOURS = [0, 6, 12, 18]                             # 6-hourly
DAYS_BEFORE = int(os.environ.get("ANALOGUE_DAYS_BEFORE", "4"))   # days before the analogue date
DAYS_AFTER = int(os.environ.get("ANALOGUE_DAYS_AFTER", "3"))     # days after the analogue date
ANALOGUE_WINDOW_DAYS = DAYS_BEFORE + DAYS_AFTER + 1              # total span (for logging/validation)
LEAD_DAYS_BEFORE_BASE = int(os.environ.get("ANALOGUE_LEAD_DAYS", "15"))   # back-trajectory length (days of data before the first arrival)

# Analogue dates come from the shared analogue_dates.txt (edit that ONE file).
from analogue_dates import ANALOGUE_DATES


def make_paths(rank, analogue_date):
    tag = f"rank{rank:02d}_MJ_{analogue_date.strftime('%Y%m%d')}"
    source = PS_BASE / tag
    if not source.is_dir():
        # the reused tree may have been built from another analogue list (other
        # rank numbers): find this analogue's P/S folder by its date instead
        hits = sorted(PS_BASE.glob(f"rank*_MJ_{analogue_date.strftime('%Y%m%d')}"))
        if len(hits) != 1:
            raise FileNotFoundError(f"P/S folder for {analogue_date.date()} not found (or ambiguous) "
                                    f"in {PS_BASE} (ANALOGUE_PS_BASE): {[h.name for h in hits]}")
        print(f"[{tag}] reading P/S from {hits[0].name} (other rank numbering)", flush=True)
        source = hits[0]
    results = PFILE_BASE / tag / "results_6h"
    tmp_base = SCRATCH_BASE / tag / "tmp_6h"
    results.mkdir(parents=True, exist_ok=True)
    tmp_base.mkdir(parents=True, exist_ok=True)
    return {"tag": tag, "source": source, "results": results, "tmp_base": tmp_base}


def start_datetimes(analogue_date):
    analogue_start = analogue_date - pd.Timedelta(days=DAYS_BEFORE)
    analogue_end = analogue_date + pd.Timedelta(days=DAYS_AFTER)
    dates = pd.date_range(analogue_start, analogue_end, freq="D")
    stamps = [f"{d.strftime('%Y%m%d')}_{hour:02d}" for d in dates for hour in HOURS]
    # Drop first and last start time (the two slots the S-file stage cannot
    # centre a DTDT on at the +/-6 h window edges).
    return stamps[1:-1]


def trajectory_file(datetime_str, paths):
    return paths["results"] / f"trajectory_{datetime_str}"


def trajectory_done(datetime_str, paths):
    return trajectory_file(datetime_str, paths).is_file()


def run_one(datetime_str, paths):
    if trajectory_done(datetime_str, paths):
        print(f"[SKIP] {datetime_str} trajectory already exists")
        return
    tmpdir = paths["tmp_base"] / f"tmp_{datetime_str}"
    result = subprocess.run(
        ["csh", str(RUNNER), datetime_str, str(paths["source"]),
         str(paths["results"]), str(paths["tmp_base"])],
        capture_output=True, text=True, env={**os.environ, "TRACEVARS_FILE": TRACEVARS_FILE},
    )
    if result.returncode != 0 or not trajectory_done(datetime_str, paths):
        tail = "\n".join((result.stdout + result.stderr).splitlines()[-25:])
        print(f"[ERROR] {datetime_str} failed (rc={result.returncode}); "
              f"scratch kept at {tmpdir}\n{tail}")
        return
    shutil.rmtree(tmpdir, ignore_errors=True)
    print(f"[OK] {datetime_str}")


def cleanup_scratch(paths):
    root = (SCRATCH_BASE / paths["tag"]).resolve()
    if SCRATCH_BASE.resolve() not in root.parents or not root.name.startswith("rank"):
        raise RuntimeError(f"Refusing to delete suspicious scratch path: {root}")
    if root.exists():
        shutil.rmtree(root)
        print(f"[CLEAN] Removed scratch tree: {root}")


def run_rank(rank, jobs):
    if not 1 <= rank <= len(ANALOGUE_DATES):
        raise ValueError(f"Rank must be 1-{len(ANALOGUE_DATES)}")
    if not RUNNER.is_file():
        raise FileNotFoundError(f"LAGRANTO runner not found: {RUNNER}")

    analogue_date = pd.Timestamp(ANALOGUE_DATES[rank - 1])
    paths = make_paths(rank, analogue_date)
    stamps = start_datetimes(analogue_date)

    print("=" * 80)
    print(f"Rank: {rank:02d}")
    print(f"Analogue: {analogue_date.date()}")
    print(f"Start times: {len(stamps)} (analogue window -{DAYS_BEFORE}/+{DAYS_AFTER} d, "
          f"hours {HOURS}, first/last trimmed)")
    print(f"P/S source: {paths['source']}")
    print(f"Results: {paths['results']}")
    print(f"Scratch tmp: {paths['tmp_base']}")
    print("=" * 80)

    if not paths["source"].is_dir():
        raise FileNotFoundError(f"P/S source folder missing: {paths['source']}")

    pending = [s for s in stamps if not trajectory_done(s, paths)]
    if not pending:
        print("[SKIP] All trajectories already computed")
        cleanup_scratch(paths)
        return

    print(f"[INFO] Pending trajectories: {len(pending)}/{len(stamps)}")
    Parallel(n_jobs=jobs, backend="loky", verbose=10)(
        delayed(run_one)(s, paths) for s in pending
    )

    missing = [s for s in stamps if not trajectory_done(s, paths)]
    if missing:
        raise RuntimeError(
            f"{len(missing)}/{len(stamps)} trajectories missing after run: {missing[:5]}..."
        )
    cleanup_scratch(paths)
    print(f"[DONE] Rank {rank:02d}: {len(stamps)} trajectories in {paths['results']}")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rank", type=int, required=True)
    parser.add_argument("--jobs", type=int, default=32)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.jobs < 1:
        raise ValueError("jobs must be positive")
    os.makedirs(SCRATCH_BASE, exist_ok=True)
    run_rank(args.rank, args.jobs)


if __name__ == "__main__":
    main()
