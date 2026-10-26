#!/usr/bin/env python3
"""Stage 06 - origin contributions for one ranked MJ analogue.

From the per-arrival lagrangian_terms_*.nc written by stage 05
(<rank>/traj_netcdfs/), build ONE origin-contribution file per analogue. Per
trajectory it stores the exact ARRIVAL point (lon0,lat0,p0), TA + its genesis,
and each PROCESS (advection/adiabatic/diabatic/seasonal) with its OWN genesis
(last sign-change of that process's cumulative contribution; coincides with or
lies downstream of the total-TA genesis).

Output: <rank>/origin_contributions.nc  (dims: time, traj_id)
"""

import argparse
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

# Paths come from config_analogues.sh (env).
#   PFILE_BASE  : INPUT  (traj_netcdfs) — the bu1561 store with the trajectories.
#   ORIGIN_BASE : OUTPUT (origin_contributions.nc) — the NEW uc1275 data dir.
_DEF_FINAL = "${LAGRANTO_RERUN}/analogues_windfix"
_DEF_ORIGIN = "${LAGRANTO_BASE}/RERUN/corrected_lagranto_pipeline_6h/analogues_windfix"
PFILE_BASE = Path(os.environ.get("ANALOGUE_FINAL_BASE", _DEF_FINAL))
ORIGIN_BASE = Path(os.environ.get("ANALOGUE_ORIGIN_DIR", _DEF_ORIGIN))
DELTA_T = 6 * 3600.0
PROCS = ["advection", "adiabatic", "diabatic", "seasonal"]

# Analogue dates come from the shared analogue_dates.txt (edit that ONE file).
from analogue_dates import ANALOGUE_DATES


# --------------------------------------------------------------------------- #
def _at(da, s):
    return da.sel(step=s).drop_vars("step", errors="ignore")


def ta_genesis_step(TA):
    """TA genesis as in Hotz et al. (2024) and the colleague's tanom files: the
    LAST point back from arrival that still has the arrival sign of TA (one step
    after the first opposite sign). If TA never changes sign, the earliest valid
    point, so the whole trajectory is integrated (her age = 360 h)."""
    TA0 = TA.sel(step=0)
    opp = (TA * TA0 < 0) & (TA.step < 0)
    first_opp = TA.step.where(opp).max("step")
    earliest = TA.step.where(TA.notnull()).min("step")
    return xr.where(first_opp.notnull(), first_opp + 1, earliest).fillna(0).astype(int)


def cumulative(ds, gen):
    """Cumulative process contribution from the TA genesis to each step (K).
    A term at step s is the interval (s-1 -> s): only steps after the genesis."""
    step2 = ds.step.broadcast_like(ds["TA"])
    mask = (step2 > gen.broadcast_like(ds["TA"])) & (step2 <= 0)
    return {p: (ds[p].where(mask) * DELTA_T).fillna(0).cumsum("step").where(mask) for p in PROCS}


def proc_genesis_step(cum_p, gen):
    """Per-process genesis: last sign-change of the cumulative contribution
    (closest to arrival); fall back to the TA genesis if it never flips."""
    sign = np.sign(cum_p)
    change = (sign.diff("step") != 0) & sign.shift(step=1).notnull().isel(step=slice(1, None))
    last = change["step"].where(change).max("step")
    return last.fillna(gen).astype(int)


def build_all(ds):
    for v in ["lon", "lat", "p"]:
        ds[v] = ds[v].where(ds[v] != -999.0)
    TA = ds["TA"]
    gen = ta_genesis_step(TA)
    cum = cumulative(ds, gen)
    d = dict(
        lon0=_at(ds["lon"], 0), lat0=_at(ds["lat"], 0), p0=_at(ds["p"], 0), arrival_p=_at(ds["p"], 0),
        TA=_at(TA, 0),
        TA_lon_genesis=_at(ds["lon"], gen), TA_lat_genesis=_at(ds["lat"], gen),
        TA_p_genesis=_at(ds["p"], gen), TA_genesis_step=gen.astype("float32"),
    )
    for p in PROCS:
        gp = proc_genesis_step(cum[p], gen)
        # genesis AT arrival (TA flips in the first 6 h back) -> nothing to integrate, 0 K
        d[p] = _at(cum[p], 0).where(gen < 0, xr.where(TA.sel(step=0).notnull(), 0.0, np.nan))
        d[f"{p}_lon_genesis"] = _at(ds["lon"], gp)
        d[f"{p}_lat_genesis"] = _at(ds["lat"], gp)
        d[f"{p}_p_genesis"] = _at(ds["p"], gp)
        d[f"{p}_genesis_step"] = gp.astype("float32")
    return xr.Dataset(d).reset_coords(drop=True)


def make_paths(rank, analogue_date):
    tag = f"rank{rank:02d}_MJ_{analogue_date.strftime('%Y%m%d')}"
    in_dir = PFILE_BASE / tag / "traj_netcdfs"       # INPUT: trajectories on bu1561
    out_dir = ORIGIN_BASE / tag                      # OUTPUT: origin file in the new data dir
    out_dir.mkdir(parents=True, exist_ok=True)
    return {"tag": tag, "in_dir": in_dir, "out": out_dir / "origin_contributions.nc"}


def run_rank(rank, jobs):
    if not 1 <= rank <= len(ANALOGUE_DATES):
        raise ValueError(f"Rank must be 1-{len(ANALOGUE_DATES)}")
    analogue_date = pd.Timestamp(ANALOGUE_DATES[rank - 1])
    paths = make_paths(rank, analogue_date)

    files = sorted(glob.glob(str(paths["in_dir"] / "lagrangian_terms_*.nc")))
    print("=" * 80)
    print(f"Rank: {rank:02d}  Analogue: {analogue_date.date()}")
    print(f"Input : {paths['in_dir']}  ({len(files)} arrival files)")
    print(f"Output: {paths['out']}")
    print("=" * 80)
    if not files:
        raise FileNotFoundError(f"No lagrangian_terms_*.nc in {paths['in_dir']} (run stage 05 first)")

    parts = []
    for f in files:
        stamp = os.path.basename(f).replace("lagrangian_terms_", "").replace(".nc", "")  # YYYYMMDD_HH
        try:
            t = pd.Timestamp(f"{stamp[:4]}-{stamp[4:6]}-{stamp[6:8]} {stamp[9:11]}:00")
        except Exception:
            print("  skip (bad stamp):", f); continue
        ds = xr.open_dataset(f)
        parts.append(build_all(ds).expand_dims(time=[np.datetime64(t)]))
        ds.close()
    pooled = xr.concat(parts, dim="time", join="outer")
    pooled.to_netcdf(paths["out"])
    print(f"[DONE] Rank {rank:02d}: {len(parts)} arrivals -> {paths['out']} | dims {dict(pooled.sizes)}")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rank", type=int, required=True)
    parser.add_argument("--jobs", type=int, default=8)   # accepted for symmetry; build is serial
    return parser.parse_args()


def main():
    args = parse_args()
    run_rank(args.rank, args.jobs)


if __name__ == "__main__":
    main()
