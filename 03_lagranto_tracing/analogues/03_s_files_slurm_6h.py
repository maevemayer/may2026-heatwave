#!/usr/bin/env python3
"""Build LAGRANTO S files for ranked MJ analogues - CORRECTED 6-hourly.

Faithful port of s_files_top20_slurm.py. The agreed 6-hourly changes:
  * SLOTS are 00/06/12/18,
  * climatology files are read as t_YYYY-MM-DD_6hHH_tbar.nc,
  * the DTDT central time difference uses +/-6 h (dt_seconds = 2*6*3600),
    because the climatology now lives on a 6 h grid,
  * the 04_02 rename phase is DROPPED: 6-hourly snapshots already carry the
    exact 00/06/12/18 valid times, so there is no +1 h offset to correct and
    no dual naming schedule to support.
THETA, PS, hyam/hybm fold-in (old
04_write_coords.sh), staging/validation and the edge-slot accounting are
UNCHANGED. (--skip-rename is still accepted but is now a no-op.)
"""

import argparse
import os
from pathlib import Path

import pandas as pd
import xarray as xr
from joblib import Parallel, delayed


# Paths come from config_analogues.sh (env); defaults = the RERUN tree.
_DEF_FINAL = "${LAGRANTO_RERUN}/analogues_windfix"
PFILE_BASE = Path(os.environ.get("ANALOGUE_FINAL_BASE", _DEF_FINAL))
CLIM_PATH = PFILE_BASE / "t_clim"

SLOTS = ["00", "06", "12", "18"]                 # was 01/04/07/10/13/16/19/22
DT_HOURS = 6                                      # was 3 (climatology cadence)
FINAL_VARIABLES = {"THETA", "TCLIM", "DTDT", "PS"}   # no DTDP: gradient along the trajectory (stage 05)
HYBRID_COEFF_VARS = ["hyam", "hybm"]

DAYS_BEFORE = int(os.environ.get("ANALOGUE_DAYS_BEFORE", "4"))   # days before the analogue date
DAYS_AFTER = int(os.environ.get("ANALOGUE_DAYS_AFTER", "3"))     # days after the analogue date
ANALOGUE_WINDOW_DAYS = DAYS_BEFORE + DAYS_AFTER + 1              # total span (for logging/validation)
LEAD_DAYS_BEFORE_BASE = int(os.environ.get("ANALOGUE_LEAD_DAYS", "15"))   # back-trajectory length (days of data before the first arrival)

# Analogue dates come from the shared analogue_dates.txt (edit that ONE file).
from analogue_dates import ANALOGUE_DATES


def file_ok(path):
    path = Path(path)
    return path.is_file() and path.stat().st_size > 0


def validate_final(path):
    if not file_ok(path):
        return False
    try:
        with xr.open_dataset(path, decode_times=False) as ds:
            if not FINAL_VARIABLES.issubset(ds.data_vars):
                return False
            if ds.sizes.get("time") != 1:
                return False
        return True
    except Exception:
        return False


def safe_remove(path):
    path = Path(path)
    if path.exists() or path.is_symlink():
        path.unlink()


def make_paths(rank, analogue_date):
    tag = f"rank{rank:02d}_MJ_{analogue_date.strftime('%Y%m%d')}"
    final = PFILE_BASE / tag
    final.mkdir(parents=True, exist_ok=True)
    return {"tag": tag, "final": final}


def p_file(date, slot, paths):
    return paths["final"] / f"P{date.strftime('%Y%m%d')}_{slot}"


def existing_p_file(date, slot, paths):
    """The readable P file for this slot (single 00/06/12/18 schedule)."""
    candidate = p_file(date, slot, paths)
    if file_ok(candidate):
        return candidate
    raise FileNotFoundError(
        f"No P file for {date.strftime('%Y-%m-%d')} slot {slot}: "
        f"{candidate.name} is not readable in {paths['final']}"
    )


def s_file(date, slot, paths):
    return paths["final"] / f"S{date.strftime('%Y%m%d')}_{slot}"


def s_slot_complete(date, slot, paths):
    return validate_final(s_file(date, slot, paths))


def compute_dp_p_theta(date, slot, paths):
    # Notebook-03 compute_dp_p_theta + hyam/hybm read (folded-in 04_write_coords.sh).
    filename = existing_p_file(date, slot, paths)
    with xr.open_dataset(filename) as ds:
        p_aligned = ds.p.transpose("time", "lev", "lat", "lon")
        p0 = 100000.0
        r_cp = 287.0 / 1004.0
        theta = ds.T * (p0 / p_aligned) ** r_cp
        theta.name = "theta"
        theta.attrs["long_name"] = "Potential Temperature"
        theta.attrs["units"] = "K"
        coeffs = {}
        for name in HYBRID_COEFF_VARS:
            if name in ds:
                coeffs[name] = ds[name].load()
            else:
                print(f"[WARN] {filename.name} has no {name}; S file will omit it")
        return ds.PS.load(), p_aligned.load(), theta.load(), coeffs


def load_clim_with_overwritten_time(date, slot):
    # Notebook-03 load_clim_with_overwritten_time; 6h filename label.
    timestamp = pd.to_datetime(date.strftime("%Y%m%d") + slot, format="%Y%m%d%H")
    date_dashed = timestamp.strftime("%Y-%m-%d")
    filename = CLIM_PATH / f"t_{date_dashed}_6h{slot}_tbar.nc"
    with xr.open_dataset(filename) as ds:
        ds_new = ds.rename({"t": "t_clim"})
        ds_new = ds_new.assign_coords(time=("time", [timestamp]))
        return ds_new.load()


def shift_slot(date, slot, offset_hours):
    base = pd.to_datetime(date.strftime("%Y%m%d") + slot, format="%Y%m%d%H")
    shifted = base + pd.Timedelta(hours=offset_hours)
    return pd.Timestamp(shifted.date()), shifted.strftime("%H")


def load_three_times(date, slot):
    # Neighbours at +/-DT_HOURS (=6 h) instead of +/-3 h.
    date_minus, slot_minus = shift_slot(date, slot, -DT_HOURS)
    date_plus, slot_plus = shift_slot(date, slot, +DT_HOURS)
    t_minus = load_clim_with_overwritten_time(date_minus, slot_minus).t_clim
    t_center = load_clim_with_overwritten_time(date, slot).t_clim
    t_plus = load_clim_with_overwritten_time(date_plus, slot_plus).t_clim
    return t_minus, t_center, t_plus


def compute_dtclim_dt_centered(t_minus, t_center, t_plus):
    dt_seconds = 2 * 3600 * DT_HOURS          # was 2 * 3600 * 3
    dtdt_values = (t_plus.drop_vars("time") - t_minus.drop_vars("time")) / dt_seconds
    dtdt = dtdt_values.assign_coords({"time": t_center.time}).transpose("time", "lev", "lat", "lon")
    dtdt.name = "dTdt"
    dtdt.attrs["units"] = "K s-1"
    dtdt.attrs["long_name"] = "Time derivative of T_clim"
    return dtdt


def save_compiled_dataset(date, slot, paths):
    destination = s_file(date, slot, paths)
    if validate_final(destination):
        return destination

    ps, p, theta, coeffs = compute_dp_p_theta(date, slot, paths)
    t_minus, t_clim, t_plus = load_three_times(date, slot)
    dtdt = compute_dtclim_dt_centered(t_minus, t_clim, t_plus)

    data_vars = {"THETA": theta, "TCLIM": t_clim, "DTDT": dtdt, "PS": ps}
    data_vars.update(coeffs)  # hyam, hybm from the P file

    ds = xr.Dataset(data_vars)
    ds.PS.attrs["long_name"] = "Pressure"; ds.PS.attrs["units"] = "hPa"
    ds.THETA.attrs["long_name"] = "Potential Temperature"; ds.THETA.attrs["units"] = "K"
    ds.DTDT.attrs["long_name"] = "Time derivative of T_clim"; ds.DTDT.attrs["units"] = "K s-1"

    temporary = Path(str(destination) + ".tmp")
    safe_remove(temporary)
    ds.to_netcdf(temporary)
    if not validate_final(temporary):
        raise RuntimeError(f"Invalid staged S file: {temporary}")
    os.replace(temporary, destination)
    print(f"[FINAL] {destination}")
    return destination


def process_one_slot(date, slot, paths):
    if s_slot_complete(date, slot, paths):
        return
    try:
        save_compiled_dataset(date, slot, paths)
    except FileNotFoundError as error:
        print(f"[SKIP] {date.date()} slot {slot}: {error}")
    except Exception as error:
        print(f"[ERROR] {date.date()} slot {slot}: {error}")


def cleanup_stale_temporaries(base_dates, paths):
    for date in base_dates:
        for slot in SLOTS:
            safe_remove(Path(str(s_file(date, slot, paths)) + ".tmp"))


def run_rank(rank, slot_jobs):
    if not 1 <= rank <= len(ANALOGUE_DATES):
        raise ValueError(f"Rank must be 1-{len(ANALOGUE_DATES)}")

    analogue_date = pd.Timestamp(ANALOGUE_DATES[rank - 1])
    analogue_start = analogue_date - pd.Timedelta(days=DAYS_BEFORE)
    analogue_end = analogue_date + pd.Timedelta(days=DAYS_AFTER)
    base_start = analogue_start - pd.Timedelta(days=LEAD_DAYS_BEFORE_BASE)
    base_end = analogue_end
    base_dates = list(pd.date_range(base_start, base_end, freq="D"))
    paths = make_paths(rank, analogue_date)
    cleanup_stale_temporaries(base_dates, paths)

    total = len(base_dates) * len(SLOTS)
    print("=" * 80)
    print(f"Rank: {rank:02d}")
    print(f"Analogue: {analogue_date.date()}")
    print(f"Full LAGRANTO input period: {base_start.date()} to {base_end.date()}")
    print(f"Days: {len(base_dates)}; expected files: {total}")
    print(f"P-file / S-file directory: {paths['final']}")
    print("=" * 80)

    tasks = [(d, slot) for d in base_dates for slot in SLOTS
             if not s_slot_complete(d, slot, paths)]
    if tasks:
        print(f"[INFO] Pending S files: {len(tasks)}/{total}")
        Parallel(n_jobs=slot_jobs, backend="loky", verbose=10)(
            delayed(process_one_slot)(d, slot, paths) for d, slot in tasks
        )
    else:
        print("[SKIP] All final S files already exist")

    present = sum(1 for d in base_dates for slot in SLOTS if s_slot_complete(d, slot, paths))
    # First day's slot 00 needs the previous day's 18 and the last day's slot 18
    # needs the next day's 00 (both outside the computed climatology window), so
    # at most 2 slots per rank may be missing at the +/-6 h edges.
    expected_minimum = total - 2
    print(f"[DONE] Rank {rank:02d}: {present}/{total} S files present "
          f"(minimum required: {expected_minimum})")
    if present < expected_minimum:
        raise RuntimeError(
            f"Rank {rank:02d} incomplete: {present}/{total} S files, need at least "
            f"{expected_minimum}. Check [SKIP]/[ERROR] lines above."
        )
    if present < total:
        print(f"[INFO] {total - present} S files missing at the +-6h window edges (expected)")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rank", type=int, required=True)
    parser.add_argument("--slot-jobs", type=int, default=16)
    parser.add_argument("--skip-rename", action="store_true",
                        help="Deprecated no-op: 6-hourly snapshots need no rename.")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.slot_jobs < 1:
        raise ValueError("slot-jobs must be positive")
    if args.skip_rename:
        print("[INFO] --skip-rename is a no-op in the 6-hourly pipeline (no rename phase).")
    run_rank(args.rank, args.slot_jobs)


if __name__ == "__main__":
    main()
