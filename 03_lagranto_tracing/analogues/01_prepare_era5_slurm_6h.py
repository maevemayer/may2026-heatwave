#!/usr/bin/env python3
"""Build LAGRANTO P files for the ranked MJ analogue dates - CORRECTED 6-hourly.

Faithful port of prepare_era5_top20_slurm.py. The ONLY changes are the ones
agreed for the 6-hourly correction:
  * SLOTS are 00/06/12/18 (was 01/04/.../22),
  * -selhour,0,6,12,18 (+ optional -sellevel) is applied BEFORE the spectral
    transform in process_era5_daily,
  * the 3-hourly mean (cdo.timselmean(3)) is REMOVED - the 6-hourly
    instantaneous snapshots are used directly,
  * intermediate suffix _0.5deg_3h.nc -> _0.5deg_6h.nc.
Everything else (validation, staging, scratch cleanup, hybrid p/dp, merge/split)
is identical. Because snapshots are instantaneous, the P files already carry the
exact 00/06/12/18 valid times (no +1 h offset -> the S-file stage does not
rename anymore).

LEVELS env var (default empty = all 137 levels; math identical to notebook):
  set e.g. LEVELS=70/137 to keep the lowest 68 levels. WARNING: that clips the
  200 hPa trajectory starts and any air above ~250 hPa - see README.
"""

import argparse
import os
import shutil
from pathlib import Path

import pandas as pd
import xarray as xr
from cdo import Cdo
from joblib import Parallel, delayed


# Paths come from config_analogues.sh (env); defaults = the RERUN tree.
_DEF_FINAL = "${LAGRANTO_RERUN}/analogues_windfix"
_DEF_SCRATCH = "${SCRATCH}/lagranto/RERUN_2026_MAY_6h_windfix"
FINAL_BASE = Path(os.environ.get("ANALOGUE_FINAL_BASE", _DEF_FINAL))
SCRATCH_BASE = Path(os.environ.get("ANALOGUE_SCRATCH_BASE", _DEF_SCRATCH)) / "prepare_era5_slurm"

DATA_PATH_ML_E5 = Path(os.environ.get("ERA5_ML_E5", "/pool/data/ERA5/E5/ml/an/1H"))
DATA_PATH_ML_ET = Path(os.environ.get("ERA5_ML_ET", "/pool/data/ERA5/ET/ml/an/1H"))
DATA_PATH_SF_E5 = Path(os.environ.get("ERA5_SF_E5", "/pool/data/ERA5/E5/sf/an/1H"))
DATA_PATH_SF_ET = Path(os.environ.get("ERA5_SF_ET", "/pool/data/ERA5/ET/sf/an/1H"))
CDO_PATH = Path(os.environ.get("CDO_BIN_DIR", "/path/to/cdo/bin"))

ERA5_VARS = {134: "sp", 130: "t", 131: "u", 132: "v", 135: "w"}

# ---------------------------------------------------------------------------
# Wind fix (Oct 2026): ERA5 spectral u/v (params 131/132) are stored as
# U*cos(lat) and V*cos(lat); cdo sp2gp does not divide that out, so the old
# P-files had winds too weak by cos(lat) (0.71 at 45N). After the transform and
# remap, U and V are divided by cos(lat) of the target grid; verified against
# independent MARS 0.5 deg winds (corr 1.000, slope 0.999). Corrected files
# carry the attribute COSLAT_FLAG, and any u/v file without it is recomputed.
# ---------------------------------------------------------------------------
COSLAT_FLAG = "u_v_divided_by_coslat"
WIND_VARS = ("u", "v")


def wind_file_ok(path, var):
    """True unless this is a u/v file that has not been cos(lat)-corrected."""
    if var not in WIND_VARS:
        return True
    try:
        with xr.open_dataset(path) as ds:
            return ds[var].attrs.get(COSLAT_FLAG) == "yes"
    except Exception:
        return False


def divide_wind_by_coslat(path, var):
    """In place: var /= cos(lat) on the (regular) target grid."""
    import numpy as np
    with xr.open_dataset(path) as ds:
        ds = ds.load()
    attrs, dtype = dict(ds[var].attrs), ds[var].dtype
    coslat = np.cos(np.deg2rad(ds["lat"].astype("float64")))
    ds[var] = (ds[var] / coslat).astype(dtype)
    ds[var].attrs.update(attrs)
    ds[var].attrs[COSLAT_FLAG] = "yes"
    tmp = str(path) + ".coslat.tmp"
    ds.to_netcdf(tmp)
    os.replace(tmp, path)

SLOTS = ["00", "06", "12", "18"]                 # was 01/04/07/10/13/16/19/22
HOURS_CSV = "0,6,12,18"
LEVELS = os.environ.get("LEVELS", "").strip()    # empty -> all levels
LATBOX = os.environ.get("DOMAIN_LATBOX", "0,360,0,90").strip()   # NH cut from the start
REMAP = os.environ.get("REMAP_GRID", "r720x360")
FINAL_VARIABLES = {"U", "V", "OMEGA", "PS", "T", "p", "dp"}

# Keep synchronized with t_clim / s_files / lagranto scripts.
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
    root = SCRATCH_BASE / tag
    work = root / "daily_0.5deg_6h"
    final = FINAL_BASE / tag
    work.mkdir(parents=True, exist_ok=True)
    final.mkdir(parents=True, exist_ok=True)
    return {"tag": tag, "root": root, "work": work, "final": final}


def variable_file(date, variable, paths):
    return paths["work"] / variable / f"{variable}_{date.strftime('%Y-%m-%d')}_0.5deg_6h.nc"


def final_file(date, slot, paths):
    return paths["final"] / f"P{date.strftime('%Y%m%d')}_{slot}"


def expected_final_files(base_dates, paths):
    return [final_file(date, slot, paths) for date in base_dates for slot in SLOTS]


def day_final_complete(date, paths):
    return all(validate_final(final_file(date, slot, paths)) for slot in SLOTS)


def era5_file_for_date(date, var_num, variable):
    date_string = date.strftime("%Y-%m-%d")
    if variable == "sp":
        candidates = [("ET", DATA_PATH_SF_ET, "sf"), ("E5", DATA_PATH_SF_E5, "sf")]
    else:
        candidates = [("ET", DATA_PATH_ML_ET, "ml"), ("E5", DATA_PATH_ML_E5, "ml")]
    for prefix, base, stream in candidates:
        source = base / str(var_num) / f"{prefix}{stream}00_1H_{date_string}_{var_num}.grb"
        if source.is_file():
            return source
    return None


def process_era5_daily(date, var_num, variable, paths, cdo_threads):
    output = variable_file(date, variable, paths)
    if file_ok(output) and wind_file_ok(output, variable):
        return output

    output.parent.mkdir(parents=True, exist_ok=True)
    date_string = date.strftime("%Y-%m-%d")
    source = era5_file_for_date(date, var_num, variable)
    if source is None:
        raise FileNotFoundError(
            f"No ERA5/ERA5T source file for {date_string} variable {variable}"
        )

    reg_file = output.parent / f"reg_{variable}_{date_string}.nc"
    gp_file = output.parent / f"gp_{variable}_{date_string}.nc"
    regrid_file = output.parent / f"{variable}_{date_string}_0.5deg.nc"
    temporary = Path(str(output) + ".tmp")
    for path in [reg_file, gp_file, regrid_file, temporary]:
        safe_remove(path)

    cdo = Cdo()
    # 6h correction: select the snapshots (+ optional levels) BEFORE the
    # spectral transform, and drop the 3-hourly mean entirely.
    sel = f"-selhour,{HOURS_CSV}"
    if LEVELS and variable != "sp":
        sel = f"-sellevel,{LEVELS} {sel}"
    cdo.setgridtype(
        "regular",
        input=f"{sel} {source}",
        output=str(reg_file),
        options="-f nc --eccodes",
    )
    spectral = False
    try:
        cdo.sp2gp(input=str(reg_file), output=str(gp_file), options=f"-P {cdo_threads}")
        remap_input = gp_file
        spectral = True
    except Exception as error:
        print(f"[WARN] sp2gp failed or unnecessary for {date_string} {variable}: {error}")
        remap_input = reg_file

    # remap straight to the staged output (no timselmean), then cut to the
    # domain box (default NH) in the same chained call.
    if LATBOX:
        cdo.sellonlatbox(LATBOX, input=f"-remapbil,{REMAP} {remap_input}", output=str(temporary))
    else:
        cdo.remapbil(REMAP, input=str(remap_input), output=str(temporary))
    # wind fix: spectral u/v are U*cos(lat) -> divide by cos(lat)
    if variable in WIND_VARS:
        if not spectral:
            raise RuntimeError(f"{variable} was not spectral - check whether the cos(lat) fix applies to {source}")
        divide_wind_by_coslat(temporary, variable)
    os.replace(temporary, output)

    for path in [reg_file, gp_file, regrid_file]:
        safe_remove(path)
    return output


def compute_pressure_and_dp(date, paths):
    """UNCHANGED from the original (hybrid p and dp)."""
    sp_file = variable_file(date, "sp", paths)
    t_file = variable_file(date, "t", paths)
    p_file = variable_file(date, "p", paths)
    dp_file = variable_file(date, "dp", paths)
    if file_ok(p_file) and file_ok(dp_file):
        return p_file, dp_file

    p_file.parent.mkdir(parents=True, exist_ok=True)
    dp_file.parent.mkdir(parents=True, exist_ok=True)

    with xr.open_dataset(t_file) as ds_t, xr.open_dataset(sp_file) as ds_sp:
        sp = ds_sp["sp"]
        hyam = ds_t["hyam"].expand_dims({"time": sp.time, "lat": sp.lat, "lon": sp.lon})
        hybm = ds_t["hybm"].expand_dims({"time": sp.time, "lat": sp.lat, "lon": sp.lon})
        hyai = ds_t["hyai"].expand_dims({"time": sp.time, "lat": sp.lat, "lon": sp.lon})
        hybi = ds_t["hybi"].expand_dims({"time": sp.time, "lat": sp.lat, "lon": sp.lon})

        p = hyam + hybm * sp
        p_half = hyai + hybi * sp
        dp = p_half.diff("nhyi")

        p = p.rename({"nhym": "lev"})
        dp = dp.rename({"nhyi": "lev"})
        p = p.assign_coords(lev=ds_t["lev"])
        dp = dp.assign_coords(lev=ds_t["lev"])
        p.name = "p"; p.attrs["units"] = "Pa"
        dp.name = "dp"; dp.attrs["units"] = "Pa"

        temporary_p = Path(str(p_file) + ".tmp")
        temporary_dp = Path(str(dp_file) + ".tmp")
        safe_remove(temporary_p); safe_remove(temporary_dp)
        p.to_dataset().to_netcdf(temporary_p)
        dp.to_dataset().to_netcdf(temporary_dp)
        os.replace(temporary_p, p_file)
        os.replace(temporary_dp, dp_file)
    return p_file, dp_file


def merge_and_split_day(date, paths):
    """UNCHANGED math; time-count check now expects len(SLOTS)=4 snapshots."""
    variables = ["t", "u", "v", "w", "sp", "dp", "p"]
    datasets = []
    try:
        missing = []
        for variable in variables:
            path = variable_file(date, variable, paths)
            if not file_ok(path):
                missing.append(variable)
            else:
                datasets.append(xr.open_dataset(path))
        if missing:
            raise FileNotFoundError(f"Missing variables for {date.date()}: {missing}")

        merged = xr.merge(datasets)
        merged["sp"] = merged["sp"] / 100.0
        merged["sp"].attrs["units"] = "hPa"
        merged["sp"].attrs["long_name"] = "Surface Pressure"
        merged = merged.rename({"u": "U", "v": "V", "w": "OMEGA", "sp": "PS", "t": "T"})

        if merged.sizes.get("time") != len(SLOTS):
            raise ValueError(
                f"Expected {len(SLOTS)} 6-hourly times for {date.date()}, "
                f"found {merged.sizes.get('time')}"
            )

        for time_index, expected_slot in enumerate(SLOTS):
            single = merged.isel(time=time_index)
            time_value = pd.to_datetime(single.time.values)
            actual_date = pd.Timestamp(time_value).normalize()
            actual_slot = pd.Timestamp(time_value).strftime("%H")
            if actual_date != date.normalize() or actual_slot != expected_slot:
                raise ValueError(
                    f"Unexpected output time {time_value}; expected "
                    f"{date.strftime('%Y-%m-%d')} {expected_slot}:00"
                )

            destination = final_file(date, expected_slot, paths)
            if validate_final(destination):
                continue
            safe_remove(destination)
            temporary = Path(str(destination) + ".tmp")
            safe_remove(temporary)
            time_coordinate = merged.time.isel(time=time_index).values
            single = single.expand_dims(time=[time_coordinate])
            single.to_netcdf(temporary)
            if not validate_final(temporary):
                raise RuntimeError(f"Invalid staged P file: {temporary}")
            os.replace(temporary, destination)
            print(f"[FINAL] {destination}")
        merged.close()
    finally:
        for dataset in datasets:
            dataset.close()


def cleanup_day_work(date, paths):
    for variable in ["t", "u", "v", "w", "sp", "p", "dp"]:
        safe_remove(variable_file(date, variable, paths))


def process_one_day(date, paths, cdo_threads):
    if day_final_complete(date, paths):
        print(f"[SKIP] {date.date()} final P files already complete")
        cleanup_day_work(date, paths)
        return
    print(f"[DAY] Processing {date.date()}")
    for var_num, variable in ERA5_VARS.items():
        process_era5_daily(date, var_num, variable, paths, cdo_threads)
    compute_pressure_and_dp(date, paths)
    merge_and_split_day(date, paths)
    if not day_final_complete(date, paths):
        raise RuntimeError(f"Final P files incomplete for {date.date()}")
    cleanup_day_work(date, paths)


def cleanup_successful_scratch(paths):
    root = paths["root"].resolve()
    if SCRATCH_BASE.resolve() not in root.parents or not root.name.startswith("rank"):
        raise RuntimeError(f"Refusing to delete suspicious scratch path: {root}")
    if root.exists():
        shutil.rmtree(root)
        print(f"[CLEAN] Removed successful scratch tree: {root}")


def cleanup_stale_temporaries(base_dates, paths):
    for date in base_dates:
        for slot in SLOTS:
            safe_remove(Path(str(final_file(date, slot, paths)) + ".tmp"))


def run_rank(rank, day_jobs, cdo_threads):
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

    print("=" * 80)
    print(f"Rank: {rank:02d}")
    print(f"Analogue: {analogue_date.date()}")
    print(f"Analogue window: {analogue_start.date()} to {analogue_end.date()}")
    print(f"Full LAGRANTO input period: {base_start.date()} to {base_end.date()}")
    print(f"Days: {len(base_dates)}; expected files: {len(base_dates) * len(SLOTS)}")
    print(f"LEVELS: {LEVELS or 'ALL (137)'}")
    print(f"Scratch: {paths['root']}")
    print(f"Final directory: {paths['final']}")
    print("=" * 80)

    pending_dates = [d for d in base_dates if not day_final_complete(d, paths)]
    if not pending_dates:
        print("[SKIP] All final P files already exist")
        cleanup_successful_scratch(paths)
        return

    print(f"[INFO] Pending days: {len(pending_dates)}/{len(base_dates)}")
    try:
        Parallel(n_jobs=day_jobs, backend="loky", verbose=10)(
            delayed(process_one_day)(d, paths, cdo_threads) for d in pending_dates
        )
        incomplete = [p for p in expected_final_files(base_dates, paths) if not validate_final(p)]
        if incomplete:
            raise RuntimeError(f"Final validation failed for {len(incomplete)} P files")
        cleanup_successful_scratch(paths)
        print(f"[DONE] Rank {rank:02d} completed successfully")
    except Exception:
        print(f"[RETAIN] Scratch retained for restart: {paths['root']}")
        raise


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rank", type=int, required=True)
    parser.add_argument("--day-jobs", type=int, default=4)
    parser.add_argument("--cdo-threads", type=int, default=8)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.day_jobs < 1 or args.cdo_threads < 1:
        raise ValueError("day-jobs and cdo-threads must be positive")
    os.environ["PATH"] = str(CDO_PATH) + os.pathsep + os.environ.get("PATH", "")
    print(Cdo().version())
    run_rank(args.rank, args.day_jobs, args.cdo_threads)


if __name__ == "__main__":
    main()
