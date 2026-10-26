#!/usr/bin/env python3
"""9-year T climatologies for ranked MJ analogues - CORRECTED 6-hourly.

Faithful port of t_clim_top20_slurm.py. ONLY the agreed 6-hourly changes:
  * SLOTS are 00/06/12/18 (was 01/04/.../22),
  * -selhour,0,6,12,18 (+ optional -sellevel) applied BEFORE the spectral
    transform in process_era5_daily,
  * the 3-hourly mean (cdo.timselmean(3)) is REMOVED (instantaneous snapshots),
  * every "_3h"/"3h" label in intermediate and final filenames becomes "_6h"/"6h",
    so the tbar files are  t_YYYY-MM-DD_6hHH_tbar.nc  (read by s_files_*_6h.py).
The centered 9-year window selection (analogue_year +/-4, clamped to 1940-2025),
the 21-day running mean, chunked cat/timmean, and the final 9-year timmean are
all UNCHANGED.
"""

import argparse
import glob
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
FINAL_PATH = Path(os.environ.get("ANALOGUE_FINAL_BASE", _DEF_FINAL)) / "t_clim"
SCRATCH_BASE = Path(os.environ.get("ANALOGUE_SCRATCH_BASE", _DEF_SCRATCH)) / "t_clim_9year_slurm"

DATA_PATH_ML_E5 = Path(os.environ.get("ERA5_ML_E5", "/pool/data/ERA5/E5/ml/an/1H"))
DATA_PATH_ML_ET = Path(os.environ.get("ERA5_ML_ET", "/pool/data/ERA5/ET/ml/an/1H"))
CDO_PATH = Path(os.environ.get("CDO_BIN_DIR", "/path/to/cdo/bin"))

VAR_NUM = 130
VAR_NAME = "t"
SLOTS = ["00", "06", "12", "18"]                 # was 01/04/07/10/13/16/19/22
HOURS_CSV = "0,6,12,18"
LEVELS = os.environ.get("LEVELS", "").strip()    # empty -> all levels
LATBOX = os.environ.get("DOMAIN_LATBOX", "0,360,0,90").strip()   # NH cut from the start
REMAP = os.environ.get("REMAP_GRID", "r720x360")

N_CLIMATOLOGY_YEARS = 9
ERA5_FIRST_YEAR = 1940
ERA5_LAST_YEAR = 2025
DAYS_BEFORE = int(os.environ.get("ANALOGUE_DAYS_BEFORE", "4"))   # days before the analogue date
DAYS_AFTER = int(os.environ.get("ANALOGUE_DAYS_AFTER", "3"))     # days after the analogue date
ANALOGUE_WINDOW_DAYS = DAYS_BEFORE + DAYS_AFTER + 1              # total span (for logging/validation)
RUNNING_WINDOW_DAYS = 21
LEAD_DAYS_BEFORE_BASE = int(os.environ.get("ANALOGUE_LEAD_DAYS", "15"))   # back-trajectory length (days of data before the first arrival)

# Analogue dates come from the shared analogue_dates.txt (edit that ONE file).
from analogue_dates import ANALOGUE_DATES


def file_ok(path):
    path = Path(path)
    return path.is_file() and path.stat().st_size > 0


def validate_netcdf(path):
    if not file_ok(path):
        return False
    try:
        with xr.open_dataset(path) as ds:
            if not ds.data_vars:
                return False
            if "time" in ds.sizes and ds.sizes["time"] != 1:
                return False
        return True
    except Exception:
        return False


def safe_remove(path):
    path = Path(path)
    if path.exists() or path.is_symlink():
        path.unlink()


def same_month_day_in_year(date, year):
    return pd.Timestamp(year=year, month=date.month, day=date.day)


def climatology_years_for_analogue(analogue_year, n_years=N_CLIMATOLOGY_YEARS,
                                   first_year=ERA5_FIRST_YEAR, last_year=ERA5_LAST_YEAR):
    """UNCHANGED: centered n-year window, clamped to the ERA5 range."""
    if n_years % 2 == 0:
        raise ValueError("The climatology window must contain an odd number of years")
    half = n_years // 2
    start_year = analogue_year - half
    end_year = analogue_year + half
    if start_year < first_year:
        start_year = first_year
        end_year = first_year + n_years - 1
    elif end_year > last_year:
        end_year = last_year
        start_year = last_year - n_years + 1
    years = list(range(start_year, end_year + 1))
    if len(years) != n_years or analogue_year not in years:
        raise RuntimeError(f"Invalid climatology years for {analogue_year}: {years}")
    return years


def make_paths(rank, analogue_date):
    tag = f"rank{rank:02d}_MJ_{analogue_date.strftime('%Y%m%d')}"
    root = SCRATCH_BASE / tag
    paths = {
        "root": root,
        "dump": root / "tbar_dump",
        "temp": root / "tbar_dump" / "tmp",
        "intermediate": root / "bar_intermediate",
        "stage": root / "final_stage",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def final_file(date, slot):
    return FINAL_PATH / f"{VAR_NAME}_{date.strftime('%Y-%m-%d')}_6h{slot}_tbar.nc"


def expected_final_files(base_dates):
    return [final_file(date, slot) for date in base_dates for slot in SLOTS]


def date_final_complete(date):
    return all(validate_netcdf(final_file(date, slot)) for slot in SLOTS)


def split_file(date, slot, paths):
    return paths["dump"] / f"{VAR_NAME}_{date.strftime('%Y-%m-%d')}_6h{slot}.nc"


def split_day_complete(date, paths):
    return all(file_ok(split_file(date, slot, paths)) for slot in SLOTS)


def running_file(date, slot, paths):
    return paths["intermediate"] / f"{VAR_NAME}_{date.strftime('%Y-%m-%d')}_6h{slot}_mean.nc"


def running_date_complete(date, paths):
    return all(file_ok(running_file(date, slot, paths)) for slot in SLOTS)


def era5_file_for_date(date):
    date_string = date.strftime("%Y-%m-%d")
    et_file = DATA_PATH_ML_ET / str(VAR_NUM) / f"ETml00_1H_{date_string}_{VAR_NUM}.grb"
    e5_file = DATA_PATH_ML_E5 / str(VAR_NUM) / f"E5ml00_1H_{date_string}_{VAR_NUM}.grb"
    if et_file.is_file():
        return et_file
    if e5_file.is_file():
        return e5_file
    return None


def process_era5_daily(date, paths, cdo_threads):
    if split_day_complete(date, paths):
        return None

    cdo = Cdo()
    date_string = date.strftime("%Y-%m-%d")
    reg_file = paths["dump"] / f"reg_{VAR_NAME}_{date_string}.nc"
    gp_file = paths["dump"] / f"gp_{VAR_NAME}_{date_string}.nc"
    regrid_file = paths["dump"] / f"{VAR_NAME}_{date_string}_0.5deg.nc"
    regrid_6h = paths["dump"] / f"{VAR_NAME}_{date_string}_0.5deg_6h.nc"

    if file_ok(regrid_6h):
        return regrid_6h

    source_file = era5_file_for_date(date)
    if source_file is None:
        raise FileNotFoundError(f"No ERA5/ERA5T model-level file for {date_string}")

    for path in [reg_file, gp_file, regrid_file]:
        safe_remove(path)

    # 6h correction: select snapshots (+ optional levels) before the transform.
    sel = f"-selhour,{HOURS_CSV}"
    if LEVELS:
        sel = f"-sellevel,{LEVELS} {sel}"
    cdo.setgridtype("regular", input=f"{sel} {source_file}",
                    output=str(reg_file), options="-f nc --eccodes")
    try:
        cdo.sp2gp(input=str(reg_file), output=str(gp_file), options=f"-P {cdo_threads}")
        remap_input = gp_file
    except Exception as error:
        print(f"[WARN] sp2gp failed or was unnecessary for {date_string}: {error}")
        remap_input = reg_file

    # remap straight to the staged 6h file (no timselmean), then cut to the
    # domain box (default NH) in the same chained call.
    temporary = Path(str(regrid_6h) + ".tmp")
    safe_remove(temporary)
    if LATBOX:
        cdo.sellonlatbox(LATBOX, input=f"-remapbil,{REMAP} {remap_input}", output=str(temporary))
    else:
        cdo.remapbil(REMAP, input=str(remap_input), output=str(temporary))
    os.replace(temporary, regrid_6h)

    for path in [reg_file, gp_file, regrid_file]:
        safe_remove(path)
    return regrid_6h


def split_into_6hourly_files(regrid_6h, paths):
    if regrid_6h is None:
        return
    date_string = Path(regrid_6h).name.split(f"{VAR_NAME}_", 1)[1].split("_0.5deg_6h.nc", 1)[0]
    date = pd.Timestamp(date_string)
    with xr.open_dataset(regrid_6h) as ds:
        if ds.sizes.get("time", 0) != len(SLOTS):
            raise ValueError(
                f"Expected {len(SLOTS)} 6-hourly times in {regrid_6h}, "
                f"found {ds.sizes.get('time', 0)}"
            )
        for time_index, slot in enumerate(SLOTS):
            output = split_file(date, slot, paths)
            if file_ok(output):
                continue
            temporary = Path(str(output) + ".tmp")
            safe_remove(temporary)
            single = ds.isel(time=time_index)
            single = single.expand_dims(time=[single.time.values])
            single.to_netcdf(temporary)
            single.close()
            os.replace(temporary, output)
    safe_remove(regrid_6h)


def process_extended_period(target_dates, paths, preprocess_jobs, cdo_threads):
    first = min(target_dates) - pd.Timedelta(days=11)
    last = max(target_dates) + pd.Timedelta(days=11)
    extended_dates = pd.date_range(first, last, freq="D")
    missing_dates = [d for d in extended_dates if not split_day_complete(d, paths)]
    print(f"[CACHE] Split daily files: "
          f"{len(extended_dates) - len(missing_dates)}/{len(extended_dates)} complete")

    daily_files = Parallel(n_jobs=preprocess_jobs, verbose=10)(
        delayed(process_era5_daily)(d, paths, cdo_threads) for d in missing_dates
    )
    daily_files = [p for p in daily_files if p is not None]
    Parallel(n_jobs=preprocess_jobs, verbose=10)(
        delayed(split_into_6hourly_files)(p, paths) for p in daily_files
    )
    incomplete = [d for d in extended_dates if not split_day_complete(d, paths)]
    if incomplete:
        raise RuntimeError(f"Incomplete split dates: {incomplete}")


def compute_running_mean(date, slot, paths):
    """UNCHANGED 21-day running mean per slot (chunked cat/timmean)."""
    output = running_file(date, slot, paths)
    if file_ok(output):
        return output
    half_window = RUNNING_WINDOW_DAYS // 2
    window = pd.date_range(date - pd.Timedelta(days=half_window),
                           date + pd.Timedelta(days=half_window), freq="D")
    inputs = [split_file(day, slot, paths) for day in window]
    missing = [p for p in inputs if not file_ok(p)]
    if missing:
        raise FileNotFoundError(f"Missing {len(missing)} inputs for {date.date()} slot {slot}")

    cdo = Cdo()
    chunk_means = []
    try:
        for chunk_index in range(0, len(inputs), 3):
            chunk = inputs[chunk_index: chunk_index + 3]
            chunk_file = paths["temp"] / f"chunk_{chunk_index}_{date.strftime('%Y-%m-%d')}_{slot}.nc"
            chunk_mean = Path(str(chunk_file).replace(".nc", "_mean.nc"))
            safe_remove(chunk_file); safe_remove(chunk_mean)
            cdo.cat(input=" ".join(map(str, chunk)), output=str(chunk_file))
            cdo.timmean(input=str(chunk_file), output=str(chunk_mean))
            safe_remove(chunk_file)
            chunk_means.append(chunk_mean)
        concatenated = paths["temp"] / f"tmp_{slot}_{date.strftime('%Y-%m-%d')}.nc"
        temporary = Path(str(output) + ".tmp")
        safe_remove(concatenated); safe_remove(temporary)
        cdo.cat(input=" ".join(map(str, chunk_means)), output=str(concatenated))
        cdo.timmean(input=str(concatenated), output=str(temporary))
        os.replace(temporary, output)
        safe_remove(concatenated)
    finally:
        for p in chunk_means:
            safe_remove(p)
    return output


def cleanup_year_dump(paths):
    for pattern in ["*.nc", "*.tmp"]:
        for filename in glob.glob(str(paths["dump"] / pattern)):
            safe_remove(filename)
    paths["temp"].mkdir(parents=True, exist_ok=True)


def process_climatology_year(year, target_dates, paths, preprocess_jobs, light_jobs, cdo_threads):
    year_dates = [same_month_day_in_year(d, year) for d in target_dates]
    if all(running_date_complete(d, paths) for d in year_dates):
        print(f"[SKIP] Running means already complete for {year}")
        cleanup_year_dump(paths)
        return
    print(f"[YEAR] Preparing climatology year {year}")
    process_extended_period(year_dates, paths, preprocess_jobs, cdo_threads)
    tasks = [(d, slot) for d in year_dates for slot in SLOTS
             if not file_ok(running_file(d, slot, paths))]
    Parallel(n_jobs=light_jobs, verbose=10)(
        delayed(compute_running_mean)(d, slot, paths) for d, slot in tasks
    )
    incomplete = [running_file(d, slot, paths) for d in year_dates for slot in SLOTS
                  if not file_ok(running_file(d, slot, paths))]
    if incomplete:
        raise RuntimeError(f"Running means incomplete for {year}: {len(incomplete)} files")
    cleanup_year_dump(paths)


def publish_final(date, slot, climatology_years, paths):
    destination = final_file(date, slot)
    if validate_netcdf(destination):
        return destination
    inputs = [running_file(same_month_day_in_year(date, year), slot, paths)
              for year in climatology_years]
    missing = [p for p in inputs if not file_ok(p)]
    if missing:
        raise FileNotFoundError(f"Missing {len(missing)}/{len(inputs)} inputs for {destination.name}")

    cdo = Cdo()
    merged = paths["stage"] / f"merged_{date.strftime('%Y-%m-%d')}_{slot}.nc"
    staged = paths["stage"] / destination.name
    safe_remove(merged); safe_remove(staged)
    cdo.cat(input=" ".join(map(str, inputs)), output=str(merged))
    cdo.timmean(input=str(merged), output=str(staged))
    safe_remove(merged)
    if not validate_netcdf(staged):
        raise RuntimeError(f"Invalid staged climatology: {staged}")

    FINAL_PATH.mkdir(parents=True, exist_ok=True)
    temporary = FINAL_PATH / f".{destination.name}.{os.getpid()}.tmp"
    safe_remove(temporary)
    try:
        shutil.copy2(staged, temporary)
        if not validate_netcdf(temporary):
            raise RuntimeError(f"Invalid copied climatology: {temporary}")
        os.replace(temporary, destination)
    finally:
        safe_remove(temporary)
    safe_remove(staged)
    print(f"[FINAL] {destination}")
    return destination


def cleanup_ephemeral(paths):
    for directory in [paths["temp"], paths["stage"]]:
        if directory.exists():
            for path in directory.iterdir():
                if path.is_file() or path.is_symlink():
                    safe_remove(path)


def cleanup_stale_work_temporaries(base_dates):
    for date in base_dates:
        for slot in SLOTS:
            destination = final_file(date, slot)
            for filename in glob.glob(str(FINAL_PATH / f".{destination.name}.*.tmp")):
                safe_remove(filename)


def cleanup_successful_scratch(paths):
    root = paths["root"].resolve()
    if SCRATCH_BASE.resolve() not in root.parents or not root.name.startswith("rank"):
        raise RuntimeError(f"Refusing to delete suspicious scratch path: {root}")
    shutil.rmtree(root)
    print(f"[CLEAN] Removed successful scratch tree: {root}")


def run_rank(rank, preprocess_jobs, light_jobs, cdo_threads):
    if not 1 <= rank <= len(ANALOGUE_DATES):
        raise ValueError(f"Rank must be 1-{len(ANALOGUE_DATES)}")

    analogue_date = pd.Timestamp(ANALOGUE_DATES[rank - 1])
    analogue_start = analogue_date - pd.Timedelta(days=DAYS_BEFORE)
    analogue_end = analogue_date + pd.Timedelta(days=DAYS_AFTER)
    base_start = analogue_start - pd.Timedelta(days=LEAD_DAYS_BEFORE_BASE)
    base_end = analogue_end
    base_dates = list(pd.date_range(base_start, base_end, freq="D"))
    climatology_years = climatology_years_for_analogue(analogue_date.year)

    FINAL_PATH.mkdir(parents=True, exist_ok=True)
    cleanup_stale_work_temporaries(base_dates)
    paths = make_paths(rank, analogue_date)

    print("=" * 80)
    print(f"Rank: {rank:02d}")
    print(f"Analogue: {analogue_date.date()}")
    print(f"Analogue window: {analogue_start.date()} to {analogue_end.date()}")
    print(f"Final period: {base_start.date()} to {base_end.date()}")
    print(f"Nine-year climatology: {climatology_years[0]}-{climatology_years[-1]}")
    print(f"LEVELS: {LEVELS or 'ALL (137)'}")
    print(f"Scratch: {paths['root']}")
    print(f"Final directory: {FINAL_PATH}")
    print("=" * 80)

    pending_dates = [d for d in base_dates if not date_final_complete(d)]
    if not pending_dates:
        print("[SKIP] All final climatology files already exist")
        cleanup_successful_scratch(paths)
        return

    print(f"[INFO] Pending target dates: {len(pending_dates)}/{len(base_dates)} "
          f"({len(pending_dates) * len(SLOTS)} files)")
    try:
        for year in climatology_years:
            process_climatology_year(year, pending_dates, paths,
                                     preprocess_jobs, light_jobs, cdo_threads)
        Parallel(n_jobs=light_jobs, verbose=10)(
            delayed(publish_final)(date, slot, climatology_years, paths)
            for date in pending_dates for slot in SLOTS
        )
        incomplete = [p for p in expected_final_files(base_dates) if not validate_netcdf(p)]
        if incomplete:
            raise RuntimeError(f"Final validation failed for {len(incomplete)} files")
        cleanup_successful_scratch(paths)
        print(f"[DONE] Rank {rank:02d} completed successfully")
    except Exception:
        cleanup_ephemeral(paths)
        print(f"[RETAIN] Scratch retained for restart after failure: {paths['root']}")
        raise


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rank", type=int, required=True)
    parser.add_argument("--preprocess-jobs", type=int, default=8)
    parser.add_argument("--light-jobs", type=int, default=28)
    parser.add_argument("--cdo-threads", type=int, default=8)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.preprocess_jobs < 1 or args.light_jobs < 1 or args.cdo_threads < 1:
        raise ValueError("Worker counts and cdo-threads must be positive")
    os.environ["PATH"] = str(CDO_PATH) + os.pathsep + os.environ.get("PATH", "")
    print(Cdo().version())
    run_rank(args.rank, args.preprocess_jobs, args.light_jobs, args.cdo_threads)


if __name__ == "__main__":
    main()
