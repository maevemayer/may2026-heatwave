#!/usr/bin/env python3
"""
Stage 01 - prepare ERA5 P-files (CORRECTED 6-hourly snapshots).

Same math as 01_prepare_ERA5.ipynb EXCEPT:
  * -selhour,0,6,12,18 is applied BEFORE the spectral transform, and
  * the 3-hourly mean (cdo timselmean(3)) is REMOVED.
The 6-hourly instantaneous snapshots are used directly, so the resulting
P-files carry the exact valid times 00/06/12/18 (no +1 h offset -> no rename).

Usage:
    python 01_prepare_era5_6h.py START_DATE END_DATE [N_JOBS]
    e.g. python 01_prepare_era5_6h.py 2026-05-15 2026-06-01 30

Paths, HOURS, LEVELS, grid come from the environment (source config.sh first).
"""
import os
import sys
import glob
import pandas as pd
import xarray as xr
from cdo import Cdo
from joblib import Parallel, delayed

# ------------------------------------------------------------------ config ---
os.environ["PATH"] = os.environ.get("CDO_BIN_DIR", "") + ":" + os.environ["PATH"]

ERA5_ML   = os.environ["ERA5_ML_ET"]          # current-year model level (ET)
ERA5_ML_E5= os.environ["ERA5_ML_E5"]
ERA5_SF   = os.environ["ERA5_SF_ET"]
SCRATCH   = os.path.join(os.environ["SCRATCH_BASE"], "prepare")
FINAL     = os.environ["DATA_DIR"]
HOURS_CSV = os.environ.get("HOURS_CSV", "0,6,12,18")
LEVELS    = os.environ.get("LEVELS", "").strip()
REMAP     = os.environ.get("REMAP_GRID", "r720x360")
LATBOX    = os.environ.get("DOMAIN_LATBOX", "0,360,0,90").strip()   # NH cut from the start
os.makedirs(SCRATCH, exist_ok=True)
os.makedirs(FINAL, exist_ok=True)

era5_vars = {134: "sp", 130: "t", 131: "u", 132: "v", 135: "w"}

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



# ------------------------------------------------------- daily regridding ----
def process_era5_daily(year, month, day, var_num, var, ml="ml"):
    """Regrid one daily ERA5 file to 0.5 deg, keeping ONLY hours 00/06/12/18.

    Identical to the notebook chain (setgridtype -> sp2gp -> remapbil) except:
      * -selhour is prepended (fewer timesteps through the spectral transform)
      * -sellevel is inserted when LEVELS is set (4D vars only)
      * NO timselmean afterwards.
    """
    cdo = Cdo()
    date_str = f"{year}-{month:02d}-{day:02d}"

    if var == "sp":
        data_path = ERA5_SF
        ml = "sf"
    else:
        data_path = ERA5_ML

    var_file = f"{data_path}{var_num}/ET{ml}00_1H_{date_str}_{var_num}.grb"
    if not os.path.exists(var_file):
        print(f"Skipped {var_file}: file does not exist", flush=True)
        return

    reg_file    = os.path.join(SCRATCH, f"reg_{var}_{date_str}.nc")
    gp_file     = os.path.join(SCRATCH, f"gp_{var}_{date_str}.nc")
    regrid_file = os.path.join(SCRATCH, f"{var}_{date_str}_0.5deg_6h.nc")

    if os.path.exists(regrid_file):
        if wind_file_ok(regrid_file, var):
            print(f"{regrid_file} already exists - skipping.", flush=True)
            return regrid_file
        print(f"{regrid_file} is an uncorrected wind file - recomputing.", flush=True)
        os.remove(regrid_file)

    # Step 1: select the 6-hourly snapshots (+ optional levels) FIRST, then
    #         convert reduced->regular Gaussian and GRIB->netCDF.
    sel = f"-selhour,{HOURS_CSV}"
    if LEVELS and var != "sp":
        sel = f"-sellevel,{LEVELS} {sel}"
    cdo.setgridtype("regular", input=f"{sel} {var_file}",
                    output=reg_file, options="-f nc --eccodes")

    # Step 2: spectral -> gridpoint (pass-through if not spectral)
    spectral = False
    try:
        cdo.sp2gp(input=reg_file, output=gp_file, options="-P 16")
        input_for_remap = gp_file
        spectral = True
    except Exception as e:
        print(f"sp2gp failed or not needed for {var}: {e}", flush=True)
        input_for_remap = reg_file

    # Step 3: bilinear remap to 0.5 deg, then cut to the domain box (default NH).
    #         Done in one chained call so nothing outside the box is ever written.
    if LATBOX:
        cdo.sellonlatbox(LATBOX, input=f"-remapbil,{REMAP} {input_for_remap}", output=regrid_file)
    else:
        cdo.remapbil(REMAP, input=input_for_remap, output=regrid_file)

    # (Step 4 timselmean REMOVED - the 4 snapshots are the output.)

    # Step 4 (wind fix): spectral u/v are U*cos(lat) -> divide by cos(lat)
    if var in WIND_VARS:
        if not spectral:
            raise RuntimeError(f"{var} was not spectral - check whether the cos(lat) fix applies to {var_file}")
        divide_wind_by_coslat(regrid_file, var)

    for f in (reg_file, gp_file):
        if os.path.exists(f):
            os.remove(f)

    print(f"Processed {var}: {regrid_file}", flush=True)
    return regrid_file


# --------------------------------------------------- pressure & dp (as-is) ---
def compute_pressure_and_dp(sp_file, t_file, output_dir, date_str):
    """UNCHANGED from the notebook (hybrid p and dp from T-file coefficients)."""
    os.makedirs(output_dir, exist_ok=True)

    ds_t = xr.open_dataset(t_file)
    ds_sp = xr.open_dataset(sp_file)
    sp = ds_sp["sp"]

    hyam = ds_t["hyam"].expand_dims({"time": sp.time, "lat": sp.lat, "lon": sp.lon})
    hybm = ds_t["hybm"].expand_dims({"time": sp.time, "lat": sp.lat, "lon": sp.lon})
    hyai = ds_t["hyai"].expand_dims({"time": sp.time, "lat": sp.lat, "lon": sp.lon})
    hybi = ds_t["hybi"].expand_dims({"time": sp.time, "lat": sp.lat, "lon": sp.lon})

    p = hyam + hybm * sp
    p_half = hyai + hybi * sp
    dp = p_half.isel(nhyi=slice(1, None)) - p_half.isel(nhyi=slice(None, -1))

    p = p.rename({"nhym": "lev"})
    dp = dp.rename({"nhyi": "lev"})
    p = p.assign_coords(lev=ds_t["lev"]); p.name = "p"; p.attrs["units"] = "Pa"
    dp = dp.assign_coords(lev=ds_t["lev"]); dp.name = "dp"; dp.attrs["units"] = "Pa"

    os.makedirs(f"{output_dir}/p", exist_ok=True)
    os.makedirs(f"{output_dir}/dp", exist_ok=True)
    p_file = os.path.join(output_dir, f"p/p_{date_str}_0.5deg_6h.nc")
    dp_file = os.path.join(output_dir, f"dp/dp_{date_str}_0.5deg_6h.nc")
    p.to_dataset().to_netcdf(p_file)
    dp.to_dataset().to_netcdf(dp_file)
    return (p_file, dp_file)


# ------------------------------------------------- merge & split (as-is) -----
def merge_and_split_daily_vars(date_str, var_list=None, input_dir=None, output_dir=None):
    """UNCHANGED math: merge vars, sp Pa->hPa, rename, split by time -> P-files.

    With 6-hourly snapshots the split naturally yields PYYYYMMDD_00/06/12/18
    with correct valid times (no centring offset, so no rename step needed).
    """
    if var_list is None:
        var_list = ["t", "u", "v", "w", "sp", "dp", "p"]
    if input_dir is None:
        input_dir = SCRATCH
    os.makedirs(output_dir, exist_ok=True)

    datasets = []
    for var in var_list:
        file_path = os.path.join(input_dir, var, f"{var}_{date_str}_0.5deg_6h.nc") \
            if var in ("p", "dp") else \
            os.path.join(input_dir, f"{var}_{date_str}_0.5deg_6h.nc")
        if os.path.exists(file_path):
            datasets.append(xr.open_dataset(file_path))
        else:
            print(f"File not found for {var}: {file_path}", flush=True)

    if not datasets:
        print("No variables found to merge!", flush=True)
        return

    merged = xr.merge(datasets)

    if "sp" in merged.data_vars:
        merged["sp"] = merged["sp"] / 100.0
        merged["sp"].attrs["units"] = "hPa"
        merged["sp"].attrs["long_name"] = "Surface Pressure"

    rename_dict = {"u": "U", "v": "V", "w": "OMEGA", "sp": "PS", "t": "T"}
    merged = merged.rename({k: v for k, v in rename_dict.items() if k in merged.data_vars})

    for t_idx in range(len(merged.time)):
        single = merged.isel(time=t_idx)
        time_value = pd.to_datetime(single.time.values)
        time_str = time_value.strftime("P%Y%m%d_%H")
        single = single.expand_dims(time=[merged.time.isel(time=t_idx).values])
        out_file = os.path.join(output_dir, time_str)
        single.to_netcdf(out_file)
        single.close()
        print(f"Saved {out_file}", flush=True)

    merged.close()
    for ds in datasets:
        ds.close()


# --------------------------------------------------------------- driver ------
def process_one_day(date):
    year, month, day = date.year, date.month, date.day
    date_str = f"{year}-{month:02d}-{day:02d}"
    print(f"\n=== Processing {date_str} ===", flush=True)

    regrid_files = {}
    for var_num, var_short in era5_vars.items():
        f = process_era5_daily(year, month, day, var_num, var_short)
        if f:
            regrid_files[var_short] = f

    if "t" in regrid_files and "sp" in regrid_files:
        compute_pressure_and_dp(regrid_files["sp"], regrid_files["t"], SCRATCH, date_str)
    else:
        print("Skipping p/dp computation (missing t or sp).", flush=True)

    merge_and_split_daily_vars(date_str, input_dir=SCRATCH, output_dir=FINAL)


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: 01_prepare_era5_6h.py START_DATE END_DATE [N_JOBS]")
    start_date, end_date = sys.argv[1], sys.argv[2]
    n_jobs = int(sys.argv[3]) if len(sys.argv) > 3 else 30
    dates = pd.date_range(start_date, end_date, freq="D")
    Parallel(n_jobs=n_jobs, backend="loky", verbose=10)(
        delayed(process_one_day)(d) for d in dates
    )


if __name__ == "__main__":
    main()
