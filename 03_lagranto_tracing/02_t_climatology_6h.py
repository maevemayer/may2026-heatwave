#!/usr/bin/env python3
"""
Stage 02 - 9-year T climatology (tbar), CORRECTED 6-hourly snapshots.

Same math as 02_S_files_T_climatology.ipynb EXCEPT:
  * -selhour,0,6,12,18 applied BEFORE the spectral transform,
  * the 3-hourly mean (timselmean) is REMOVED (instantaneous snapshots),
  * slots are 00/06/12/18 instead of 01/04/.../22.
The 9-year selection, the 21-day running mean, the chunked cat/timmean, and the
final multi-year timmean are all UNCHANGED.

Usage:
    python 02_t_climatology_6h.py BASE_START BASE_END [N_JOBS]
    e.g. python 02_t_climatology_6h.py 2026-05-14 2026-06-02 28

tbar files are written to $TBAR_DIR as t_YYYY-MM-DD_6hHH_tbar.nc
NOTE: choose BASE_START/END to extend +/-1 day beyond the S-file period
(stage 03 needs the +/-6 h climatology neighbours across day boundaries).
"""
import os
import sys
import glob
import pandas as pd
import xarray as xr
from cdo import Cdo
from joblib import Parallel, delayed

os.environ["PATH"] = os.environ.get("CDO_BIN_DIR", "") + ":" + os.environ["PATH"]

ERA5_ML_E5 = os.environ["ERA5_ML_E5"]
ERA5_ML_ET = os.environ["ERA5_ML_ET"]
data_path_ml = ERA5_ML_E5

SCRATCH_BASE = os.path.join(os.environ["SCRATCH_BASE"], "t_clim")
dump_path = os.path.join(SCRATCH_BASE, "tbar_dump")
temp_dir = os.path.join(dump_path, "tmp")
intermediate_path = os.path.join(SCRATCH_BASE, "bar_intermediate")
final_path = os.environ["TBAR_DIR"]
for p in [dump_path, intermediate_path, temp_dir, final_path]:
    os.makedirs(p, exist_ok=True)

HOURS_CSV = os.environ.get("HOURS_CSV", "0,6,12,18")
LEVELS = os.environ.get("LEVELS", "").strip()
REMAP = os.environ.get("REMAP_GRID", "r720x360")
LATBOX = os.environ.get("DOMAIN_LATBOX", "0,360,0,90").strip()   # NH cut from the start
SLOTS = ["00", "06", "12", "18"]          # was 01/04/07/10/13/16/19/22

era5_vars = {130: "t"}

# base period is set by main(); referenced by year-selection helpers
base_period_start = None
base_period_end = None


def process_era5_daily(year, month, day, var_num, var, ml="ml"):
    """Regrid one daily T file, 6-hourly snapshots only (no 3h mean)."""
    cdo = Cdo()
    date_str = f"{year}-{month:02d}-{day:02d}"

    if year == base_period_start.year:
        et_file = f"{ERA5_ML_ET}{var_num}/ET{ml}00_1H_{date_str}_{var_num}.grb"
        e5_file = f"{ERA5_ML_E5}{var_num}/E5{ml}00_1H_{date_str}_{var_num}.grb"
        var_file = et_file if os.path.exists(et_file) else e5_file
    else:
        var_file = f"{data_path_ml}{var_num}/E5{ml}00_1H_{date_str}_{var_num}.grb"

    if not os.path.exists(var_file):
        print(f"Skipped {var_file}: file does not exist", flush=True)
        return

    reg_file = os.path.join(dump_path, f"reg_{var}_{date_str}.nc")
    gp_file = os.path.join(dump_path, f"gp_{var}_{date_str}.nc")
    regrid_6h = os.path.join(dump_path, f"{var}_{date_str}_0.5deg_6h.nc")

    if os.path.exists(regrid_6h):
        print(f"{regrid_6h} already exists - skipping.", flush=True)
        return regrid_6h

    sel = f"-selhour,{HOURS_CSV}"
    if LEVELS:
        sel = f"-sellevel,{LEVELS} {sel}"
    cdo.setgridtype("regular", input=f"{sel} {var_file}",
                    output=reg_file, options="-f nc --eccodes")

    try:
        cdo.sp2gp(input=reg_file, output=gp_file, options="-P 16")
        input_for_remap = gp_file
    except Exception as e:
        print(f"sp2gp failed or not needed for {var}: {e}", flush=True)
        input_for_remap = reg_file

    # remap to 0.5 deg, then cut to the domain box (default NH) in one chained call
    if LATBOX:
        cdo.sellonlatbox(LATBOX, input=f"-remapbil,{REMAP} {input_for_remap}", output=regrid_6h)
    else:
        cdo.remapbil(REMAP, input=input_for_remap, output=regrid_6h)
    # (timselmean REMOVED)

    for f in (reg_file, gp_file):
        if os.path.exists(f):
            os.remove(f)
    return regrid_6h


def process_extended_period(start_date, end_date, var_num=130, var="t", n_jobs=8):
    extended_range = pd.date_range(start_date - pd.Timedelta(days=11),
                                   end_date + pd.Timedelta(days=11), freq="D")
    results = Parallel(n_jobs=n_jobs, verbose=10)(
        delayed(process_era5_daily)(d.year, d.month, d.day, var_num, var)
        for d in extended_range
    )
    return [r for r in results if r is not None]


def split_into_6hourly_files(regrid_6h_file, var, date_str):
    """Split the daily 6-hourly file (4 steps) into one file per slot."""
    ds = xr.open_dataset(regrid_6h_file)
    output_files = []
    for t_idx, slot in enumerate(SLOTS):
        single = ds.isel(time=t_idx)
        single = single.expand_dims(time=[single.time.values])
        out_file = os.path.join(dump_path, f"{var}_{date_str}_6h{slot}.nc")
        single.to_netcdf(out_file)
        single.close()
        output_files.append(out_file)
    ds.close()
    if os.path.exists(regrid_6h_file):
        os.remove(regrid_6h_file)
    return output_files


def compute_running_mean_split(center_date, var="t", window_days=21, chunk_size=3):
    """UNCHANGED 21-day running mean per slot (chunked cat/timmean)."""
    half_window = window_days // 2
    window = pd.date_range(center_date - pd.Timedelta(days=half_window),
                           center_date + pd.Timedelta(days=half_window), freq="D")
    outputs = []
    for slot in SLOTS:
        files = [os.path.join(dump_path, f"{var}_{d.strftime('%Y-%m-%d')}_6h{slot}.nc")
                 for d in window
                 if os.path.exists(os.path.join(dump_path, f"{var}_{d.strftime('%Y-%m-%d')}_6h{slot}.nc"))]
        if len(files) < window_days:
            continue
        out_file = os.path.join(intermediate_path,
                                f"{var}_{center_date.strftime('%Y-%m-%d')}_6h{slot}_mean.nc")
        if os.path.exists(out_file):
            continue
        tmp_chunks = []
        cdo = Cdo()
        for i in range(0, len(files), chunk_size):
            sub = files[i:i + chunk_size]
            tmp = os.path.join(temp_dir, f"chunk_{i}_{center_date.strftime('%Y-%m-%d')}_{slot}.nc")
            cdo.cat(input=" ".join(sub), output=tmp)
            tmp_m = tmp.replace(".nc", "_mean.nc")
            cdo.timmean(input=tmp, output=tmp_m)
            os.remove(tmp)
            tmp_chunks.append(tmp_m)
        tmp_cat = os.path.join(temp_dir, f"tmp_{slot}_{center_date.strftime('%Y-%m-%d')}.nc")
        cdo.cat(input=" ".join(tmp_chunks), output=tmp_cat)
        cdo.timmean(input=tmp_cat, output=out_file)
        os.remove(tmp_cat)
        for f in tmp_chunks:
            os.remove(f)
        outputs.append(out_file)
    return outputs


def era5_file_exists_for_date(d, var_num=130, ml="ml"):
    date_str = d.strftime("%Y-%m-%d")
    if d.year == base_period_start.year:
        et_file = f"{ERA5_ML_ET}{var_num}/ET{ml}00_1H_{date_str}_{var_num}.grb"
        e5_file = f"{ERA5_ML_E5}{var_num}/E5{ml}00_1H_{date_str}_{var_num}.grb"
        return os.path.exists(et_file) or os.path.exists(e5_file)
    e5_file = f"{ERA5_ML_E5}{var_num}/E5{ml}00_1H_{date_str}_{var_num}.grb"
    return os.path.exists(e5_file)


def paper_climatology_years(year, n_years=9, first_year=1940, last_year=2025):
    """The paper's climatology (Hotz et al. 2024): a centred n-year window of
    21-day running means, clamped to first_year..last_year, so it never
    extends past 2025 (2021-2026 all get 2017-2025). Identical to the
    analogue pipeline's climatology_years_for_analogue for 1940-2025."""
    start, end = year - n_years // 2, year + n_years // 2
    if end > last_year:
        start, end = last_year - n_years + 1, last_year
    elif start < first_year:
        start, end = first_year, first_year + n_years - 1
    return list(range(start, end + 1))


def choose_climatology_years(base_start, base_end, var_num=130, ml="ml"):
    """CLIM_YEARS (config.sh) if set, otherwise the paper rule for the event
    year. Never any year after 2025: the old automatic rule kept 2026 in its own
    climatology (2018-2026) once ERA5T for 2026 was complete."""
    spec = os.environ.get("CLIM_YEARS", "").replace(",", " ").strip()
    if spec:
        if "-" in spec:
            y0, y1 = (int(x) for x in spec.split("-"))
            years = list(range(y0, y1 + 1))
        else:
            years = [int(y) for y in spec.split()]
        how = "CLIM_YEARS"
    else:
        years = paper_climatology_years(base_start.year)
        how = "paper rule (centred 9 years, clamped to 1940-2025)"
    if max(years) > 2025:
        raise SystemExit(f"climatology years {years} extend past 2025 - not the paper definition")
    print(f"Climatology years ({how}): {years}", flush=True)
    return years


def compute_all_years_split(base_start, base_end, var="t", var_num=130, n_jobs=8):
    base_years = choose_climatology_years(base_start, base_end, var_num=var_num)
    for year in base_years:
        print(f"\n=== Processing year {year} ===", flush=True)
        start = base_start.replace(year=year)
        end = base_end.replace(year=year)
        date_range = pd.date_range(start, end, freq="D")

        existing = [os.path.join(intermediate_path, f"{var}_{d.strftime('%Y-%m-%d')}_6h{s}_mean.nc")
                    for d in date_range for s in SLOTS]
        existing = [f for f in existing if os.path.exists(f)]
        if len(existing) == len(date_range) * len(SLOTS):
            print(f"All running-mean files for {year} already exist. Skipping.", flush=True)
            continue

        for f in glob.glob(os.path.join(dump_path, f"{var}_{year}-*_0.5deg_6h.nc")):
            os.remove(f)

        daily_files = process_extended_period(start, end, var_num, var, n_jobs=n_jobs)
        print(f"Splitting {len(daily_files)} daily files...", flush=True)
        Parallel(n_jobs=n_jobs, verbose=10)(
            delayed(split_into_6hourly_files)(
                df, var,
                date_str=os.path.basename(df).split(f"{var}_")[1].split("_0.5deg_6h.nc")[0]
            ) for df in daily_files
        )
        print(f"Computing 21-day running means for {len(date_range)} days...", flush=True)
        Parallel(n_jobs=n_jobs, verbose=10)(
            delayed(compute_running_mean_split)(d, var) for d in date_range
        )
        for f in glob.glob(os.path.join(dump_path, f"{var}_*.nc")):
            os.remove(f)
        print(f"Cleaned hourly files for {year}.", flush=True)
    print("\nAll years processed.", flush=True)
    return base_years


def process_finalclim(d, var, climatology_years):
    cdo = Cdo()
    for slot in SLOTS:
        hourly_files = []
        for year in climatology_years:
            try:
                file_date = d.replace(year=year)
            except ValueError:
                continue
            fp = os.path.join(intermediate_path,
                              f"{var}_{file_date.strftime('%Y-%m-%d')}_6h{slot}_mean.nc")
            if os.path.exists(fp):
                hourly_files.append(fp)
        if not hourly_files:
            print(f"No files for {d.date()} slot {slot}, skipping", flush=True)
            continue
        out_file = os.path.join(final_path, f"{var}_{d.strftime('%Y-%m-%d')}_6h{slot}_tbar.nc")
        if os.path.exists(out_file):
            continue
        tmp_merged = os.path.join(final_path, f"tmp_merged_{d.strftime('%Y-%m-%d')}_6h{slot}.nc")
        cdo.cat(input=" ".join(hourly_files), output=tmp_merged)
        cdo.timmean(input=tmp_merged, output=out_file)
        if os.path.exists(tmp_merged):
            os.remove(tmp_merged)
        print(f"Climatology written: {out_file}", flush=True)


def compute_final_climatology_parallel(base_start, base_end, climatology_years, var="t", n_jobs=8):
    date_range = pd.date_range(base_start, base_end, freq="D")
    Parallel(n_jobs=n_jobs, verbose=10)(
        delayed(process_finalclim)(d, var, climatology_years) for d in date_range
    )


def main():
    global base_period_start, base_period_end
    if len(sys.argv) < 3:
        sys.exit("usage: 02_t_climatology_6h.py BASE_START BASE_END [N_JOBS]")
    base_period_start = pd.Timestamp(sys.argv[1])
    base_period_end = pd.Timestamp(sys.argv[2])
    n_jobs = int(sys.argv[3]) if len(sys.argv) > 3 else 28

    years = compute_all_years_split(base_period_start, base_period_end,
                                    var="t", var_num=130, n_jobs=n_jobs)
    compute_final_climatology_parallel(base_period_start, base_period_end,
                                       years, var="t", n_jobs=n_jobs)


if __name__ == "__main__":
    main()
