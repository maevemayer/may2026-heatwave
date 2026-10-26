#!/usr/bin/env python3
"""Download the ERA5 daily fields used for Figure 1 (and the climatology/detrending throughout),
from the Copernicus Climate Data Store (needs a free CDS account + `pip install cdsapi` and a
~/.cdsapirc key: https://cds.climate.copernicus.eu/how-to-api).

Produces, after the cdo daily-mean step below:
  era5_t2m_daily_NA_<Y0>_<Y1>_FULL_0p5deg.nc   (2 m temperature, 0.5 deg)
  era5_z500_daily_NA_<Y0>_<Y1>_FULL_1deg.nc    (500 hPa geopotential height, 1 deg)
matching ${ERA5_DERIVED}/${ANALOGUE_BASE}. For Figure 1 alone you only need the event month
and the 1995-2025 May climatology window - set Y0=1995 to keep the download small; set Y0=1940
to reproduce the climatology/trend analyses (large, ~decades of hourly data).
"""
import cdsapi, calendar, os, subprocess
c = cdsapi.Client()
OUTDIR = "."                       # set your output folder
Y0, Y1 = 1995, 2026                # 1940..2026 for the full series used elsewhere in the paper
AREA   = [80, -130, 10, 80]        # N, W, S, E  (North-Atlantic / Europe crop)
HOURS  = [f"{h:02d}:00" for h in range(24)]
os.makedirs(OUTDIR, exist_ok=True)

def grab(dataset, req, fname):
    fp = os.path.join(OUTDIR, fname)
    if os.path.exists(fp) and os.path.getsize(fp) > 0:
        print("skip", fname); return fp
    c.retrieve(dataset, req, fp); print("done", fname); return fp

for year in range(Y0, Y1 + 1):
    # --- 2 m temperature, 0.5 deg, single levels ---
    grab("reanalysis-era5-single-levels",
         {"product_type": "reanalysis", "variable": "2m_temperature",
          "year": str(year), "month": [f"{m:02d}" for m in range(1,13)],
          "day": [f"{d:02d}" for d in range(1,32)], "time": HOURS,
          "area": AREA, "grid": [0.5, 0.5], "data_format": "netcdf"},
         f"era5_t2m_hourly_NA_{year}_0p5deg.nc")
    # --- 500 hPa geopotential, 1 deg, pressure levels ---
    grab("reanalysis-era5-pressure-levels",
         {"product_type": "reanalysis", "variable": "geopotential", "pressure_level": "500",
          "year": str(year), "month": [f"{m:02d}" for m in range(1,13)],
          "day": [f"{d:02d}" for d in range(1,32)], "time": HOURS,
          "area": AREA, "grid": [1.0, 1.0], "data_format": "netcdf"},
         f"era5_z500_hourly_NA_{year}_1deg.nc")

# --- daily means + concatenation (needs cdo on PATH) ---
# per-year hourly -> daily mean, then mergetime into the FULL files:
#   for f in era5_t2m_hourly_NA_*_0p5deg.nc; do cdo daymean "$f" "${f/hourly/daymean}"; done
#   cdo mergetime era5_t2m_daymean_NA_*_0p5deg.nc era5_t2m_daily_NA_%d_%d_FULL_0p5deg.nc
# (geopotential is m^2 s^-2; divide by g0=9.80665 for height in m, as the figure code does.)
print("downloaded hourly ERA5; run the cdo daymean+mergetime step in the comment above.")
