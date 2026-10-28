# DATA MANIFEST - obtaining the large inputs

Everything the code needs that is **too large to ship** is listed here with where to get
it. Small derived files are already in `06_derived_data/`.

The full derived-data archive (~11 GB) is deposited on Zenodo:
**DOI [10.5281/zenodo.23666704](https://doi.org/10.5281/zenodo.23666704)**.
Download and unzip it into `06_derived_data/` to reproduce the figures without
rebuilding the heavy intermediate products.

## Shipped (06_derived_data/)
- `col_features_AMJJA_1940_2025_largebox.nc`      ERA5 500 hPa cut-off/trough catalogue 1940-2025 (Fig 8, S2)
- `col_features_AMJJA_1940_2025_largebox_z200.nc` ERA5 200 hPa cut-off catalogue (S2 matching)
- `ridge_control_t2m_z500_{AMJJA,MJJA,JJA}.json`  learned Z500->T2m ridge-control coefficients (Fig 3/4)
- `target_region_daily_means_1940_2026.nc`        daily target-box t2m/z500 + RMST
- `ifs_box_daily_t2m_z500_20260525.csv`           IFS member box-mean t2m/z500 on 25 May 2026 (Fig 3)
- `composite_all_members_t2m_z500_20260525.csv`   all IFS members for the Fig 1 composites
- `figure5b_correlations_ridge.csv`, `figure5b_robin_attribution.csv`  Fig 4 numbers

## ERA5 (Copernicus CDS / your archive)  ->  ${ERA5_DERIVED} and ${ANALOGUE_BASE}/rawdat*
Daily fields, North-Atlantic/European crop, 1940-2026:
- `era5_t2m_daily_NA_1940_2026_FULL_0p5deg.nc`  (2 m temperature, 0.5 deg)  **use the 0.5 deg file only**
- `era5_z500_daily_NA_1940_2026_FULL_1deg.nc`   (500 hPa geopotential, 1 deg)
ERA5 is freely available from the Copernicus Climate Data Store (reanalysis-era5-single-levels /
-pressure-levels). Fig 2 additionally needs model-level PV/theta.

## ECMWF/TIGGE forecast GRIBs  ->  ${IFS_FORECAST_DIR}
Download with the notebooks in `01_data_acquisition/` (needs a CDS account + `cdsapi`):
- `download_ifs_z500_t2m.ipynb`   -> gh500 (1 deg & 0.5 deg) and 2t (0.5 deg), 6-hourly, May 2026 inits
Dataset: `tigge-forecasts`, origin ecmwf, perturbed_forecast. AIFS forecasts: ECMWF AIFS-ENS.

## IFS Lagrangian temperature decomposition  ->  ${IFS_TANOM_DECOMP}
Per-member trajectory files `tanom_YYYYMMDD_HH` with the decomposition columns
(tanom0, tanom_advective, tanom_adiab1, tanom_adiab2, tanom_diabatic, genesis coords).
Produced by the LAGRANTO workflow in `03_lagranto_tracing/` from ERA5; feeds Fig 5/6/7 and S6-S10.

## Analogue trajectory origin products  ->  ${ANALOGUE_BASE}
`rank*/origin_contributions.nc` and pooled origin files, produced by `03_lagranto_tracing`
over the selected analogue dates. Used by Fig 5/6/7, S8-S11.

## AIFS derived (${DERIVED_DATA}/aifs)
`regional_metrics_Init_*.csv`, `COL_d4_index_Init_*_UP_EX1.csv`, `cols_index_Init_*.csv`,
`pdf_data_with_threshold.csv`, `significance_cutoff_vs_no_cutoff.nc`,
`vertical_temperature_profiles_Init_*.nc` - produced by the AIFS notebooks
(fig3_aifs_cutoff_matched.ipynb, aifs_col_detection.ipynb). Z500-matched pool: daily-mean,
cos(lat)-weighted box mean (5.5W-5E, 44.5-55N) > 5860 gpm (g0=9.80665).

## Suggested tiny path loader
Each file refers to data as e.g. `${ERA5_DERIVED}/era5_z500...nc`. To resolve at runtime:
```python
import yaml, os, re
CFG = yaml.safe_load(open("config.yaml"))
def P(s): return re.sub(r"\$\{(\w+)\}", lambda m: str(CFG[m.group(1)]), s)
# e.g.  xr.open_dataarray(P("${ERA5_DERIVED}/era5_z500_daily_NA_1940_2026_FULL_1deg.nc"))
```
or simply search each script for `${` and hard-set the two or three paths it uses.
