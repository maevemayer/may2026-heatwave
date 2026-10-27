# Mayer et al. (2026) - The May 2026 Western European Heatwave: code & data

Reproducibility package for the manuscript and supplement. Organised by workflow
section; every figure maps to one script/notebook. Personal absolute paths have been
replaced by `${KEYS}` set in **`config.yaml`**. Large inputs are not shipped - see
**`DATA_MANIFEST.md`**; small derived inputs are in **`06_derived_data/`**.

## Requirements: LAGRANTO

The trajectory tracing in `03_lagranto_tracing/` requires the **shell-script
version of LAGRANTO.ECMWF**, which is not redistributed here. It can be requested
from Michael Sprenger (https://iacweb.ethz.ch/staff/sprenger/lagranto/). If you
use it, please cite:

- Sprenger, M. and Wernli, H. (2015): The LAGRANTO Lagrangian analysis tool -
  version 2.0. *Geosci. Model Dev.*, 8, 2569-2586. doi:10.5194/gmd-8-2569-2015.
- Wernli, H. and Davies, H. C. (1997): A Lagrangian-based analysis of
  extratropical cyclones. I: The method and some applications. *Q. J. R.
  Meteorol. Soc.*, 123, 467-489. doi:10.1002/qj.49712353811.

## Layout
```
01_data_acquisition/   download the ECMWF/TIGGE forecast fields from CDS
02_col_catalogue/      ERA5 cut-off-low detection catalogue (colindex2) + trends   -> Fig 8, S2
03_lagranto_tracing/   LAGRANTO trajectory tracing + Lagrangian T-decomposition    -> feeds Fig 5/6/7, S6-S10
04_analogue_selection/ flow-analogue selection + learned Z500->T2m ridge control   -> Fig 3(ridge), Fig 4, S11
05_figures/            one folder per main figure (+ supplements/)
06_derived_data/       small derived files used directly by the figure code
config.yaml            YOUR data paths
environment.yml        conda environment
DATA_MANIFEST.md       how to obtain the large inputs, per figure
```

## Figure -> code -> main data
| Fig | code | key input(s) |
|-----|------|--------------|
| 1 Overview + IFS composites        | 05_figures/fig1_overview_and_predictability/fig1_overview_ifs_composites.ipynb | ERA5 daily t2m/z500 (`${ERA5_DERIVED}`); IFS forecast GRIBs (`${IFS_FORECAST_DIR}`) |
| 2 PV / Rossby-wave-breaking        | 05_figures/fig2_pv_rwb_dynamics/fig2_pv_rossby_wave_breaking.ipynb | ERA5 PV / theta (model-level) |
| 3 Cut-off effect in IFS & AIFS     | 05_figures/fig3_cutoff_ifs_aifs/fig3_ifs_cutoff_ridgecontrol.py  +  fig3_aifs_cutoff_matched.ipynb | 06_derived_data (ifs_box csv, ridge JSON); AIFS regional metrics + COL d4 (`${DERIVED_DATA}/aifs`) |
| 4 Historical analogues             | 04_analogue_selection/fig4_analogues_and_ridge_control.ipynb | ERA5 daily z500/t2m (`${ANALOGUE_BASE}`); COL catalogue (06_derived_data) |
| 5 Vertical Lagrangian decomposition| 05_figures/fig5_vertical_decomposition/fig5_vertical_lagrangian_decomposition.py | IFS tanom_decomp (`${IFS_TANOM_DECOMP}`); analogue origin (`${ANALOGUE_BASE}`) |
| 6 COL-activity slope regression    | 05_figures/fig6_slope_regression/fig6_col_vs_decomposition_slope.py | IFS tanom_decomp; analogue origin; COL catalogue |
| 7 Genesis / advective origin maps  | 05_figures/fig7_genesis_origin_maps/fig7_genesis_advective_origin.py | IFS tanom_decomp; analogue origin_contributions.nc |
| 8 Iberian COL trend 1940-2025      | 02_col_catalogue/era5_col_catalogue_and_trend.ipynb | COL catalogue (06_derived_data) + ERA5 z500 1deg |
| S1 study regions map               | 05_figures/supplements/figS1_study_regions_map.ipynb | ERA5 daily t2m/z500 |
| S2 200/500 hPa matching            | 02_col_catalogue/col_200_500hPa_threshold_matching.ipynb | COL z500 + z200 catalogues (06_derived_data) |
| S6 IFS composite bias              | 05_figures/supplements/figS6_ifs_composite_bias.py | IFS tanom_decomp + ERA5 origin |
| S7 land-only decomposition         | 05_figures/supplements/figS7_vertical_decomposition_land.py | as Fig 5 |
| S8 origin maps IFS/analogue        | 05_figures/fig7_genesis_origin_maps/figS8_origin_maps_ifs_analogue.py | IFS tanom_decomp + analogue origin |
| S9 analogue origin, all processes  | 05_figures/supplements/figS9_analogue_origin_allprocess.py | analogue origin_contributions.nc |
| S10 ERA5 origin, all processes     | 05_figures/supplements/figS10_era5_origin_allprocess.py | ERA5 origin_contributions.nc |
| S11 2022 / June 2026 analogues     | cell in 04_analogue_selection/fig4_analogues_and_ridge_control.ipynb | analogue data |

(S3/S4/S5 predictability & AIFS composites are produced by the Fig 1 notebook and the AIFS
notebooks in fig3_cutoff_ifs_aifs / 02_col_catalogue.)

## Running
1. `conda env create -f environment.yml && conda activate may2026` (then `pip install colindex2`).
2. Edit `config.yaml` with your paths; fetch the large inputs per `DATA_MANIFEST.md`.
3. Scripts/notebooks reference data via the `${KEYS}` - substitute your `config.yaml`
   values (a tiny loader is suggested in DATA_MANIFEST.md; search each file for `${` to see
   which keys it needs). Figure scripts are standalone (`python 05_figures/.../fig*.py`).

## Notes
- Cut-off detection uses **colindex2** (`pip install colindex2`).
- Trajectory tracing uses **LAGRANTO** (see 03_lagranto_tracing/README.md); the heavy
  trajectory outputs are not shipped - the decomposition products needed by the figures are
  the `tanom_*` files (`${IFS_TANOM_DECOMP}`) and `origin_contributions.nc` (`${ANALOGUE_BASE}`).
- Notebook outputs were cleared for size/privacy; re-run to regenerate.
