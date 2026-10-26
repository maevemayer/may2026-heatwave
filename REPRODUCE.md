# REPRODUCE — run order for Mayer et al. (2026)

1. **Environment**
   ```
   conda env create -f environment.yml && conda activate may2026
   pip install colindex2            # cut-off-low detector
   ```
2. **Paths** — edit `config.yaml`. Most keys already point at `./06_derived_data` (shipped);
   set the external ones only if you re-derive from scratch (see `DATA_MANIFEST.md`).
   Figure code refers to data via `${KEY}` placeholders — substitute your `config.yaml`
   values (the loader snippet in `DATA_MANIFEST.md`), or set the 2-3 paths each file uses.
3. **Download the raw inputs that are not shipped** (optional — only Fig 1 & 2 need raw; a
   free CDS account + `~/.cdsapirc` required):
   - Fig 1 forecast: `01_data_acquisition/download_ifs_z500_t2m.ipynb`
   - Fig 1 ERA5 daily: already shipped; regenerate with `download_era5_daily_t2m_z500.py` if needed
   - Fig 2 PV/θ/MSLP: `05_figures/fig2_pv_rwb_dynamics/download_preprocess_pv_theta_mslp.ipynb`
4. **Make the figures** (each is standalone):

   | figure | run |
   |---|---|
   | 1 overview + predictability | `05_figures/fig1_overview_and_predictability/fig1_overview_ifs_composites.ipynb` |
   | 2 PV / Rossby-wave-breaking | `05_figures/fig2_pv_rwb_dynamics/fig2_pv_rossby_wave_breaking.ipynb` (needs `find_rwb.py`) |
   | 3 cut-off in IFS & AIFS | `fig3_cutoff_ifs_aifs/fig3_ifs_cutoff_ridgecontrol.py` + `fig3_aifs_cutoff_matched.ipynb` |
   | 4 historical analogues | `04_analogue_selection/fig4_analogues_and_ridge_control.ipynb` |
   | 5 vertical decomposition | `05_figures/fig5_vertical_decomposition/fig5_vertical_lagrangian_decomposition.py` |
   | 6 COL-slope regression | `05_figures/fig6_slope_regression/fig6_col_vs_decomposition_slope.py` |
   | 7 genesis / origin maps | `05_figures/fig7_genesis_origin_maps/fig7_genesis_advective_origin.py` |
   | 8 Iberian COL trend | `02_col_catalogue/era5_col_catalogue_and_trend.ipynb` |
   | S1–S11 | `05_figures/supplements/*` (+ S2 `02_col_catalogue/col_200_500hPa_threshold_matching.ipynb`, S3/S4/S5 in the Fig 1 / AIFS notebooks, S11 in the Fig 4 notebook) |

   Ridge-control coefficients (Fig 3/4) are shipped; regenerate with
   `04_analogue_selection/learn_ridge_control_z500_t2m_*.py` if desired.

Notes: the LAGRANTO trajectory outputs (Fig 5/6/7, S6–S10) are shipped in
`06_derived_data/` (`ifs_tanom_decomp/`, `analogue_origin_6h_windfix/`, `origin_*`), so you do
not need to re-run `03_lagranto_tracing/` unless starting from raw ERA5.
