# 07_aifs_simulation — AIFS-ENS ensemble boosting (initial-condition ensemble)

Scripts that generate the AIFS-ENS forecast ensemble used in the AIFS analysis
(Fig 3 and the AIFS composites). Contributed by **Paolo Chevalier**; adapted
here from the repository
`rare-event-simulation-using-ai-weather-forecasting` (folders `DL_INPUTS_ENS`
and `SIMPLE_BOOSTING`). Personal cluster paths have been replaced with
`/path/to/...` placeholders; the analysis logic is unchanged.

## Pipeline

```
DL_INPUTS_ENS/              1. download + preprocess the initial conditions
  aifs_inputs_dl+preprocess.py   ERA5 at (t-6h, t) from the CDS, regridded 0.25deg -> N320,
                                 saved as the AIFS-ENS v1 init-state .pkl
  example_use.sh                 SLURM wrapper, one init date per block (15-25 May 2026)

SIMPLE_BOOSTING/            2. run the ensemble
  run_aifs_crps_simple_boosting.py   drives AIFS-ENS CRPS (anemoi SimpleRunner, GPU);
                                     N_MEMBERS perturbed members + control, saves
                                     selected variables/levels every 6 h as per-member .pkl
  example_use.sh                     SLURM (A100) wrapper

  POSTPROCESS/             3. assemble to NetCDF
    simple_boosting_postprocess.py   stacks the per-member .pkl into one NetCDF per
                                     variable, dims (number, step, values) on the N320 grid
    example_use.sh                   SLURM wrapper
```

## Requirements

- **AIFS-ENS CRPS v1.0 model weights** (`aifs-ens-crps-1.0.ckpt`) — obtained from
  ECMWF; see the data-availability statement in the manuscript. The weights are
  *not* redistributed here.
- Python: `anemoi-inference` (`SimpleRunner`), `earthkit-data` / `earthkit-regrid`,
  `ecmwf-opendata`, `cdsapi`, `xarray`, `numpy`, `pandas`. A CDS account is needed
  for the initial conditions (`~/.cdsapirc`).
- An N320 lat/lon lookup pickle for the postprocessing step (`N320_COORDS`
  env var, default `n320_coordinates.pkl`).

## Data availability note

Only the **model weights** are needed to reproduce the method; the initial
conditions are simply ERA5 at `(t-6h, t)` (freely available from the CDS, see
`../DATA_MANIFEST.md`). The raw AIFS-ENS simulations (~1 TB) are therefore **not
archived**. The ensemble was run without a fixed seed or GPU determinism, so a
re-run reproduces the results *in distribution* (statistically equivalent), not
bit-for-bit.
