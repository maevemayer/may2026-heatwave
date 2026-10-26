#!/bin/bash
# ===========================================================================
# Shared configuration for the corrected 6-hourly LAGRANTO pipeline.
# Source this from every stage:  source config.sh
#
# CORRECTION vs the old notebooks:
#   * NO 3-hourly means (cdo timselmean removed everywhere).
#   * The 6-hourly INSTANTANEOUS snapshots (00/06/12/18 UTC) are selected
#     directly from the hourly ERA5 files, BEFORE the spectral transform
#     (-selhour applied ahead of sp2gp/remapbil -> big speed-up).
#   * Because snapshots are instantaneous (not centred 3h-bin means), the
#     old +1 h timestamp offset is gone -> the 04_02 rename step is dropped.
#   * The DTDT climatology time-stencil is 6 h instead of 3 h (forced by the
#     new cadence). All other math is unchanged from the notebooks.
# ===========================================================================

# ---- CDO ----

# ===== USER CONFIG (Zenodo deposit): set these to YOUR paths; see ../config.yaml =====
: "${LAGRANTO_BASE:=/EDIT/lagranto_base}"
: "${LAGRANTO_RERUN:=/EDIT/lagranto_rerun}"
: "${WORK:=/EDIT/work}"
: "${SCRATCH:=/EDIT/scratch}"
: "${ANALOGUE_BASE:=/EDIT/analogue_base}"
: "${IFS_TANOM_DECOMP:=/EDIT/ifs_tanom_decomp}"
: "${MSE_BASE:=/EDIT/mse_base}"
: "${DERIVED_DATA:=/EDIT/derived_data}"
: "${HOME:=/EDIT/home}"
# =====================================================================================
export CDO_BIN_DIR="/path/to/cdo/bin"
export PATH="${CDO_BIN_DIR}:${PATH}"

# ---- Slurm ----
export SLURM_ACCOUNT="uc1275"
export SLURM_PARTITION="shared"

# ---- ERA5 input pools (fixed DKRZ locations) ----
export ERA5_ML_ET="/pool/data/ERA5/ET/ml/an/1H/"   # model-level, current-year (ET)
export ERA5_ML_E5="/pool/data/ERA5/E5/ml/an/1H/"   # model-level, reanalysis (E5)
export ERA5_SF_ET="/pool/data/ERA5/ET/sf/an/1H/"   # surface, current-year (ET)
export ERA5_SF_E5="/pool/data/ERA5/E5/sf/an/1H/"   # surface, reanalysis (E5)

# ---- Output tree (NEW, so nothing overwrites the old 3h data) ----
# _tc10hpa: NEW tree for the along-trajectory dT_clim/dp method, so the
# results of the old DTDP-field method in .../corrected_6h stay untouched
# _clim2017_2025: the 2026 event with a FIXED 2017-2025 climatology (the old run
# kept 2026 in its own climatology, 2018-2026)
# _windfix: new P-files with the cos(lat) wind fix (stage 01); the S-files of the
# _clim2017_2025 run are correct (no U/V in them) and are linked in, not rebuilt
export BASE_OUT="${LAGRANTO_BASE}/RERUN/corrected_6h_tc10hpa_clim2017_2025_windfix"
export SCRATCH_BASE="${SCRATCH}/lagranto/2026_MAY_6h_tc10hpa_clim2017_2025_windfix"   # fresh: no stale wind files
export DATA_DIR="${BASE_OUT}/data"        # P- and S-files land here (stages 01-03)
# P/S files that stage 04 READS (here: the new P-files + the linked S-files)
export PS_DATA_DIR="${PS_DATA_DIR:-${DATA_DIR}}"
# Stages submit_all.sh submits. Without 03, the S-files of S_REUSE_DIR are linked
# into DATA_DIR first. Full run: STAGES="01 02 03 04 05 06".
export STAGES="${STAGES:-01 04 05 06}"
export S_REUSE_DIR="${S_REUSE_DIR:-${LAGRANTO_BASE}/RERUN/corrected_6h_tc10hpa_clim2017_2025/data}"

# Climatology years (stage 02). Empty = the paper rule: centred 9-year window,
# clamped to 1940-2025 (2026 -> 2017-2025, as for 2021-2025; Hotz et al. 2024).
# An explicit list/range overrides it; years after 2025 are refused.
export CLIM_YEARS="${CLIM_YEARS:-}"
export TBAR_DIR="${BASE_OUT}/tbar"        # final T climatology (tbar) files
export RESULTS_DIR="${BASE_OUT}/results_6h"
mkdir -p "${DATA_DIR}" "${TBAR_DIR}" "${RESULTS_DIR}" "${SCRATCH_BASE}"

# ---- The 6-hourly cadence ----
export HOURS_CSV="0,6,12,18"              # cdo -selhour argument
export SLOTS="00 06 12 18"               # human-readable slot labels
export STEPS_PER_DAY=4

# ---- Optional vertical restriction (see README) ----
#   empty  -> keep ALL 137 model levels (math identical to notebooks) [DEFAULT]
#   "70/137" -> lowest 68 levels (~250 hPa -> surface); WILL clip 200 hPa
#               trajectory starts and any air above ~250 hPa. Enable only
#               after confirming your trajectories stay below that top.
export LEVELS=""

# ---- ERA5 remap target grid ----
export REMAP_GRID="r720x360"             # 0.5 deg global (unchanged)

# ---- Domain cut (applied right after remap, in stages 01 & 02) ----
#   ERA5 is cut to this lon,lat box the moment it is regridded, so every
#   downstream file (p, dp, P-files, tbar) is domain-only -> less I/O and RAM.
#   Format: lon1,lon2,lat1,lat2 (CDO -sellonlatbox).
#     "0,360,0,90"  -> whole Northern Hemisphere [DEFAULT]
#     "-15,90"      caveat: 15-day back-trajectories from W Europe stay well
#                   inside the NH; any parcel that would leave the box is flagged
#                   -999 by LAGRANTO and masked in stage 06. Widen the southern
#                   edge (e.g. "0,360,-15,90") if you want an equator buffer.
#     empty         -> keep the global grid (no cut).
export DOMAIN_LATBOX="0,360,0,90"

# ---- LAGRANTO ----
export LAGRANTO_ROOT="${WORK}/01_LAGRANTO/lagranto.era5.test/lagranto.shellscript"
# tracevars with TCLIM traced +/-10 hPa (vertical gradient along the trajectory)
export TRACEVARS_FILE="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}/tracevars_tc10hpa.shellscript"
# IMPORTANT: set  set interval = "6"  in
#   ${LAGRANTO_ROOT}/goodies/datelist.sh
# before running stage 05 (was "3" for the old 3-hourly data).
