#!/bin/bash
# ===========================================================================
# Shared configuration for the ANALOGUE 6-hourly LAGRANTO pipeline.
# Source this from every submit script:  source config_analogues.sh
#
# This mirrors ../config.sh (the 2026-event pipeline). The analogue stages use
# the SAME setup as the event (6-hourly snapshots, NH cut, r720x360 remap,
# 21-day running mean, +/-6 h DTDT stencil, identical LAGRANTO start box /
# @profile(1000,200,17) / 15-day back trajectories, identical stage-05/06 math).
# The ONE intentional difference is the climatology window:
#   * event    -> trailing 9 years (2026 is incomplete, so the previous 9),
#   * analogue -> CENTERED 9-year window (analogue_year +/-4, clamped 1940-2025).
# Each analogue is referenced to its own era's climate; that is by design.
#
# Everything below is read by the Python stages (os.environ) and the Slurm
# wrappers, so paths live in ONE place. Override any line before sourcing, or
# export it on the sbatch command line.
# ===========================================================================

# ---- Slurm ----
export ANALOGUE_ACCOUNT="${ANALOGUE_ACCOUNT:-uc1275}"
export ANALOGUE_PARTITION="${ANALOGUE_PARTITION:-shared}"

# ---- Storage layout (all NEW paths for this rerun, so nothing collides) ----
#   ANALOGUE_HOME        : scripts + Slurm logs + the ORIGIN output (uc1275).
#   ANALOGUE_ORIGIN_DIR  : where stage 06 writes origin_contributions.nc
#                          (rankNN_MJ_<date>/origin_contributions.nc) — the "new data dir".
#   ANALOGUE_FINAL_BASE  : P/S files, trajectories and the t_clim/ climatology,
#                          on a NEW bu1561 path (kept off /work, large & persistent).
#   ANALOGUE_SCRATCH_BASE: NEW scratch root (Python appends the per-stage subdir).
export ANALOGUE_HOME="${ANALOGUE_HOME:-${LAGRANTO_BASE}/RERUN/corrected_lagranto_pipeline_6h/analogues}"
# _tc10hpa: NEW output trees for the along-trajectory dT_clim/dp method, so the
# results of the old DTDP-field method stay untouched
# _windfix: new P-files with the cos(lat) wind fix (stage 01) + the old, correct
# S-files linked in (no U/V in them), then stages 04-06; all output in new trees
export ANALOGUE_ORIGIN_DIR="${ANALOGUE_ORIGIN_DIR:-${ANALOGUE_HOME}_windfix}"
export ANALOGUE_LOG_DIR="${ANALOGUE_LOG_DIR:-${ANALOGUE_ORIGIN_DIR}/slurm_logs}"
export ANALOGUE_FINAL_BASE="${ANALOGUE_FINAL_BASE:-${LAGRANTO_RERUN}/analogues_windfix}"
# P/S files that stage 04 READS: the new tree (new P-files + linked S-files).
export ANALOGUE_PS_BASE="${ANALOGUE_PS_BASE:-${ANALOGUE_FINAL_BASE}}"
# Stages submit_all_analogues.sh submits. Without 03, each rank's S-files are linked
# from ANALOGUE_S_REUSE_BASE (found by date: rank*_MJ_<date>). Full run: "01 02 03 04 05 06".
export STAGES="${STAGES:-01 04 05 06}"
export ANALOGUE_S_REUSE_BASE="${ANALOGUE_S_REUSE_BASE:-${LAGRANTO_RERUN}/analogues}"
export ANALOGUE_SCRATCH_BASE="${ANALOGUE_SCRATCH_BASE:-${SCRATCH}/lagranto/RERUN_2026_MAY_6h_windfix}"   # fresh: no stale wind files

# ---- Python kernel (the analogues-2026 env). EDIT if your env moved. ----
export ANALOGUE_PYTHON="${ANALOGUE_PYTHON:-${ANALOGUE_BASE}/kernel/analogues-2026/bin/python}"

# ---- CDO ----
export CDO_BIN_DIR="${CDO_BIN_DIR:-/path/to/cdo/bin}"
export PATH="${CDO_BIN_DIR}:${PATH}"

# ---- ERA5 input pools (fixed DKRZ locations; same as the event) ----
export ERA5_ML_ET="${ERA5_ML_ET:-/pool/data/ERA5/ET/ml/an/1H}"
export ERA5_ML_E5="${ERA5_ML_E5:-/pool/data/ERA5/E5/ml/an/1H}"
export ERA5_SF_ET="${ERA5_SF_ET:-/pool/data/ERA5/ET/sf/an/1H}"
export ERA5_SF_E5="${ERA5_SF_E5:-/pool/data/ERA5/E5/sf/an/1H}"

# ---- LAGRANTO install (same as the event) ----
export LAGRANTO_ROOT="${LAGRANTO_ROOT:-${WORK}/01_LAGRANTO/lagranto.era5.test/lagranto.shellscript}"
# IMPORTANT (one-time): set  set interval = "6"  in
#   ${LAGRANTO_ROOT}/goodies/datelist.sh   before stage 04/05.

# ---- Physics knobs (identical meaning to the event pipeline) ----
export HOURS_CSV="${HOURS_CSV:-0,6,12,18}"       # cdo -selhour (6-hourly snapshots)
export REMAP_GRID="${REMAP_GRID:-r720x360}"      # 0.5 deg global remap target (same as event)
export DOMAIN_LATBOX="${DOMAIN_LATBOX:-0,360,0,90}"   # NH cut from the start ("" = global)
export LEVELS="${LEVELS:-}"                       # ""=all 137 levels; e.g. 70/137 (see README caveat)

# ---- Arrival window around each analogue date (asymmetric allowed) ----
#   Trajectories are started from (date - BEFORE) to (date + AFTER), 6-hourly.
#   Default: 4 days before, 3 days after (an 8-day window, date-4 .. date+3).
export ANALOGUE_DAYS_BEFORE="${ANALOGUE_DAYS_BEFORE:-4}"
export ANALOGUE_DAYS_AFTER="${ANALOGUE_DAYS_AFTER:-3}"
# Back-trajectory length (days of P/S data needed before the earliest arrival):
export ANALOGUE_LEAD_DAYS="${ANALOGUE_LEAD_DAYS:-15}"

# ---- Default worker counts per stage (overridable per sbatch) ----
export DAY_JOBS="${DAY_JOBS:-4}"                  # stage 01
export CDO_THREADS="${CDO_THREADS:-8}"            # stages 01/02
export PREPROCESS_JOBS="${PREPROCESS_JOBS:-8}"    # stage 02
export LIGHT_JOBS="${LIGHT_JOBS:-28}"            # stage 02
export SLOT_JOBS="${SLOT_JOBS:-16}"              # stage 03
export LAGRANTO_JOBS="${LAGRANTO_JOBS:-32}"       # stage 04 (passed as JOBS)
export POST_JOBS="${POST_JOBS:-16}"              # stage 05 (passed as JOBS)
export ORIGIN_JOBS="${ORIGIN_JOBS:-8}"           # stage 06 (passed as JOBS)

mkdir -p "${ANALOGUE_FINAL_BASE}" "${ANALOGUE_ORIGIN_DIR}" "${ANALOGUE_LOG_DIR}" "${ANALOGUE_SCRATCH_BASE}" 2>/dev/null || true
