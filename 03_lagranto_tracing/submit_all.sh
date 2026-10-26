#!/bin/bash
# ===========================================================================
# Master submit - chains the 6-hourly pipeline with Slurm dependencies.
#
#   01 prepare P-files        \
#                       >-- both --> 03 S-files --> 04 LAGRANTO --> 05 decomp --> 06 origin
#   02 T climatology (tbar)   /
#
# DATA_START/DATA_END : period for which P/S files are built. Must cover
#                       [TRAJ_START - 15 days , TRAJ_END] (15-day back traj).
# TRAJ_START/TRAJ_END : the 6-hourly trajectory start times (stage 04).
#
# Usage:
#   ./submit_all.sh DATA_START DATA_END TRAJ_START TRAJ_END
# Example (ERA5T 2026 event, starts 2026-05-10..06-01, 15 d back => data from 2026-04-25):
#   ./submit_all.sh 2026-04-25 2026-06-01 2026-05-10 2026-06-01
# ===========================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_DIR="${HERE}"   # stage scripts locate config.sh / *.py via this

if [[ $# -ne 4 ]]; then
    echo "usage: $0 DATA_START DATA_END TRAJ_START TRAJ_END"; exit 1
fi
DATA_START="$1"; DATA_END="$2"; TRAJ_START="$3"; TRAJ_END="$4"

source "${HERE}/config.sh"
has() { [[ " ${STAGES} " == *" $1 "* ]]; }
echo "stages: ${STAGES} | P/S read from ${PS_DATA_DIR}"
if [[ "${PS_DATA_DIR}" != "${DATA_DIR}" ]] && { has 01 || has 03; }; then
    echo "PS_DATA_DIR (${PS_DATA_DIR}) differs from DATA_DIR: stages 01/03 would not feed stage 04."; exit 1
fi

# climatology base period = data period +/- 1 day (for the +/-6 h neighbours)
BASE_START="$(date -d "${DATA_START} -1 day" +%Y-%m-%d)"
BASE_END="$(date -d "${DATA_END} +1 day" +%Y-%m-%d)"

echo "P/S data period : ${DATA_START} .. ${DATA_END}"
echo "climatology     : ${BASE_START} .. ${BASE_END}"
echo "trajectory runs : ${TRAJ_START} .. ${TRAJ_END}"

# ---- without stage 03: link the existing (correct) S-files and check coverage ----
if ! has 03; then
    if [[ -n "${S_REUSE_DIR}" && "${S_REUSE_DIR}" != "${PS_DATA_DIR}" ]]; then
        n=0
        for f in "${S_REUSE_DIR}"/S[0-9]*; do
            [[ -e "$f" ]] || continue                       # no match: pattern stays literal
            [[ -e "${PS_DATA_DIR}/$(basename "$f")" ]] || { ln -s "$f" "${PS_DATA_DIR}/"; n=$((n + 1)); }
        done
        echo "linked ${n} S-files from ${S_REUSE_DIR}"
    fi
    missing=0; d="${DATA_START}"
    while [[ ! "$d" > "${DATA_END}" ]]; do
        for hh in 00 06 12 18; do [[ -e "${PS_DATA_DIR}/S$(date -d "$d" +%Y%m%d)_${hh}" ]] || missing=$((missing + 1)); done
        d="$(date -d "$d +1 day" +%Y-%m-%d)"
    done
    if (( missing > 0 )); then
        echo "${missing} S-files for ${DATA_START}..${DATA_END} missing in ${PS_DATA_DIR} (and stage 03 is not submitted)"; exit 1
    fi
fi

j01=""; j02=""; j03=""
if has 01; then
    j01=$(sbatch --parsable --export=ALL,START="${DATA_START}",END="${DATA_END}" "${HERE}/01_prepare_era5.slurm")
    echo "  01 prepare       -> ${j01}"
fi
if has 02; then
    j02=$(sbatch --parsable --export=ALL,BASE_START="${BASE_START}",BASE_END="${BASE_END}" "${HERE}/02_t_climatology.slurm")
    echo "  02 climatology   -> ${j02}"
fi
if has 03; then
    dep=$(echo "${j01} ${j02}" | xargs | tr ' ' ':')
    j03=$(sbatch --parsable ${dep:+--dependency=afterok:${dep}} \
          --export=ALL,START="${DATA_START}",END="${DATA_END}" "${HERE}/03_build_s_files.slurm")
    echo "  03 S-files       -> ${j03}  (after ${dep:-nothing})"
fi
prev="${j03:-${j01}}"
if has 04; then
    j04=$(sbatch --parsable ${prev:+--dependency=afterok:${prev}} \
          --export=ALL,START="${TRAJ_START}",END="${TRAJ_END}" "${HERE}/04_run_lagranto.slurm")
    echo "  04 LAGRANTO      -> ${j04}  (after ${prev:-nothing})"; prev="${j04}"
fi
if has 05; then
    j05=$(sbatch --parsable ${prev:+--dependency=afterok:${prev}} "${HERE}/05_traj_postprocess.slurm")
    echo "  05 decomposition -> ${j05}  (after ${prev:-nothing})"; prev="${j05}"
fi
if has 06; then
    j06=$(sbatch --parsable ${prev:+--dependency=afterok:${prev}} "${HERE}/06_origin_contributions.slurm")
    echo "  06 origin        -> ${j06}  (after ${prev:-nothing})"
fi

echo "Submitted. Watch with:  squeue -u \$USER"
