#!/bin/bash
# ===========================================================================
# Global submit for the ANALOGUE 6-hourly pipeline — the analogue counterpart
# of ../submit_all.sh. For every rank it chains all six stages with Slurm
# dependencies (01 & 02 independent -> 03 -> 04 -> 05 -> 06), exactly the same
# stage graph as the 2026 event.
#
#   source config_analogues.sh  (all paths / account / workers live there)
#
# Usage:
#   ./submit_all_analogues.sh                 # all analogues in analogue_dates.txt
#   START_RANK=1 END_RANK=5 ./submit_all_analogues.sh
#   SERIAL_RANKS=1 ./submit_all_analogues.sh  # chain rank N after rank N-1 (bounds scratch)
# ===========================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${CONFIG_FILE:-${SCRIPT_DIR}/config_analogues.sh}"
# shellcheck source=/dev/null
source "${CONFIG_FILE}"

START_RANK="${START_RANK:-1}"
END_RANK="${END_RANK:-0}"          # 0 = all analogues (count from analogue_dates.txt)
SERIAL_RANKS="${SERIAL_RANKS:-0}"  # 1 = make each rank wait for the previous one's stage 06
has() { [[ " ${STAGES} " == *" $1 "* ]]; }
if [[ "${ANALOGUE_PS_BASE}" != "${ANALOGUE_FINAL_BASE}" ]] && { has 01 || has 03; }; then
    echo "ANALOGUE_PS_BASE differs from ANALOGUE_FINAL_BASE: stages 01/03 would not feed stage 04." >&2; exit 2
fi
echo "stages: ${STAGES} | P/S read from ${ANALOGUE_PS_BASE}"
DATES=($(grep -E '^[[:space:]]*[0-9]{4}-[0-9]{2}-[0-9]{2}' "${SCRIPT_DIR}/analogue_dates.txt" | tr -d ' '))

N_ANALOGUES=$(grep -cE '^[[:space:]]*[0-9]{4}-[0-9]{2}-[0-9]{2}' "${SCRIPT_DIR}/analogue_dates.txt")
(( END_RANK == 0 )) && END_RANK="${N_ANALOGUES}"
if (( START_RANK < 1 || END_RANK > N_ANALOGUES || START_RANK > END_RANK )); then
    echo "Ranks must satisfy 1 <= START_RANK <= END_RANK <= ${N_ANALOGUES}" >&2
    exit 2
fi

mkdir -p "${ANALOGUE_LOG_DIR}"

# stage -> (python workflow script, slurm wrapper)
PY01="${SCRIPT_DIR}/01_prepare_era5_slurm_6h.py";     SL01="${SCRIPT_DIR}/01_run_prepare_era5_analogue.slurm"
PY02="${SCRIPT_DIR}/02_t_clim_slurm_6h.py";           SL02="${SCRIPT_DIR}/02_run_t_clim_analogue.slurm"
PY03="${SCRIPT_DIR}/03_s_files_slurm_6h.py";          SL03="${SCRIPT_DIR}/03_run_s_files_analogue.slurm"
PY04="${SCRIPT_DIR}/04_lagranto_slurm_6h.py";         SL04="${SCRIPT_DIR}/04_run_lagranto_analogue.slurm"
PY05="${SCRIPT_DIR}/05_traj_postprocess_slurm_6h.py"; SL05="${SCRIPT_DIR}/05_run_traj_postprocess_analogue.slurm"
PY06="${SCRIPT_DIR}/06_origin_contributions_slurm_6h.py"; SL06="${SCRIPT_DIR}/06_run_origin_contributions_analogue.slurm"
for f in "$PY01" "$SL01" "$PY02" "$SL02" "$PY03" "$SL03" "$PY04" "$SL04" "$PY05" "$SL05" "$PY06" "$SL06" \
         "${SCRIPT_DIR}/04_run_lagranto_analogue.sh"; do
    [[ -f "$f" ]] || { echo "Missing required file: $f" >&2; exit 2; }
done

COMMON="ALL,CONFIG_FILE=${CONFIG_FILE}"
sb() {  # sb <dependency-or-empty> <jobname> <extra-export> <slurm>  -> prints job id
    local dep="$1" name="$2" extra="$3" slurm="$4"
    local args=(--parsable --account="${ANALOGUE_ACCOUNT}" --partition="${ANALOGUE_PARTITION}"
                --job-name="${name}" --output="${ANALOGUE_LOG_DIR}/${name}-%j.out"
                --export="${COMMON},${extra}")
    [[ -n "${dep}" ]] && args+=(--dependency="afterok:${dep}")
    sbatch "${args[@]}" "${slurm}"
}

prev_rank_tail=""
for rank in $(seq "${START_RANK}" "${END_RANK}"); do
    echo "=== rank ${rank} ==="
    gate="${prev_rank_tail}"   # empty unless SERIAL_RANKS
    j01=""; j02=""; j03=""; j04=""; j05=""; j06=""
    if ! has 03; then        # link this rank's existing S-files (found by date) into the new tree
        ymd="$(echo "${DATES[$((rank - 1))]}" | tr -d '-')"
        tag="$(printf 'rank%02d_MJ_%s' "${rank}" "${ymd}")"
        src=( "${ANALOGUE_S_REUSE_BASE}"/rank*_MJ_"${ymd}" )
        if [[ ! -d "${src[0]}" || ${#src[@]} -ne 1 ]]; then
            echo "  no unique S-file folder for ${ymd} in ${ANALOGUE_S_REUSE_BASE}: ${src[*]}" >&2; exit 2
        fi
        mkdir -p "${ANALOGUE_PS_BASE}/${tag}"; n=0
        for f in "${src[0]}"/S[0-9]*; do
            [[ -e "$f" ]] || continue
            [[ -e "${ANALOGUE_PS_BASE}/${tag}/$(basename "$f")" ]] || { ln -s "$f" "${ANALOGUE_PS_BASE}/${tag}/"; n=$((n + 1)); }
        done
        nS=$(ls "${ANALOGUE_PS_BASE}/${tag}" | grep -c '^S[0-9]')
        echo "  ${tag}: linked ${n} S-files from $(basename "${src[0]}") (${nS} present)"
        (( nS > 0 )) || { echo "  no S-files for ${tag}" >&2; exit 2; }
    fi
    if has 01; then j01=$(sb "${gate}" "an01-r${rank}" "ANALOGUE_RANK=${rank},WORKFLOW_SCRIPT=${PY01},DAY_JOBS=${DAY_JOBS},CDO_THREADS=${CDO_THREADS}" "${SL01}"); fi
    if has 02; then j02=$(sb "${gate}" "an02-r${rank}" "ANALOGUE_RANK=${rank},WORKFLOW_SCRIPT=${PY02},PREPROCESS_JOBS=${PREPROCESS_JOBS},LIGHT_JOBS=${LIGHT_JOBS},CDO_THREADS=${CDO_THREADS}" "${SL02}"); fi
    if has 03; then
        dep=$(echo "${j01} ${j02}" | xargs | tr ' ' ':'); dep="${dep:-${gate}}"
        j03=$(sb "${dep}" "an03-r${rank}" "ANALOGUE_RANK=${rank},WORKFLOW_SCRIPT=${PY03},SLOT_JOBS=${SLOT_JOBS}" "${SL03}")
    fi
    prev="${j03:-${j01:-${gate}}}"
    if has 04; then j04=$(sb "${prev}" "an04-r${rank}" "ANALOGUE_RANK=${rank},WORKFLOW_SCRIPT=${PY04},JOBS=${LAGRANTO_JOBS}" "${SL04}"); prev="${j04}"; fi
    if has 05; then j05=$(sb "${prev}" "an05-r${rank}" "ANALOGUE_RANK=${rank},WORKFLOW_SCRIPT=${PY05},JOBS=${POST_JOBS}" "${SL05}"); prev="${j05}"; fi
    if has 06; then j06=$(sb "${prev}" "an06-r${rank}" "ANALOGUE_RANK=${rank},WORKFLOW_SCRIPT=${PY06},JOBS=${ORIGIN_JOBS}" "${SL06}"); prev="${j06}"; fi
    echo "  01=${j01} 02=${j02} 03=${j03} 04=${j04} 05=${j05} 06=${j06}"
    (( SERIAL_RANKS == 1 )) && prev_rank_tail="${prev}"
done

echo
echo "Submitted ranks ${START_RANK}..${END_RANK}."
echo "P/S read: ${ANALOGUE_PS_BASE}/rank<NN>_MJ_<date>/"
echo "traj    : ${ANALOGUE_FINAL_BASE}/rank<NN>_MJ_<date>/   (results_6h, traj_netcdfs)"
echo "Origin  : ${ANALOGUE_ORIGIN_DIR}/rank<NN>_MJ_<date>/origin_contributions.nc"
echo "Logs    : ${ANALOGUE_LOG_DIR}"
echo "Watch   : squeue -u \$USER"
