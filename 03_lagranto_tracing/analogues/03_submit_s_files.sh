#!/bin/bash
set -euo pipefail

# Submit one analogue per job in a strict dependency chain, mirroring
# 01_submit_prepare_era5.sh / 02_submit_t_clim.sh. Each S file needs a
# finished P file and a finished t_clim file for the same date/slot, so run
# this only after both of those chains have completed for the ranks below
# (or pass AFTER_JOBID to chain onto the last job of whichever ran last).

ACCOUNT="${ACCOUNT:-uc1275}"
START_RANK="${START_RANK:-1}"
END_RANK="${END_RANK:-0}"   # 0 = all analogues (count from analogue_dates.txt)
SLOT_JOBS="${SLOT_JOBS:-16}"
AFTER_JOBID="${AFTER_JOBID:-}"

_SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${CONFIG_FILE:-${_SD}/config_analogues.sh}"   # central paths/account/workers
PROJECT_ROOT="${ANALOGUE_HOME}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_SCRIPT="${WORKFLOW_SCRIPT:-${SCRIPT_DIR}/03_s_files_slurm_6h.py}"
SLURM_SCRIPT="${SLURM_SCRIPT:-${SCRIPT_DIR}/03_run_s_files_analogue.slurm}"
LOG_DIR="${ANALOGUE_LOG_DIR}"
FINAL_DIR="${ANALOGUE_FINAL_BASE}" 
#"${PROJECT_ROOT}/lagranto"

N_ANALOGUES=$(grep -cE '^[[:space:]]*[0-9]{4}-[0-9]{2}-[0-9]{2}' "${SCRIPT_DIR}/analogue_dates.txt")
(( END_RANK == 0 )) && END_RANK="${N_ANALOGUES}"
if (( START_RANK < 1 || END_RANK > N_ANALOGUES || START_RANK > END_RANK )); then
    echo "Ranks must satisfy 1 <= START_RANK <= END_RANK <= ${N_ANALOGUES}" >&2
    exit 2
fi
if [[ ! -f "${WORKFLOW_SCRIPT}" || ! -f "${SLURM_SCRIPT}" ]]; then
    echo "Keep the submit, Slurm, and Python files in the same directory." >&2
    exit 2
fi

mkdir -p "${LOG_DIR}"

dependency="${AFTER_JOBID}"
for rank in $(seq "${START_RANK}" "${END_RANK}"); do
    arguments=(
        --parsable
        --account="${ACCOUNT}"
        --job-name="sfiles-r${rank}"
        --output="${LOG_DIR}/sfiles-rank${rank}-%j.out"
        --export="ALL,ANALOGUE_RANK=${rank},WORKFLOW_SCRIPT=${WORKFLOW_SCRIPT},SLOT_JOBS=${SLOT_JOBS}"
    )
    if [[ -n "${dependency}" ]]; then
        arguments+=(--dependency="afterok:${dependency}")
    fi

    job_id="$(sbatch "${arguments[@]}" "${SLURM_SCRIPT}")"
    echo "Submitted rank ${rank}: ${job_id}"
    dependency="${job_id}"
done

echo
echo "Final chained job: ${dependency}"
echo "Final S/P-file directory: ${FINAL_DIR}"
echo "Logs: ${LOG_DIR}"
