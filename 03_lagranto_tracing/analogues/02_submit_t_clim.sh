#!/bin/bash
set -euo pipefail

# Submit one analogue per job in a strict dependency chain. This keeps only
# one active scratch tree and stops the sequence if any analogue fails.

ACCOUNT="${ACCOUNT:-uc1275}"
START_RANK="${START_RANK:-1}"
END_RANK="${END_RANK:-0}"   # 0 = all analogues (count from analogue_dates.txt)
PREPROCESS_JOBS="${PREPROCESS_JOBS:-8}"
LIGHT_JOBS="${LIGHT_JOBS:-28}"
CDO_THREADS="${CDO_THREADS:-8}"

_SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${CONFIG_FILE:-${_SD}/config_analogues.sh}"   # central paths/account/workers
PROJECT_ROOT="${ANALOGUE_HOME}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_SCRIPT="${WORKFLOW_SCRIPT:-${SCRIPT_DIR}/02_t_clim_slurm_6h.py}"
SLURM_SCRIPT="${SLURM_SCRIPT:-${SCRIPT_DIR}/02_run_t_clim_analogue.slurm}"
LOG_DIR="${ANALOGUE_LOG_DIR}"
FINAL_DIR="${ANALOGUE_FINAL_BASE}/t_clim"

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

mkdir -p "${LOG_DIR}" "${FINAL_DIR}"

dependency=""
for rank in $(seq "${START_RANK}" "${END_RANK}"); do
    arguments=(
        --parsable
        --account="${ACCOUNT}"
        --job-name="tclim-r${rank}"
        --output="${LOG_DIR}/tclim-rank${rank}-%j.out"
        --export="ALL,ANALOGUE_RANK=${rank},WORKFLOW_SCRIPT=${WORKFLOW_SCRIPT},PREPROCESS_JOBS=${PREPROCESS_JOBS},LIGHT_JOBS=${LIGHT_JOBS},CDO_THREADS=${CDO_THREADS}"
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
echo "Final climatologies: ${FINAL_DIR}"
echo "Logs: ${LOG_DIR}"
