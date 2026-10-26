#!/bin/bash
set -euo pipefail

# Submit the origin-contribution build for one analogue per job in a dependency
# chain. Lightweight; run after (or chained onto) 05_submit_traj_postprocess.sh.

ACCOUNT="${ACCOUNT:-uc1275}"
START_RANK="${START_RANK:-1}"
END_RANK="${END_RANK:-0}"   # 0 = all analogues (count from analogue_dates.txt)
JOBS="${JOBS:-8}"

_SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${CONFIG_FILE:-${_SD}/config_analogues.sh}"   # central paths/account/workers
PROJECT_ROOT="${ANALOGUE_HOME}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_SCRIPT="${WORKFLOW_SCRIPT:-${SCRIPT_DIR}/06_origin_contributions_slurm_6h.py}"
SLURM_SCRIPT="${SLURM_SCRIPT:-${SCRIPT_DIR}/06_run_origin_contributions_analogue.slurm}"
LOG_DIR="${ANALOGUE_LOG_DIR}"

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

dependency=""
for rank in $(seq "${START_RANK}" "${END_RANK}"); do
    arguments=(
        --parsable
        --account="${ACCOUNT}"
        --job-name="origin-r${rank}"
        --output="${LOG_DIR}/origin-rank${rank}-%j.out"
        --export="ALL,ANALOGUE_RANK=${rank},WORKFLOW_SCRIPT=${WORKFLOW_SCRIPT},JOBS=${JOBS}"
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
echo "Origin-contribution files: ${ANALOGUE_ORIGIN_DIR}/rank<NN>_MJ_<date>/origin_contributions.nc"
echo "Logs: ${LOG_DIR}"
