#!/bin/bash
#SBATCH --account=uc1275
#SBATCH --job-name=lagranto-analogue
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=120G
#SBATCH --time=08:00:00
#SBATCH --mail-type=FAIL

set -euo pipefail

: "${ANALOGUE_RANK:?ANALOGUE_RANK must be set (rank into analogue_dates.txt)}"
: "${WORKFLOW_SCRIPT:?WORKFLOW_SCRIPT must point to lagranto_slurm_6h.py}"

# Config is inherited from config_analogues.sh via the submit script (--export=ALL);
# if CONFIG_FILE is set we (re)source it so a standalone sbatch also works.
[[ -n "${CONFIG_FILE:-}" && -f "${CONFIG_FILE}" ]] && source "${CONFIG_FILE}"
PYTHON="${ANALOGUE_PYTHON:-${ANALOGUE_BASE}/kernel/analogues-2026/bin/python}"

if [[ ! -x "${PYTHON}" ]]; then
    echo "Python environment not found: ${PYTHON}" >&2
    exit 2
fi
if [[ ! -f "${WORKFLOW_SCRIPT}" ]]; then
    echo "Workflow script not found: ${WORKFLOW_SCRIPT}" >&2
    exit 2
fi

# Each LAGRANTO trajectory runs serial; parallelism is across start times via
# joblib, so pin every worker (and its caltra/trace children) to one thread.
export OMP_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

JOBS="${JOBS:-32}"
if (( JOBS > SLURM_CPUS_PER_TASK )); then
    echo "JOBS exceeds allocated CPUs" >&2
    exit 2
fi

echo "Job ID: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Analogue rank: ${ANALOGUE_RANK}"
echo "Parallel trajectories: ${JOBS}"

"${PYTHON}" "${WORKFLOW_SCRIPT}" \
    --rank "${ANALOGUE_RANK}" \
    --jobs "${JOBS}"
