#!/bin/bash
#SBATCH --account=uc1275
#SBATCH --job-name=sfiles-analogue
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --mail-type=FAIL

set -euo pipefail

: "${ANALOGUE_RANK:?ANALOGUE_RANK must be set (rank into analogue_dates.txt)}"
: "${WORKFLOW_SCRIPT:?WORKFLOW_SCRIPT must point to s_files_slurm_6h.py}"

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

export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

SLOT_JOBS="${SLOT_JOBS:-16}"
SKIP_RENAME="${SKIP_RENAME:-0}"

if (( SLOT_JOBS > SLURM_CPUS_PER_TASK )); then
    echo "SLOT_JOBS exceeds allocated CPUs" >&2
    exit 2
fi

RENAME_ARGS=()
if [[ "${SKIP_RENAME}" != "0" ]]; then
    RENAME_ARGS+=(--skip-rename)
fi

echo "Job ID: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Analogue rank: ${ANALOGUE_RANK}"
echo "Slot workers: ${SLOT_JOBS}"
echo "Filename rename to synoptic schedule: $([[ "${SKIP_RENAME}" != "0" ]] && echo disabled || echo enabled)"

"${PYTHON}" "${WORKFLOW_SCRIPT}" \
    --rank "${ANALOGUE_RANK}" \
    --slot-jobs "${SLOT_JOBS}" \
    ${RENAME_ARGS[@]+"${RENAME_ARGS[@]}"}
