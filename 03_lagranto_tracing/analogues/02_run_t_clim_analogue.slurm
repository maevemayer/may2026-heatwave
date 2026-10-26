#!/bin/bash
#SBATCH --account=uc1275
#SBATCH --job-name=tclim-analogue
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=58G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=FAIL

set -euo pipefail

: "${ANALOGUE_RANK:?ANALOGUE_RANK must be set (rank into analogue_dates.txt)}"
: "${WORKFLOW_SCRIPT:?WORKFLOW_SCRIPT must point to t_clim_slurm_6h.py}"

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

PREPROCESS_JOBS="${PREPROCESS_JOBS:-8}"
LIGHT_JOBS="${LIGHT_JOBS:-28}"
CDO_THREADS="${CDO_THREADS:-8}"

if (( PREPROCESS_JOBS * CDO_THREADS > SLURM_CPUS_PER_TASK )); then
    echo "PREPROCESS_JOBS * CDO_THREADS exceeds allocated CPUs" >&2
    exit 2
fi
if (( LIGHT_JOBS > SLURM_CPUS_PER_TASK )); then
    echo "LIGHT_JOBS exceeds allocated CPUs" >&2
    exit 2
fi

echo "Job ID: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Analogue rank: ${ANALOGUE_RANK}"
echo "Preprocess workers: ${PREPROCESS_JOBS} x ${CDO_THREADS} CDO threads"
echo "Light workers: ${LIGHT_JOBS}"

"${PYTHON}" "${WORKFLOW_SCRIPT}" \
    --rank "${ANALOGUE_RANK}" \
    --preprocess-jobs "${PREPROCESS_JOBS}" \
    --light-jobs "${LIGHT_JOBS}" \
    --cdo-threads "${CDO_THREADS}"
