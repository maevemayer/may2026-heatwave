#!/bin/bash

# ==========================================
#              SBATCH SETTINGS
# ==========================================

#SBATCH --job-name=aifs_run          # name of job
#SBATCH -C a100
#SBATCH -A udt@a100
##SBATCH --qos=qos_gpu_a100-dev      # uncomment for fast debug, time<2h
#SBATCH --nodes=1                    # we request one node
#SBATCH --ntasks-per-node=1          # with one task per node (= number of GPUs here)
#SBATCH --gres=gpu:1                 # number of GPUs per node (max 8 with gpu_p2, gpu_p5)
#SBATCH --cpus-per-task=8           # number of cores per task (1/4 of the 4-GPUs V100 node)
#SBATCH --hint=nomultithread         # hyperthreading is deactivated
#SBATCH --time=04:00:00              # maximum execution time requested (HH:MM:SS)
#SBATCH --output=out%j.out    # name of output file
#SBATCH --error=out%j.out     # name of error file (here, in common with the output file)

# ==========================================
#              PREPARE ENV
# ==========================================

# Cleans out the modules loaded in interactive and inherited by default 
module purge
module load pytorch-gpu/py3/2.8.0

# Echo of launched commands
set -x

# ==========================================
#          EXPERIMENT PARAMETERS
# ==========================================

EXP_NAME="may_2026_heatwave"
DATE_STR="2026051500"

N_MEMBERS=200
LEAD=480

# Surface variables
SURF_VARS="2t 2d sp tp"
# Temperature
T_VARS="t_200 t_250 t_300 t_400 t_500 t_600 t_700 t_850 t_925 t_1000"
# Geopotential
Z_VARS="z_200 z_250 z_300 z_400 z_500 z_600 z_700 z_850 z_925 z_1000"
# U-component of wind
U_VARS="u_200 u_250 u_300 u_400 u_500 u_600 u_700 u_850 u_925 u_1000"
# V-component of wind
V_VARS="v_200 v_250 v_300 v_400 v_500 v_600 v_700 v_850 v_925 v_1000"
# Vertical velocity (omega)
W_VARS="w_200 w_250 w_300 w_400 w_500 w_600 w_700 w_850 w_925 w_1000"

VARIABLES="$SURF_VARS $T_VARS $Z_VARS $U_VARS $V_VARS $W_VARS" #variables to save

SAVETIMES="0 6 12 18" #at which times to save variables

RUN_IN_PATH="/path/to/inputs_aifs_ensemble_boosting"
RUN_OUT_PATH="/path/to/outputs_aifs_ensemble_boosting"
CHECKPOINT="/path/to/aifs-ens-crps-1.0.ckpt"   # AIFS-ENS CRPS v1.0 weights (see README)
MASK="/path/to/mask_eu.pkl"


# ==========================================
#              RUN AIFS
# ==========================================

python run_aifs_crps_simple_boosting.py \
    --input_path "$RUN_IN_PATH" \
    --output_path "$RUN_OUT_PATH" \
    --date_str "$DATE_STR" \
    --exp_name "$EXP_NAME" \
    --checkpoint "$CHECKPOINT" \
    --n_members $N_MEMBERS \
    --variables $VARIABLES \
    --savetimes $SAVETIMES \
    --lead_time $LEAD \
    --mask "$MASK"
