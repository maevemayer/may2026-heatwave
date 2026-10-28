#!/bin/bash

# ==========================================
#              SBATCH SETTINGS
# ==========================================

#SBATCH --partition=zen4
#SBATCH --time=04:00:00
#SBATCH --mem=64G

# ==========================================
#              PREPARE ENV
# ==========================================

module purge
module load pangeo-meso/2026.01.21

conda activate /path/to/your/pangeo-earthkit-env

# =========================================================================================================================

# ==========================================
#          EXPERIMENT PARAMETERS
# ==========================================

EXP_NAME="may_2026_heatwave"
DATE_INPUT="2026-05-15 00:00:00"
INPUTS_SAVE="/path/to/AIFS_INPUTS/"

# ==========================================
#              RUN PREPROCESS
# ==========================================

python aifs_inputs_dl+preprocess.py \
	--date_input "$DATE_INPUT" \
	--exp_name "$EXP_NAME" \
	--output_save "$INPUTS_SAVE"


# ==========================================
#          EXPERIMENT PARAMETERS
# ==========================================

EXP_NAME="may_2026_heatwave"
DATE_INPUT="2026-05-16 00:00:00"
INPUTS_SAVE="/path/to/AIFS_INPUTS/"

# ==========================================
#              RUN PREPROCESS
# ==========================================

python aifs_inputs_dl+preprocess.py \
	--date_input "$DATE_INPUT" \
	--exp_name "$EXP_NAME" \
	--output_save "$INPUTS_SAVE"



# ==========================================
#          EXPERIMENT PARAMETERS
# ==========================================

EXP_NAME="may_2026_heatwave"
DATE_INPUT="2026-05-17 00:00:00"
INPUTS_SAVE="/path/to/AIFS_INPUTS/"

# ==========================================
#              RUN PREPROCESS
# ==========================================

python aifs_inputs_dl+preprocess.py \
	--date_input "$DATE_INPUT" \
	--exp_name "$EXP_NAME" \
	--output_save "$INPUTS_SAVE"



# ==========================================
#          EXPERIMENT PARAMETERS
# ==========================================

EXP_NAME="may_2026_heatwave"
DATE_INPUT="2026-05-18 00:00:00"
INPUTS_SAVE="/path/to/AIFS_INPUTS/"

# ==========================================
#              RUN PREPROCESS
# ==========================================

python aifs_inputs_dl+preprocess.py \
	--date_input "$DATE_INPUT" \
	--exp_name "$EXP_NAME" \
	--output_save "$INPUTS_SAVE"



# ==========================================
#          EXPERIMENT PARAMETERS
# ==========================================

EXP_NAME="may_2026_heatwave"
DATE_INPUT="2026-05-19 00:00:00"
INPUTS_SAVE="/path/to/AIFS_INPUTS/"

# ==========================================
#              RUN PREPROCESS
# ==========================================

python aifs_inputs_dl+preprocess.py \
	--date_input "$DATE_INPUT" \
	--exp_name "$EXP_NAME" \
	--output_save "$INPUTS_SAVE"



# ==========================================
#          EXPERIMENT PARAMETERS
# ==========================================

EXP_NAME="may_2026_heatwave"
DATE_INPUT="2026-05-20 00:00:00"
INPUTS_SAVE="/path/to/AIFS_INPUTS/"

# ==========================================
#              RUN PREPROCESS
# ==========================================

python aifs_inputs_dl+preprocess.py \
	--date_input "$DATE_INPUT" \
	--exp_name "$EXP_NAME" \
	--output_save "$INPUTS_SAVE"



# ==========================================
#          EXPERIMENT PARAMETERS
# ==========================================

EXP_NAME="may_2026_heatwave"
DATE_INPUT="2026-05-21 00:00:00"
INPUTS_SAVE="/path/to/AIFS_INPUTS/"

# ==========================================
#              RUN PREPROCESS
# ==========================================

python aifs_inputs_dl+preprocess.py \
	--date_input "$DATE_INPUT" \
	--exp_name "$EXP_NAME" \
	--output_save "$INPUTS_SAVE"



# ==========================================
#          EXPERIMENT PARAMETERS
# ==========================================

EXP_NAME="may_2026_heatwave"
DATE_INPUT="2026-05-22 00:00:00"
INPUTS_SAVE="/path/to/AIFS_INPUTS/"

# ==========================================
#              RUN PREPROCESS
# ==========================================

python aifs_inputs_dl+preprocess.py \
	--date_input "$DATE_INPUT" \
	--exp_name "$EXP_NAME" \
	--output_save "$INPUTS_SAVE"



# ==========================================
#          EXPERIMENT PARAMETERS
# ==========================================

EXP_NAME="may_2026_heatwave"
DATE_INPUT="2026-05-23 00:00:00"
INPUTS_SAVE="/path/to/AIFS_INPUTS/"

# ==========================================
#              RUN PREPROCESS
# ==========================================

python aifs_inputs_dl+preprocess.py \
	--date_input "$DATE_INPUT" \
	--exp_name "$EXP_NAME" \
	--output_save "$INPUTS_SAVE"



# ==========================================
#          EXPERIMENT PARAMETERS
# ==========================================

EXP_NAME="may_2026_heatwave"
DATE_INPUT="2026-05-24 00:00:00"
INPUTS_SAVE="/path/to/AIFS_INPUTS/"

# ==========================================
#              RUN PREPROCESS
# ==========================================

python aifs_inputs_dl+preprocess.py \
	--date_input "$DATE_INPUT" \
	--exp_name "$EXP_NAME" \
	--output_save "$INPUTS_SAVE"



# ==========================================
#          EXPERIMENT PARAMETERS
# ==========================================

EXP_NAME="may_2026_heatwave"
DATE_INPUT="2026-05-25 00:00:00"
INPUTS_SAVE="/path/to/AIFS_INPUTS/"

# ==========================================
#              RUN PREPROCESS
# ==========================================

python aifs_inputs_dl+preprocess.py \
	--date_input "$DATE_INPUT" \
	--exp_name "$EXP_NAME" \
	--output_save "$INPUTS_SAVE"
