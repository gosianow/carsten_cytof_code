#!/bin/bash

###############################################################################################################
## Define paths to software and reference files

RWD_MAIN=/Users/gosia/Dropbox/UZH/carsten_cytof
RCODE=/Users/gosia/Dropbox/UZH/carsten_cytof_code

###############################################################################################################
# Analysis of Patientdata
###############################################################################################################

data_dir="Patientdata"

RWD=$RWD_MAIN/${data_dir}
ROUT=$RWD/Rout
mkdir -p $ROUT


echo ">>> 10_patient_data"
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' outdir='ck_analysis' path_data='ck_orig_files/Patient_data_cytof_FACS_09_12_16_CK_GN.xlsx' path_variables='ck_orig_files/variables_to_use.xlsx' path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_1dataset_2responses_base.R' path_fun_plot_heatmaps='$RCODE/00_plot_heatmaps_for_sign_expr.R'" $RCODE/10_patient_data.R $ROUT/10_patient_data.Rout
tail $ROUT/10_patient_data.Rout


















#
