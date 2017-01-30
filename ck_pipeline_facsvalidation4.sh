#!/bin/bash

shopt -s expand_aliases
source ~/.bash_aliases

###############################################################################################################
## Define paths to software and reference files

# RWD_MAIN=/Users/gosia/Dropbox/UZH/carsten_cytof
# RCODE=/Users/gosia/Dropbox/UZH/carsten_cytof_code
RWD_MAIN=/home/Shared/data/cytof/carsten_cytof
RCODE=/home/gosia/R/carsten_cytof_code

###############################################################################################################
# Analysis of Patientdata
###############################################################################################################

data_dir="FACSvalidation4"

RWD=$RWD_MAIN/${data_dir}
ROUT=$RWD/Rout
mkdir -p $ROUT

prefix='facs_valid3_'
metadata='CK_2017-01-27_ValidationDataset2_metadata.xlsx'
freqs='CK_2017-01-27_ValidationDataset2_frequencies.xlsx'

echo ">>> 10_facs_validation"
R33 CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${prefix}' freq_outdir='ck_analysis' path_metadata='ck_orig_files/${metadata}' path_freqs='ck_orig_files/${freqs}' path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_1dataset_3responses_base.R' path_fun_plot_heatmaps='$RCODE/00_plot_heatmaps_for_sign_freqs.R' pdf_hight=4" $RCODE/10_facs_validation.R $ROUT/10_facs_validation.Rout
tail $ROUT/10_facs_validation.Rout




















#
