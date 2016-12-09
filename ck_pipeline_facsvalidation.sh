#!/bin/bash

###############################################################################################################
## Define paths to software and reference files

RWD_MAIN=/Users/gosia/Dropbox/UZH/carsten_cytof
RCODE=/Users/gosia/Dropbox/UZH/carsten_cytof_code

###############################################################################################################
# Analysis of Patientdata
###############################################################################################################

data_dir="FACS_validation"

RWD=$RWD_MAIN/${data_dir}
ROUT=$RWD/Rout
mkdir -p $ROUT


echo ">>> 10_facs_validation"
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='facs_valid_' freq_outdir='ck_analysis' path_metadata='ck_orig_files/FACSvalidation_01_metadata.xlsx'  path_freqs='ck_orig_files/FACSvalidation_01_frequencies.xlsx'   path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_1dataset_3responses_base.R' pdf_hight=4" $RCODE/10_facs_validation.R $ROUT/10_facs_validation.Rout
tail $ROUT/10_facs_validation.Rout




















#
