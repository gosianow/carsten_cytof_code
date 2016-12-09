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
R CMD BATCH --no-save --no-restore $RCODE/10_patient_data.R $ROUT/10_patient_data.Rout
tail $ROUT/10_patient_data.Rout




















#
