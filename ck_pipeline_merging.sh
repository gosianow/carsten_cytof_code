#!/bin/bash

###############################################################################################################
## Define paths to software and reference files

RWD_MAIN=/Users/gosia/Dropbox/UZH/carsten_cytof
RCODE=/Users/gosia/Dropbox/UZH/carsten_cytof_code
METADATA=$RWD_MAIN/CK_metadata
PANELS=$RWD_MAIN/CK_panels

## Define which analysis to re-run
frequencies_merged=true


###############################################################################################################
# Merging data 23 and 29 for panel 1
###############################################################################################################

PANEL=1

new_data_dir="CK_2016-06-merged_23_29/01"

data_name=("data23" "data29")
data_dir=("CK_2016-06-23_01" "CK_2016-06-29_01")
file_metadata=("metadata_23_01.xlsx" "metadata_29_01.xlsx")

prefix_data=("23_" "29_")
prefix_panel="01_"
prefix_pca="pca1_"
prefix_merging=("merging6_" "merging4_")

prefix_freq_merging="23m6_29m4_"


RWD=$RWD_MAIN/${new_data_dir}
ROUT=$RWD/Rout
mkdir -p $ROUT
echo "$RWD"


### Get cluster frequencies
if ${frequencies_merged}; then
  echo "08_frequencies_merged"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${prefix_freq_merging}' freq_outdir='08_frequencies_merged' path_metadata=c('${METADATA}/${file_metadata[0]}','${METADATA}/${file_metadata[1]}')  path_counts=c('${RWD_MAIN}/${data_dir[0]}/050_frequencies/${prefix_data[0]}${prefix_panel}${prefix_pca}${prefix_merging[0]}counts.xls','${RWD_MAIN}/${data_dir[1]}/050_frequencies/${prefix_data[1]}${prefix_panel}${prefix_pca}${prefix_merging[1]}counts.xls') data_name=c('${data_name[0]}','${data_name[1]}') path_fun_models='$RCODE/00_models_merged.R'" $RCODE/08_frequencies_merged.R $ROUT/08_frequencies_merged.Rout
  tail $ROUT/08_frequencies_merged.Rout
fi



###############################################################################################################
# Merging data 23 and 29 for panel 3
###############################################################################################################

PANEL=3

new_data_dir="CK_2016-06-merged_23_29/03"

data_name=("data23" "data29")
data_dir=("CK_2016-06-23_03" "CK_2016-06-29_03")
file_metadata=("metadata_23_03.xlsx" "metadata_29_03.xlsx")

prefix_data=("23_" "29_")
prefix_panel="03_"
prefix_pca="pca1_"
prefix_merging=("merging4_" "merging5_")

prefix_freq_merging="23m4_29m5_"


RWD=$RWD_MAIN/${new_data_dir}
ROUT=$RWD/Rout
mkdir -p $ROUT
echo "$RWD"


# ### Get cluster frequencies
# if ${frequencies_merged}; then
#   echo "08_frequencies_merged"
#   R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${prefix_freq_merging}' freq_outdir='08_frequencies_merged' path_metadata=c('${METADATA}/${file_metadata[0]}','${METADATA}/${file_metadata[1]}')  path_counts=c('${RWD_MAIN}/${data_dir[0]}/050_frequencies/${prefix_data[0]}${prefix_panel}${prefix_pca}${prefix_merging[0]}counts.xls','${RWD_MAIN}/${data_dir[1]}/050_frequencies/${prefix_data[1]}${prefix_panel}${prefix_pca}${prefix_merging[1]}counts.xls') data_name=c('${data_name[0]}','${data_name[1]}') path_fun_models='$RCODE/00_models_merged.R'" $RCODE/08_frequencies_merged.R $ROUT/08_frequencies_merged.Rout
#   tail $ROUT/08_frequencies_merged.Rout
# fi


###############################################################################################################
# Merging data 23 and 29 for panel 1 CD4
###############################################################################################################

PANEL=1

new_data_dir="CK_2016-06-merged_23_29/01_CD4"

data_name=("data23" "data29")
data_dir=("CK_2016-06-23_01_CD4_mergingNEW2" "CK_2016-06-29_01_CD4_merging2")
file_metadata=("metadata_23_01.xlsx" "metadata_29_01.xlsx")

prefix_data=("23CD4_" "29CD4_")
prefix_panel="01CD4_"
prefix_pca="pca1_"
prefix_merging=("merging4_" "merging4_")

prefix_freq_merging="23CD4m4_29CD4m4_"


RWD=$RWD_MAIN/${new_data_dir}
ROUT=$RWD/Rout
mkdir -p $ROUT
echo "$RWD"


# ### Get cluster frequencies
# if ${frequencies_merged}; then
#   echo "08_frequencies_merged"
#   R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${prefix_freq_merging}' freq_outdir='08_frequencies_merged' path_metadata=c('${METADATA}/${file_metadata[0]}','${METADATA}/${file_metadata[1]}')  path_counts=c('${RWD_MAIN}/${data_dir[0]}/050_frequencies/${prefix_data[0]}${prefix_panel}${prefix_pca}${prefix_merging[0]}counts.xls','${RWD_MAIN}/${data_dir[1]}/050_frequencies/${prefix_data[1]}${prefix_panel}${prefix_pca}${prefix_merging[1]}counts.xls') data_name=c('${data_name[0]}','${data_name[1]}') path_fun_models='$RCODE/00_models_merged.R'" $RCODE/08_frequencies_merged.R $ROUT/08_frequencies_merged.Rout
#   tail $ROUT/08_frequencies_merged.Rout
# fi



###############################################################################################################
# Merging data 23 and 29 for panel 1 CD8
###############################################################################################################

PANEL=1

new_data_dir="CK_2016-06-merged_23_29/01_CD8"

data_name=("data23" "data29")
data_dir=("CK_2016-06-23_01_CD8_mergingNEW2" "CK_2016-06-29_01_CD8_merging2")
file_metadata=("metadata_23_01.xlsx" "metadata_29_01.xlsx")

prefix_data=("23CD8_" "29CD8_")
prefix_panel="01CD8_"
prefix_pca="pca1_"
prefix_merging=("merging3_" "merging3_")

prefix_freq_merging="23CD8m3_29CD8m3_"


RWD=$RWD_MAIN/${new_data_dir}
ROUT=$RWD/Rout
mkdir -p $ROUT
echo "$RWD"


# ### Get cluster frequencies
# if ${frequencies_merged}; then
#   echo "08_frequencies_merged"
#   R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${prefix_freq_merging}' freq_outdir='08_frequencies_merged' path_metadata=c('${METADATA}/${file_metadata[0]}','${METADATA}/${file_metadata[1]}')  path_counts=c('${RWD_MAIN}/${data_dir[0]}/050_frequencies/${prefix_data[0]}${prefix_panel}${prefix_pca}${prefix_merging[0]}counts.xls','${RWD_MAIN}/${data_dir[1]}/050_frequencies/${prefix_data[1]}${prefix_panel}${prefix_pca}${prefix_merging[1]}counts.xls') data_name=c('${data_name[0]}','${data_name[1]}') path_fun_models='$RCODE/00_models_merged.R'" $RCODE/08_frequencies_merged.R $ROUT/08_frequencies_merged.Rout
#   tail $ROUT/08_frequencies_merged.Rout
# fi













#
