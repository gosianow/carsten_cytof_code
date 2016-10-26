#!/bin/bash

###############################################################################################################
## Define paths to software and reference files

RWD_MAIN=/Users/gosia/Dropbox/UZH/carsten_cytof
RCODE=/Users/gosia/Dropbox/UZH/carsten_cytof_code
METADATA=$RWD_MAIN/CK_metadata
PANELS=$RWD_MAIN/CK_panels

## Define which analysis to re-run
frequencies_merged=false
expression_merged=false
runtsne_merged=false
plottsne_merged=false
cytokines_merging=true

###############################################################################################################
# Merging data 23 and 29 for panel 1
###############################################################################################################

PANEL=1

new_data_dir="CK_2016-06-merged_23_29/01"

data_name1="data23"
data_name2="data29"
data_dir1="CK_2016-06-23_01"
data_dir2="CK_2016-06-29_01"
file_metadata1="metadata_23_01.xlsx"
file_metadata2="metadata_29_01.xlsx"

prefix_data1="23_"
prefix_data2="29_"
prefix_panel="01_"
prefix_pca="pca1_"
prefix_merging1="merging6_"
prefix_merging2="merging4_"

prefix_data_merging="23m6_29m4_"


./Analysis_block_8_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --new_data_dir ${new_data_dir} --frequencies_merged ${frequencies_merged} --expression_merged ${expression_merged} --runtsne_merged ${runtsne_merged} --plottsne_merged ${plottsne_merged} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --METADATA ${METADATA}


###############################################################################################################
# Merging data 23 and 29 for panel 3
###############################################################################################################

PANEL=3

new_data_dir="CK_2016-06-merged_23_29/03"

data_name1="data23"
data_name2="data29"
data_dir1="CK_2016-06-23_03"
data_dir2="CK_2016-06-29_03"
file_metadata1="metadata_23_03.xlsx"
file_metadata2="metadata_29_03.xlsx"

prefix_data1="23_"
prefix_data2="29_"
prefix_panel="03_"
prefix_pca="pca1_"
prefix_merging1="merging4_"
prefix_merging2="merging2_"

prefix_data_merging="23m4_29m2_"


./Analysis_block_8_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --new_data_dir ${new_data_dir} --frequencies_merged ${frequencies_merged} --expression_merged ${expression_merged} --runtsne_merged ${runtsne_merged} --plottsne_merged ${plottsne_merged} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --METADATA ${METADATA}



###############################################################################################################
# Merging data 23 and 29 for panel 1 CD4
###############################################################################################################

PANEL=1

new_data_dir="CK_2016-06-merged_23_29/01_CD4"

data_name1="data23"
data_name2="data29"
data_dir1="CK_2016-06-23_01_CD4_mergingNEW2"
data_dir2="CK_2016-06-29_01_CD4_merging2"
file_metadata1="metadata_23_01.xlsx"
file_metadata2="metadata_29_01.xlsx"

prefix_data1="23CD4_"
prefix_data2="29CD4_"
prefix_panel="01CD4_"
prefix_pca="pca1_"
prefix_merging1="merging4_"
prefix_merging2="merging4_"

prefix_data_merging="23CD4m4_29CD4m4_"


./Analysis_block_8_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --new_data_dir ${new_data_dir} --frequencies_merged ${frequencies_merged} --expression_merged ${expression_merged} --runtsne_merged ${runtsne_merged} --plottsne_merged ${plottsne_merged} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --METADATA ${METADATA}


###############################################################################################################
# Merging data 23 and 29 for panel 1 CD8
###############################################################################################################

PANEL=1

new_data_dir="CK_2016-06-merged_23_29/01_CD8"

data_name1="data23"
data_name2="data29"
data_dir1="CK_2016-06-23_01_CD8_mergingNEW2"
data_dir2="CK_2016-06-29_01_CD8_merging2"
file_metadata1="metadata_23_01.xlsx"
file_metadata2="metadata_29_01.xlsx"

prefix_data1="23CD8_"
prefix_data2="29CD8_"
prefix_panel="01CD8_"
prefix_pca="pca1_"
prefix_merging1="merging3_"
prefix_merging2="merging3_"

prefix_data_merging="23CD8m3_29CD8m3_"


./Analysis_block_8_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --new_data_dir ${new_data_dir} --frequencies_merged ${frequencies_merged} --expression_merged ${expression_merged} --runtsne_merged ${runtsne_merged} --plottsne_merged ${plottsne_merged} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --METADATA ${METADATA}

###############################################################################################################
# Merging data 23 and 29 for panel 2 CD4 - cytokines
###############################################################################################################

PANEL=2

new_data_dir="CK_2016-06-merged_23_29/02_CD4"

data_name1="data23"
data_name2="data29"
data_dir1="CK_2016-06-23_02_CD4_merging2"
data_dir2="CK_2016-06-29_02_CD4_merging"
file_metadata1="metadata_23_02.xlsx"
file_metadata2="metadata_29_02.xlsx"

prefix_data1="23CD4_"
prefix_data2="29CD4_"
prefix_panel="02CD4_"
prefix_pca="pca1_"
prefix_merging1="merging2_"
prefix_merging2="merging3_"

prefix_clsubset='Tmem_'
prefix_cytokines_cutoffs='cytCM_'
prefix_merging_cyt1='cytmerging3_'
prefix_merging_cyt2='cytmerging2_'

prefix_data_merging="23CD4m2cm3_29CD4m3cm2_"

RWD=$RWD_MAIN/${new_data_dir}
ROUT=$RWD/Rout
mkdir -p $ROUT

### Get cluster frequencies
if ${cytokines_merging}; then
  echo "08_frequencies_merged"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${prefix_data_merging}' freq_outdir='08_cytokines_merged' path_metadata=c('${METADATA}/${file_metadata1}','${METADATA}/${file_metadata2}')  path_counts=c('${RWD_MAIN}/${data_dir1}/060_cytokines_bimatrix/${prefix_data1}${prefix_panel}${prefix_pca}${prefix_merging1}${prefix_clsubset}${prefix_cytokines_cutoffs}raw2_${prefix_merging_cyt1}counts.xls','${RWD_MAIN}/${data_dir2}/060_cytokines_bimatrix/${prefix_data2}${prefix_panel}${prefix_pca}${prefix_merging2}${prefix_clsubset}${prefix_cytokines_cutoffs}raw2_${prefix_merging_cyt2}counts.xls') data_name=c('${data_name1}','${data_name2}') path_fun_models='$RCODE/00_models.R'" $RCODE/08_frequencies_merged.R $ROUT/08_frequencies_merged.Rout
  tail $ROUT/08_frequencies_merged.Rout
fi


###############################################################################################################
# Merging data 23 and 29 for panel 2 CD8 - cytokines
###############################################################################################################

PANEL=2

new_data_dir="CK_2016-06-merged_23_29/02_CD8"

data_name1="data23"
data_name2="data29"
data_dir1="CK_2016-06-23_02_CD8_merging2"
data_dir2="CK_2016-06-29_02_CD8_merging"
file_metadata1="metadata_23_02.xlsx"
file_metadata2="metadata_29_02.xlsx"

prefix_data1="23CD8_"
prefix_data2="29CD8_"
prefix_panel="02CD8_"
prefix_pca="pca1_"
prefix_merging1="merging2_"
prefix_merging2="merging3_"

prefix_clsubset='Tmem_'
prefix_cytokines_cutoffs='cytCM_'
prefix_merging_cyt1='cytmerging3_'
prefix_merging_cyt2='cytmerging2_'

prefix_data_merging="23CD8m2cm3_29CD8m3cm2_"

RWD=$RWD_MAIN/${new_data_dir}
ROUT=$RWD/Rout
mkdir -p $ROUT

### Get cluster frequencies
if ${cytokines_merging}; then
  echo "08_frequencies_merged"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${prefix_data_merging}' freq_outdir='08_cytokines_merged' path_metadata=c('${METADATA}/${file_metadata1}','${METADATA}/${file_metadata2}')  path_counts=c('${RWD_MAIN}/${data_dir1}/060_cytokines_bimatrix/${prefix_data1}${prefix_panel}${prefix_pca}${prefix_merging1}${prefix_clsubset}${prefix_cytokines_cutoffs}raw2_${prefix_merging_cyt1}counts.xls','${RWD_MAIN}/${data_dir2}/060_cytokines_bimatrix/${prefix_data2}${prefix_panel}${prefix_pca}${prefix_merging2}${prefix_clsubset}${prefix_cytokines_cutoffs}raw2_${prefix_merging_cyt2}counts.xls') data_name=c('${data_name1}','${data_name2}') path_fun_models='$RCODE/00_models.R'" $RCODE/08_frequencies_merged.R $ROUT/08_frequencies_merged.Rout
  tail $ROUT/08_frequencies_merged.Rout
fi























#
