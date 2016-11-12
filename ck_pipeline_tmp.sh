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
cytokines_merging_v2=true # based on clustering the two joint bimatrices


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
prefix_cytokines_cutoffs1='cytCM_'
prefix_cytokines_cutoffs2='cytCM_'


prefix_merging_cyt1='cytmerging3_'
prefix_merging_cyt2='cytmerging2_'

prefix_data_merging="23CD8m2cm3_29CD8m3cm2_"

./Analysis_block_8_cytokine_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --METADATA ${METADATA} --new_data_dir ${new_data_dir} --cytokines_merging ${cytokines_merging} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --prefix_clsubset ${prefix_clsubset} --prefix_cytokines_cutoffs1 ${prefix_cytokines_cutoffs1} --prefix_cytokines_cutoffs2 ${prefix_cytokines_cutoffs2} --prefix_merging_cyt1 ${prefix_merging_cyt1} --prefix_merging_cyt2 ${prefix_merging_cyt2}

















#
