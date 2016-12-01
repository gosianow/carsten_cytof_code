#!/bin/bash

###############################################################################################################
## Define paths to software and reference files

RWD_MAIN=/Users/gosia/Dropbox/UZH/carsten_cytof
RCODE=/Users/gosia/Dropbox/UZH/carsten_cytof_code
METADATA=$RWD_MAIN/CK_metadata
PANELS=$RWD_MAIN/CK_panels

## Define which analysis to re-run
frequencies_merged=true
expression_merged=true
runtsne_merged=false
plottsne_merged=false
cytokines_merging=false
cytokines_merging_v2=true # based on clustering the two joint bimatrices
cytokines_merging_v3=true # based on selecting the most occuring combinations
cytokines_merging_v4=true
pd1_merging_v2=true
pd1_merging_v3=true
cd69_merging_v2=true
cd69_merging_v3=true


###############################################################################################################
# Merging data 23 and 29 for panel 2 CD4 --- cytokines --- cytokines_merging_v2 / cytokines_merging_v3
###############################################################################################################

# -----------------------------
# Tmem cluster
# -----------------------------

PANEL=2

new_data_dir="CK_2016-06-merged_23_29/02_CD4"

data_name1="data23"
data_name2="data29"
data_dir1="CK_2016-06-23_02_CD4_merging2/060_cytokines_bimatrix"
data_dir2="CK_2016-06-29_02_CD4_merging/060_cytokines_bimatrix"
file_metadata1="metadata_23_02.xlsx"
file_metadata2="metadata_29_02.xlsx"

prefix_data1="23CD4_"
prefix_data2="29CD4_"
prefix_panel="02CD4_"
prefix_pca="pca1_"
prefix_merging1="merging2_"
prefix_merging2="merging3_"

prefix_clsubset='Tmem_'
prefix_cytokines_cutoffs='cytCM_raw2_'

prefix_data_merging="23_29_CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_"


### cytokines_merging_v4

./Analysis_block_10_cytokine_merging_2datasets_overall_freqs.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --METADATA ${METADATA} --new_data_dir ${new_data_dir} --cytokines_merging_v4 ${cytokines_merging_v4} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --prefix_clsubset ${prefix_clsubset} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs} --outdir "10_cytokines_merged_overall_freqs"
















#
