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

prefix_freq_merging="23m6_29m4_"


./Analysis_block_8_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --new_data_dir ${new_data_dir} --frequencies_merged ${frequencies_merged} --expression_merged ${expression_merged} --runtsne_merged ${runtsne_merged} --plottsne_merged ${plottsne_merged} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_freq_merging ${prefix_freq_merging} --METADATA ${METADATA}



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

prefix_freq_merging="23m4_29m2_"


./Analysis_block_8_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --new_data_dir ${new_data_dir} --frequencies_merged ${frequencies_merged} --expression_merged ${expression_merged} --runtsne_merged ${runtsne_merged} --plottsne_merged ${plottsne_merged} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_freq_merging ${prefix_freq_merging} --METADATA ${METADATA}



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

prefix_freq_merging="23CD4m4_29CD4m4_"


./Analysis_block_8_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --new_data_dir ${new_data_dir} --frequencies_merged ${frequencies_merged} --expression_merged ${expression_merged} --runtsne_merged ${runtsne_merged} --plottsne_merged ${plottsne_merged} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_freq_merging ${prefix_freq_merging} --METADATA ${METADATA}


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

prefix_freq_merging="23CD8m3_29CD8m3_"


./Analysis_block_8_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --new_data_dir ${new_data_dir} --frequencies_merged ${frequencies_merged} --expression_merged ${expression_merged} --runtsne_merged ${runtsne_merged} --plottsne_merged ${plottsne_merged} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_freq_merging ${prefix_freq_merging} --METADATA ${METADATA}









#
