#!/bin/bash

###############################################################################################################
## Define paths to software and reference files

RWD_MAIN=/Users/gosia/Dropbox/UZH/carsten_cytof/PD1_project
RCODE=/Users/gosia/Dropbox/UZH/carsten_cytof_code
# RWD_MAIN=/home/Shared/data/cytof/carsten_cytof
# RCODE=/home/gosia/R/carsten_cytof_code
METADATA=$RWD_MAIN/CK_metadata
PANELS=$RWD_MAIN/CK_panels

## Define which analysis to re-run
frequencies_merged=false
expression_merged=false

cd69_merging_v2=false
cd69_merging_v3=false


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


./Analysis_block_8_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --new_data_dir ${new_data_dir} --frequencies_merged ${frequencies_merged} --expression_merged ${expression_merged} --runtsne_merged false --plottsne_merged false --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --METADATA ${METADATA}


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


./Analysis_block_8_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --new_data_dir ${new_data_dir} --frequencies_merged ${frequencies_merged} --expression_merged ${expression_merged} --runtsne_merged false --plottsne_merged false --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --METADATA ${METADATA}


###############################################################################################################
# Merging data 23all and 29all3 for panel 3v2
###############################################################################################################

new_data_dir="CK_2016-06-merged_23_29/03all"

data_name1="data23"
data_name2="data29"
data_dir1="CK_2016-06-23_03all"
data_dir2="CK_2016-06-29_03all3"
file_metadata1="metadata_23_03all.xlsx"
file_metadata2="metadata_29_03all3.xlsx"

prefix_data1="23all_"
prefix_data2="29all3_"
prefix_panel="03v2_"
prefix_pca="pca1_"
prefix_merging1="merging4_"
prefix_merging2="merging_"

prefix_data_merging="23allm4_29all3m_"


./Analysis_block_8_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --new_data_dir ${new_data_dir} --frequencies_merged ${frequencies_merged} --expression_merged ${expression_merged} --runtsne_merged false --plottsne_merged false --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --METADATA ${METADATA}


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


# ----------------------------------------------------
# 23 merging 5 + 29 merging 5
# ----------------------------------------------------

prefix_merging1="merging5_"
prefix_merging2="merging5_"

prefix_data_merging="23CD4m5_29CD4m5_"


./Analysis_block_8_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --new_data_dir ${new_data_dir} --frequencies_merged ${frequencies_merged} --expression_merged ${expression_merged} --runtsne_merged false --plottsne_merged false --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --METADATA ${METADATA}




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

# ----------------------------------------------------
# 23 merging 6 + 29 merging 6
# ----------------------------------------------------

prefix_merging1="merging6_"
prefix_merging2="merging6_"

prefix_data_merging="23CD8m6_29CD8m6_"


./Analysis_block_8_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --new_data_dir ${new_data_dir} --frequencies_merged ${frequencies_merged} --expression_merged ${expression_merged} --runtsne_merged false --plottsne_merged false --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --METADATA ${METADATA}



###############################################################################################################
# Merging data 23 and 29 for panel 1 CD4 --- no HD10
###############################################################################################################

PANEL=1

new_data_dir="CK_2016-06-merged_23_29/01_CD4"

data_name1="data23"
data_name2="data29"
data_dir1="CK_2016-06-23_01_CD4_mergingNEW2"
data_dir2="CK_2016-06-29_01_CD4_merging2"
file_metadata1="metadata_23_01.xlsx"
file_metadata2="metadata_29_01_noHD10.xlsx"

prefix_data1="23CD4_"
prefix_data2="29CD4_"
prefix_panel="01CD4_"
prefix_pca="pca1_"

# ----------------------------------------------------
# 23 merging 5 + 29 merging 5
# ----------------------------------------------------

prefix_merging1="merging5_"
prefix_merging2="merging5_"

prefix_data_merging="noHD10_23CD4m5_29CD4m5_"


./Analysis_block_8_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --new_data_dir ${new_data_dir} --frequencies_merged ${frequencies_merged} --expression_merged ${expression_merged} --runtsne_merged false --plottsne_merged false --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --METADATA ${METADATA}



###############################################################################################################
# Merging data 23 and 29 for panel 1 CD8 --- no HD10
###############################################################################################################

PANEL=1

new_data_dir="CK_2016-06-merged_23_29/01_CD8"

data_name1="data23"
data_name2="data29"
data_dir1="CK_2016-06-23_01_CD8_mergingNEW2"
data_dir2="CK_2016-06-29_01_CD8_merging2"
file_metadata1="metadata_23_01.xlsx"
file_metadata2="metadata_29_01_noHD10.xlsx"

prefix_data1="23CD8_"
prefix_data2="29CD8_"
prefix_panel="01CD8_"
prefix_pca="pca1_"

# ----------------------------------------------------
# 23 merging 6 + 29 merging 6
# ----------------------------------------------------

prefix_merging1="merging6_"
prefix_merging2="merging6_"

prefix_data_merging="noHD10_23CD8m6_29CD8m6_"


./Analysis_block_8_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --new_data_dir ${new_data_dir} --frequencies_merged ${frequencies_merged} --expression_merged ${expression_merged} --runtsne_merged false --plottsne_merged false --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --METADATA ${METADATA}




###############################################################################################################
# Merging data 23 and 29 for panel 2 CD4 --- cd69 --- cytokines_merging_v2 / cytokines_merging_v3
###############################################################################################################

# -----------------------------
# Tmem cluster
# -----------------------------

PANEL=2

new_data_dir="CK_2016-06-merged_23_29/02_CD4"

data_name1="data23"
data_name2="data29"
data_dir1="CK_2016-06-23_02_CD4_merging2/070_cd69_bimatrix"
data_dir2="CK_2016-06-29_02_CD4_merging/070_cd69_bimatrix"
file_metadata1="metadata_23_02.xlsx"
file_metadata2="metadata_29_02.xlsx"

prefix_data1="23CD4_"
prefix_data2="29CD4_"
prefix_panel="02CD4_"
prefix_pca="pca1_"
prefix_merging1="merging2_"
prefix_merging2="merging3_"

prefix_clsubset='Tmem_'
prefix_cytokines_cutoffs='cytCM_raw2_cd69_positive_'

prefix_data_merging="23_29_CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_cd69_positive_"


### cytokines_merging_v4

./Analysis_block_10_cytokine_merging_2datasets_overall_freqs.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --METADATA ${METADATA} --new_data_dir ${new_data_dir} --cytokines_merging_v4 ${cytokines_merging_v4} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --prefix_clsubset ${prefix_clsubset} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs} --outdir "10_cd69_merged_overall_freqs"


for i in 7
do
  som_dim=$i
  nmetaclusts=$((som_dim*som_dim))
  prefix_clust="cl${nmetaclusts}_"
  ./Analysis_block_10_cytokine_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --METADATA ${METADATA} --new_data_dir ${new_data_dir} --cytokines_merging_v2 ${cd69_merging_v2} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --prefix_clsubset ${prefix_clsubset} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs} --prefix_clust ${prefix_clust} --som_dim ${som_dim} --nmetaclusts ${nmetaclusts} --outdir "10_cd69_merged"
done


for i in 40
do
  top_combinations=$i
  prefix_clust="cl${i}_"
  ./Analysis_block_10_cytokine_merging_2datasets_top_selection.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --METADATA ${METADATA} --new_data_dir ${new_data_dir} --cytokines_merging_v3 ${cd69_merging_v3} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --prefix_clsubset ${prefix_clsubset} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs} --prefix_clust ${prefix_clust} --top_combinations ${top_combinations} --outdir "10_cd69_merged_top_combinations"
done


###############################################################################################################
# Merging data 23 and 29 for panel 2 CD8 --- cd69 --- cytokines_merging_v2 / cytokines_merging_v3
###############################################################################################################

# -----------------------------
# Tmem cluster
# -----------------------------

PANEL=2

new_data_dir="CK_2016-06-merged_23_29/02_CD8"

data_name1="data23"
data_name2="data29"
data_dir1="CK_2016-06-23_02_CD8_merging2/070_cd69_bimatrix"
data_dir2="CK_2016-06-29_02_CD8_merging/070_cd69_bimatrix"
file_metadata1="metadata_23_02.xlsx"
file_metadata2="metadata_29_02.xlsx"

prefix_data1="23CD8_"
prefix_data2="29CD8_"
prefix_panel="02CD8_"
prefix_pca="pca1_"
prefix_merging1="merging2_"
prefix_merging2="merging3_"

prefix_clsubset='Tmem_'
prefix_cytokines_cutoffs='cytCM_raw2_cd69_positive_'

prefix_data_merging="23_29_CD8_02CD8_pca1_merging2_Tmem_cytCM_raw2_cd69_positive_"


### cytokines_merging_v4

./Analysis_block_10_cytokine_merging_2datasets_overall_freqs.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --METADATA ${METADATA} --new_data_dir ${new_data_dir} --cytokines_merging_v4 ${cytokines_merging_v4} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --prefix_clsubset ${prefix_clsubset} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs} --outdir "10_cd69_merged_overall_freqs"



for i in 7
do
  som_dim=$i
  nmetaclusts=$((som_dim*som_dim))
  prefix_clust="cl${nmetaclusts}_"
  ./Analysis_block_10_cytokine_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --METADATA ${METADATA} --new_data_dir ${new_data_dir} --cytokines_merging_v2 ${cd69_merging_v2} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --prefix_clsubset ${prefix_clsubset} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs} --prefix_clust ${prefix_clust} --som_dim ${som_dim} --nmetaclusts ${nmetaclusts} --outdir "10_cd69_merged"
done


for i in 40
do
  top_combinations=$i
  prefix_clust="cl${i}_"
  ./Analysis_block_10_cytokine_merging_2datasets_top_selection.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --METADATA ${METADATA} --new_data_dir ${new_data_dir} --cytokines_merging_v3 ${cd69_merging_v3} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --prefix_clsubset ${prefix_clsubset} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs} --prefix_clust ${prefix_clust} --top_combinations ${top_combinations} --outdir "10_cd69_merged_top_combinations"
done









#
