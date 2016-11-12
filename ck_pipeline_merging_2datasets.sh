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
cytokines_merging=false
cytokines_merging_v2=false # based on clustering the two joint bimatrices

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

# ----------------------------------------------------
# 23 merging 5 + 29 merging 5
# ----------------------------------------------------

prefix_merging1="merging5_"
prefix_merging2="merging5_"

prefix_data_merging="23CD4m5_29CD4m5_"


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

# ----------------------------------------------------
# 23 merging 4 + 29 merging 4
# ----------------------------------------------------

prefix_merging1="merging4_"
prefix_merging2="merging4_"

prefix_data_merging="23CD8m4_29CD8m4_"


./Analysis_block_8_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --new_data_dir ${new_data_dir} --frequencies_merged ${frequencies_merged} --expression_merged ${expression_merged} --runtsne_merged ${runtsne_merged} --plottsne_merged ${plottsne_merged} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --METADATA ${METADATA}


# ----------------------------------------------------
# 23 merging 5 + 29 merging 4
# ----------------------------------------------------

prefix_merging1="merging5_"
prefix_merging2="merging4_"

prefix_data_merging="23CD8m5_29CD8m4_"


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

prefix_cytokines_cutoffs1='cytCMbak_'
prefix_cytokines_cutoffs2='cytCM_'

prefix_merging_cyt1='cytmerging3_'
prefix_merging_cyt2='cytmerging2_'

prefix_data_merging="23CD4m2cm3_29CD4m3cm2_"


./Analysis_block_8_cytokine_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --METADATA ${METADATA} --new_data_dir ${new_data_dir} --cytokines_merging ${cytokines_merging} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --prefix_clsubset ${prefix_clsubset} --prefix_cytokines_cutoffs1 ${prefix_cytokines_cutoffs1} --prefix_cytokines_cutoffs2 ${prefix_cytokines_cutoffs2} --prefix_merging_cyt1 ${prefix_merging_cyt1} --prefix_merging_cyt2 ${prefix_merging_cyt2}


prefix_cytokines_cutoffs1='cytCM_'
prefix_cytokines_cutoffs2='cytCM_'

prefix_merging_cyt1='cytmerging4_'
prefix_merging_cyt2='cytmerging4_'

prefix_data_merging="23CD4m2cm4_29CD4m3cm4_"

./Analysis_block_8_cytokine_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --METADATA ${METADATA} --new_data_dir ${new_data_dir} --cytokines_merging ${cytokines_merging} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --prefix_clsubset ${prefix_clsubset} --prefix_cytokines_cutoffs1 ${prefix_cytokines_cutoffs1} --prefix_cytokines_cutoffs2 ${prefix_cytokines_cutoffs2} --prefix_merging_cyt1 ${prefix_merging_cyt1} --prefix_merging_cyt2 ${prefix_merging_cyt2}



prefix_merging_cyt1='cytmerging5_'
prefix_merging_cyt2='cytmerging5_'

prefix_data_merging="23CD4m2cm5_29CD4m3cm5_"

./Analysis_block_8_cytokine_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --METADATA ${METADATA} --new_data_dir ${new_data_dir} --cytokines_merging ${cytokines_merging} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --prefix_clsubset ${prefix_clsubset} --prefix_cytokines_cutoffs1 ${prefix_cytokines_cutoffs1} --prefix_cytokines_cutoffs2 ${prefix_cytokines_cutoffs2} --prefix_merging_cyt1 ${prefix_merging_cyt1} --prefix_merging_cyt2 ${prefix_merging_cyt2}



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



prefix_merging_cyt1='cytmerging4_'
prefix_merging_cyt2='cytmerging4_'

prefix_data_merging="23CD8m2cm4_29CD8m3cm4_"

./Analysis_block_8_cytokine_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --METADATA ${METADATA} --new_data_dir ${new_data_dir} --cytokines_merging ${cytokines_merging} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --prefix_clsubset ${prefix_clsubset} --prefix_cytokines_cutoffs1 ${prefix_cytokines_cutoffs1} --prefix_cytokines_cutoffs2 ${prefix_cytokines_cutoffs2} --prefix_merging_cyt1 ${prefix_merging_cyt1} --prefix_merging_cyt2 ${prefix_merging_cyt2}




###############################################################################################################
# Merging data 23 and 29 for panel 2 CD4 - cytokines - cytokines_merging_v2
###############################################################################################################

# -----------------------------
# Tmem cluster
# -----------------------------

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

prefix_data_merging="23_29_CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_"

prefix_clust='cl25_'
som_dim=5

./Analysis_block_10_cytokine_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --METADATA ${METADATA} --new_data_dir ${new_data_dir} --cytokines_merging_v2 ${cytokines_merging_v2} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --prefix_clsubset ${prefix_clsubset} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs} --prefix_clust ${prefix_clust} --som_dim ${som_dim}


prefix_clust='cl36_'
som_dim=6

./Analysis_block_10_cytokine_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --METADATA ${METADATA} --new_data_dir ${new_data_dir} --cytokines_merging_v2 ${cytokines_merging_v2} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --prefix_clsubset ${prefix_clsubset} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs} --prefix_clust ${prefix_clust} --som_dim ${som_dim}


prefix_clust='cl49_'
som_dim=7

./Analysis_block_10_cytokine_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --METADATA ${METADATA} --new_data_dir ${new_data_dir} --cytokines_merging_v2 ${cytokines_merging_v2} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --prefix_clsubset ${prefix_clsubset} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs} --prefix_clust ${prefix_clust} --som_dim ${som_dim}


# -----------------------------
# EM+CM cluster
# -----------------------------

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

prefix_clsubset='EM_CM_'
prefix_cytokines_cutoffs='cytCM_'

prefix_data_merging="23_29_CD4_02CD4_pca1_merging2_EM_CM_cytCM_raw2_"

prefix_clust='cl25_'
som_dim=5

./Analysis_block_10_cytokine_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --METADATA ${METADATA} --new_data_dir ${new_data_dir} --cytokines_merging_v2 ${cytokines_merging_v2} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --prefix_clsubset ${prefix_clsubset} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs} --prefix_clust ${prefix_clust} --som_dim ${som_dim}


prefix_clust='cl36_'
som_dim=6

./Analysis_block_10_cytokine_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --METADATA ${METADATA} --new_data_dir ${new_data_dir} --cytokines_merging_v2 ${cytokines_merging_v2} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --prefix_clsubset ${prefix_clsubset} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs} --prefix_clust ${prefix_clust} --som_dim ${som_dim}


###############################################################################################################
# Merging data 23 and 29 for panel 2 CD8 - cytokines - cytokines_merging_v2
###############################################################################################################

# -----------------------------
# Tmem cluster
# -----------------------------

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

prefix_data_merging="23_29_CD8_02CD8_pca1_merging2_Tmem_cytCM_raw2_"

prefix_clust='cl25_'
som_dim=5

./Analysis_block_10_cytokine_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --METADATA ${METADATA} --new_data_dir ${new_data_dir} --cytokines_merging_v2 ${cytokines_merging_v2} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --prefix_clsubset ${prefix_clsubset} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs} --prefix_clust ${prefix_clust} --som_dim ${som_dim}


prefix_clust='cl36_'
som_dim=6

./Analysis_block_10_cytokine_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --METADATA ${METADATA} --new_data_dir ${new_data_dir} --cytokines_merging_v2 ${cytokines_merging_v2} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --prefix_clsubset ${prefix_clsubset} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs} --prefix_clust ${prefix_clust} --som_dim ${som_dim}


prefix_clust='cl49_'
som_dim=7

./Analysis_block_10_cytokine_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --METADATA ${METADATA} --new_data_dir ${new_data_dir} --cytokines_merging_v2 ${cytokines_merging_v2} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --prefix_clsubset ${prefix_clsubset} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs} --prefix_clust ${prefix_clust} --som_dim ${som_dim}


# -----------------------------
# CM cluster
# -----------------------------

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

prefix_clsubset='CM_'
prefix_cytokines_cutoffs='cytCM_'

prefix_data_merging="23_29_CD8_02CD8_pca1_merging2_CM_cytCM_raw2_"


prefix_clust='cl25_'
som_dim=5

./Analysis_block_10_cytokine_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --METADATA ${METADATA} --new_data_dir ${new_data_dir} --cytokines_merging_v2 ${cytokines_merging_v2} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --prefix_clsubset ${prefix_clsubset} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs} --prefix_clust ${prefix_clust} --som_dim ${som_dim}


prefix_clust='cl36_'
som_dim=6

./Analysis_block_10_cytokine_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --METADATA ${METADATA} --new_data_dir ${new_data_dir} --cytokines_merging_v2 ${cytokines_merging_v2} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --prefix_clsubset ${prefix_clsubset} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs} --prefix_clust ${prefix_clust} --som_dim ${som_dim}


# -----------------------------
# CM+naive cluster
# -----------------------------

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

prefix_clsubset='CM_naive_'
prefix_cytokines_cutoffs='cytCM_'

prefix_data_merging="23_29_CD8_02CD8_pca1_merging2_CM_naive_cytCM_raw2_"


prefix_clust='cl25_'
som_dim=5

./Analysis_block_10_cytokine_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --METADATA ${METADATA} --new_data_dir ${new_data_dir} --cytokines_merging_v2 ${cytokines_merging_v2} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --prefix_clsubset ${prefix_clsubset} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs} --prefix_clust ${prefix_clust} --som_dim ${som_dim}


prefix_clust='cl36_'
som_dim=6

./Analysis_block_10_cytokine_merging_2datasets.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --METADATA ${METADATA} --new_data_dir ${new_data_dir} --cytokines_merging_v2 ${cytokines_merging_v2} --data_name1 ${data_name1} --data_name2 ${data_name2} --data_dir1 ${data_dir1} --data_dir2 ${data_dir2} --file_metadata1 ${file_metadata1} --file_metadata2 ${file_metadata2} --prefix_data1 ${prefix_data1} --prefix_data2 ${prefix_data2} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging1 ${prefix_merging1} --prefix_merging2 ${prefix_merging2} --prefix_data_merging ${prefix_data_merging} --prefix_clsubset ${prefix_clsubset} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs} --prefix_clust ${prefix_clust} --som_dim ${som_dim}






















#
