#!/bin/bash

###############################################################################################################
## Define paths to software and reference files

RWD_MAIN=/Users/gosia/Dropbox/UZH/carsten_cytof
RCODE=/Users/gosia/Dropbox/UZH/carsten_cytof_code
METADATA=$RWD_MAIN/CK_metadata
PANELS=$RWD_MAIN/CK_panels


## Define which analysis to re-run
data_normalization=false
pcascores=false
select_observables=false
flowsom=false
heatmaps=false
runtsne=false
plottsne=false
frequencies=false
expression=true
cluster_merging=false
cluster_extracting=false
fcs_saving=false
cytokines=false
pd1=false

###############################################################################################################
# Analysis of CK_2016-06-23_01 data
# Use Analysis block 1
###############################################################################################################

DATA=23
PANEL=1
data_dir="CK_2016-06-23_01"

file_panel="panel1.xlsx"
file_metadata="metadata_23_01.xlsx"

pca_score_cutoff=3
rand_seed_consensus=123
nmetaclusts=20

prefix_data="23_"
prefix_panel="01_"
prefix_pca="pca1_"
prefix_clust="cl20_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts}

rand_seed_consensus=123
nmetaclusts=9
prefix_clust="cl9_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts}


# --------------------------------------------------
# Analysis of CK_2016-06-23_01 cluster_merging
# Use Analysis block 2
# --------------------------------------------------

file_merging="${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_mergingNEW.xlsx"
prefix_merging="mergingNEW_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging}


file_merging="${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_mergingNEW2.xlsx"
prefix_merging="mergingNEW2_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging}


file_merging="${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_merging5.xlsx"
prefix_merging="merging5_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging}


# --------------------------------------------------
# CK_2016-06-23_01 - CD4 and CD8 cluster extracting from mergingNEW2
# Use Analysis block 3
# --------------------------------------------------

prefix_merging="mergingNEW2_"

extract_cluster=('CD4' 'CD8')
extract_dir=('CK_2016-06-23_01_CD4_mergingNEW2' 'CK_2016-06-23_01_CD8_mergingNEW2')

for i in 0 1
do
  ./Analysis_block_3_cluster_extracting.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_extracting ${cluster_extracting} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging} --extract_cluster ${extract_cluster[$i]} --extract_dir ${extract_dir[$i]}
done


# --------------------------------------------------
# Analysis of CK_2016-06-23_01_CD4_mergingNEW2 using panel1CD4.xlsx
# and CK_2016-06-23_01_CD8_mergingNEW2 using panel1CD8.xlsx
# Use Analysis block 1
# --------------------------------------------------

DATA=23
PANEL=1

file_metadata="metadata_23_01.xlsx"

rand_seed_consensus=1234
nmetaclusts=20

prefix_pca="pca1_"
prefix_clust="cl20_"

data_dir=('CK_2016-06-23_01_CD4_mergingNEW2' 'CK_2016-06-23_01_CD8_mergingNEW2')
prefix_data=('23CD4_' '23CD8_')
file_panel=('panel1CD4.xlsx' 'panel1CD8.xlsx')
prefix_panel=('01CD4_' '01CD8_')
pca_score_cutoff=(2 2)

for i in 0 1
do
  ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff[$i]} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts}
done


rand_seed_consensus=1234
nmetaclusts=(6 7)
prefix_clust=("cl6_" "cl7_")

for i in 0 1
do
  ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust[$i]} --pca_score_cutoff ${pca_score_cutoff[$i]} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts[$i]}
done


# --------------------------------------------------
# Analysis of CK_2016-06-23_01_CD4_mergingNEW2 using panel1CD4.xlsx for cluster_merging
# and CK_2016-06-23_01_CD8_mergingNEW2 using panel1CD8.xlsx for cluster_merging
# Use Analysis block 2
# --------------------------------------------------

file_merging=("${prefix_data[0]}${prefix_panel[0]}${prefix_pca}${prefix_clust}cluster_merging.xlsx" "${prefix_data[1]}${prefix_panel[1]}${prefix_pca}${prefix_clust}cluster_merging.xlsx")
prefix_merging=('merging_' 'merging_')

for i in 0 1
do
  ./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging[$i]} --file_merging ${file_merging[$i]}
done


###############################################################################################################
# Analysis of CK_2016-06-23_02 data
# Use Analysis block 1
###############################################################################################################

DATA=23
PANEL=2
data_dir="CK_2016-06-23_02"

file_panel="panel2.xlsx"
file_metadata="metadata_23_02.xlsx"

pca_score_cutoff=1
rand_seed_consensus=1234
nmetaclusts=20

prefix_data="23_"
prefix_panel="02_"
prefix_pca="pca1_"
prefix_clust="cl20_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts}


rand_seed_consensus=1234
nmetaclusts=9
prefix_clust="cl9_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts}


# --------------------------------------------------
# Analysis of CK_2016-06-23_02 merging2
# Use Analysis block 2
# --------------------------------------------------

file_merging="${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_merging2.xlsx"
prefix_merging="merging2_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging}

# --------------------------------------------------
# CK_2016-06-23_02 - CD4 and CD8 cluster extracting from merging2
# Use Analysis block 3
# --------------------------------------------------

prefix_merging="merging2_"

extract_cluster=('CD4' 'CD8')
extract_dir=('CK_2016-06-23_02_CD4_merging2' 'CK_2016-06-23_02_CD8_merging2')

for i in 0 1
do
  ./Analysis_block_3_cluster_extracting.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_extracting ${cluster_extracting} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging} --extract_cluster ${extract_cluster[$i]} --extract_dir ${extract_dir[$i]}
done

# --------------------------------------------------
# Analysis of CK_2016-06-23_02_CD4_merging2 using panel2CD4.xlsx
# and CK_2016-06-23_02_CD8_merging2 using panel2CD8.xlsx
# Use Analysis block 1
# --------------------------------------------------

DATA=23
PANEL=2
file_metadata="metadata_23_02.xlsx"

rand_seed_consensus=1234
nmetaclusts=20

prefix_pca="pca1_"
prefix_clust="cl20_"

data_dir=('CK_2016-06-23_02_CD4_merging2' 'CK_2016-06-23_02_CD8_merging2')
prefix_data=('23CD4_' '23CD8_')
file_panel=('panel2CD4.xlsx' 'panel2CD8.xlsx')
prefix_panel=('02CD4_' '02CD8_')
pca_score_cutoff=(0.96 1)

for i in 0 1
do
  ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff[$i]} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts}
done

rand_seed_consensus=1234
nmetaclusts=(6 11)
prefix_clust=("cl6_" "cl11_")

for i in 0 1
do
  ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust[$i]} --pca_score_cutoff ${pca_score_cutoff[$i]} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts[$i]}
done


# --------------------------------------------------
# Analysis of CK_2016-06-23_02_CD4_merging2 using panel2CD4.xlsx for cluster_merging
# and CK_2016-06-23_02_CD8_merging2 using panel2CD8.xlsx for cluster_merging
# Use Analysis block 2
# --------------------------------------------------

file_merging=("${prefix_data[0]}${prefix_panel[0]}${prefix_pca}${prefix_clust}cluster_merging.xlsx" "${prefix_data[1]}${prefix_panel[1]}${prefix_pca}${prefix_clust}cluster_merging.xlsx")
prefix_merging=('merging_' 'merging_')

for i in 0 1
do
  ./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging[$i]} --file_merging ${file_merging[$i]}
done


# --------------------------------------------------
# Analysis of cytokines
# in CK_2016-06-23_02_CD4_merging2 and CK_2016-06-23_02_CD8_merging2
# Use Analysis block 4
# --------------------------------------------------

prefix_merging=('merging_' 'merging_')
clsubset=("c('CM','EM')" "c('CM','EM','TE')")
prefix_clsubset=('Tmem_' 'Tmem_')
file_cytokines_cutoffs=('panel2CD4_cytokines_CM.xlsx' 'panel2CD8_cytokines_CM.xlsx')
prefix_cytokines_cutoffs=('cytCM_' 'cytCM_')
nmetaclusts=(40 20)

file_merging_cyt=("${prefix_data[0]}${prefix_panel[0]}${prefix_pca}${prefix_merging[0]}${prefix_clsubset[0]}${prefix_cytokines_cutoffs[0]}raw2_cl${nmetaclusts[0]}_cytokine_merging.xlsx" "${prefix_data[1]}${prefix_panel[1]}${prefix_pca}${prefix_merging[1]}${prefix_clsubset[1]}${prefix_cytokines_cutoffs[1]}raw2_cl${nmetaclusts[1]}_cytokine_merging.xlsx")
prefix_merging_cyt=("cytmerging_" "cytmerging_")


for i in 0
do
  ./Analysis_block_4_cytokines.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cytokines ${cytokines} --PANELS ${PANELS} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging[$i]} --prefix_clsubset ${prefix_clsubset[$i]} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs[$i]} --file_cytokines_cutoffs ${file_cytokines_cutoffs[$i]} --clsubset ${clsubset[$i]} --nmetaclusts ${nmetaclusts[$i]} --file_merging_cyt ${file_merging_cyt[$i]} --prefix_merging_cyt ${prefix_merging_cyt[$i]}
done




# --------------------------------------------------
# Analysis of PD1+ and PD1- cells
# in CK_2016-06-23_02_CD4_merging2 and CK_2016-06-23_02_CD8_merging2
# Use Analysis block 5
# --------------------------------------------------

prefix_merging=('merging_' 'merging_')
clsubset=("c('CM','EM')" "c('CM','EM','TE')")
prefix_clsubset=('Tmem_' 'Tmem_')
file_cytokines_cutoffs=('panel2CD4_cytokines_CM.xlsx' 'panel2CD8_cytokines_CM.xlsx')
nmetaclusts=(40 20)
prefix_cytokines_cutoffs=('cytCM_' 'cytCM_')

file_merging_pd1=("${prefix_data[0]}${prefix_panel[0]}${prefix_pca}${prefix_merging[0]}${prefix_clsubset[0]}${prefix_cytokines_cutoffs[0]}raw2_pd1_cl${nmetaclusts[0]}_pd1_merging.xlsx" "${prefix_data[1]}${prefix_panel[1]}${prefix_pca}${prefix_merging[1]}${prefix_clsubset[1]}${prefix_cytokines_cutoffs[1]}raw2_pd1_cl${nmetaclusts[1]}_pd1_merging.xlsx")
prefix_merging_pd1=("pd1merging_" "pd1merging_")


for i in 0 1
do
  ./Analysis_block_5_pd1.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --pd1 ${pd1} --PANELS ${PANELS} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging[$i]} --prefix_clsubset ${prefix_clsubset[$i]} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs[$i]} --file_cytokines_cutoffs ${file_cytokines_cutoffs[$i]} --clsubset ${clsubset[$i]} --nmetaclusts ${nmetaclusts[$i]} --file_merging_pd1 ${file_merging_pd1[$i]} --prefix_merging_pd1 ${prefix_merging_pd1[$i]}
done


###############################################################################################################
# Analysis of CK_2016-06-23_03 data
# Use Analysis block 1
###############################################################################################################

DATA=23
PANEL=3
data_dir="CK_2016-06-23_03"

file_panel="panel3.xlsx"
file_metadata="metadata_23_03.xlsx"

pca_score_cutoff=0.9
rand_seed_consensus=1234
nmetaclusts=20

prefix_data="23_"
prefix_panel="03_"
prefix_pca="pca1_"
prefix_clust="cl20_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts}


rand_seed_consensus=1234
nmetaclusts=6
prefix_clust="cl6_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts}


# --------------------------------------------------
# Analysis of CK_2016-06-23_03 merging2
# Use Analysis block 2
# --------------------------------------------------

file_merging="${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_merging2.xlsx"
prefix_merging="merging2_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging}

file_merging="${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_merging3.xlsx"
prefix_merging="merging3_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging}

# --------------------------------------------------
# Analysis of CK_2016-06-23_03all data using panel3_v2.xlsx
# Use Analysis block 1
# --------------------------------------------------

DATA=23
PANEL=3
data_dir="CK_2016-06-23_03all"

file_panel="panel3_v2.xlsx"
file_metadata="metadata_23_03all.xlsx"

pca_score_cutoff=1
nmetaclusts=20
rand_seed_consensus=1234

prefix_data="23_"
prefix_panel="03v2_"
prefix_pca="pca1_"
prefix_clust="cl20_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts}


rand_seed_consensus=1234
nmetaclusts=8
prefix_clust="cl8_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts}


# --------------------------------------------------
# Analysis of CK_2016-06-23_03all merging2
# Use Analysis block 2
# --------------------------------------------------

file_merging="${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_merging2.xlsx"
prefix_merging="merging2_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging}

# --------------------------------------------------
# Analysis of CK_2016-06-23_03all data using panel3_v3.xlsx
# Use Analysis block 1
# --------------------------------------------------

DATA=23
PANEL=3
data_dir="CK_2016-06-23_03all"

file_panel="panel3_v3.xlsx"
file_metadata="metadata_23_03all.xlsx"

pca_score_cutoff=1
nmetaclusts=20
rand_seed_consensus=1234

prefix_data="23_"
prefix_panel="03v3_"
prefix_pca="pca1_"
prefix_clust="cl20_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts}


rand_seed_consensus=1234
nmetaclusts=8
prefix_clust="cl8_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts}



###############################################################################################################
# Analysis of CK_2016-06-29_01 data
# Use Analysis block 1
###############################################################################################################

DATA=29
PANEL=1
data_dir="CK_2016-06-29_01"

file_panel="panel1.xlsx"
file_metadata="metadata_29_01.xlsx"

pca_score_cutoff=3.5
rand_seed_consensus=1234
nmetaclusts=20

prefix_data="29_"
prefix_panel="01_"
prefix_pca="pca1_"
prefix_clust="cl20_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts}

# --------------------------------------------------
# Analysis of CK_2016-06-29_01 merging3
# Use Analysis block 2
# --------------------------------------------------

file_merging="${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_merging2.xlsx"
prefix_merging="merging2_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging}


file_merging="${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_merging3.xlsx"
prefix_merging="merging3_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging}

# --------------------------------------------------
# CK_2016-06-29_01 - CD4 and CD8 cluster extracting from merging2
# Use Analysis block 3
# --------------------------------------------------

prefix_merging="merging2_"

extract_cluster=('CD4' 'CD8')
extract_dir=('CK_2016-06-29_01_CD4_merging2' 'CK_2016-06-29_01_CD8_merging2')

for i in 0 1
do
  ./Analysis_block_3_cluster_extracting.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_extracting ${cluster_extracting} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging} --extract_cluster ${extract_cluster[$i]} --extract_dir ${extract_dir[$i]}
done

# --------------------------------------------------
# Analysis of CK_2016-06-29_01_CD4_merging2 using panel1CD4.xlsx
# and CK_2016-06-29_01_CD8_merging2 using panel1CD8.xlsx
# Use Analysis block 1
# --------------------------------------------------

DATA=29
PANEL=1

file_metadata="metadata_29_01.xlsx"

rand_seed_consensus=1234
nmetaclusts=20

prefix_pca="pca1_"
prefix_clust="cl20_"

data_dir=('CK_2016-06-29_01_CD4_merging2' 'CK_2016-06-29_01_CD8_merging2')
prefix_data=('29CD4_' '29CD8_')
file_panel=('panel1CD4.xlsx' 'panel1CD8.xlsx')
prefix_panel=('01CD4_' '01CD8_')
pca_score_cutoff=(1.9 1.9)

for i in 0 1
do
  ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff[$i]} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts}
done


# --------------------------------------------------
# Analysis of CK_2016-06-29_01_CD4_merging2 using panel1CD4.xlsx of cluster_merging2
# and CK_2016-06-29_01_CD8_merging2 using panel1CD8.xlsx of cluster_merging2
# Use Analysis block 2
# --------------------------------------------------

file_merging=("${prefix_data[0]}${prefix_panel[0]}${prefix_pca}${prefix_clust}cluster_merging2.xlsx" "${prefix_data[1]}${prefix_panel[1]}${prefix_pca}${prefix_clust}cluster_merging2.xlsx")
prefix_merging=('merging2_' 'merging2_')

for i in 0 1
do
  ./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging[$i]} --file_merging ${file_merging[$i]}
done


###############################################################################################################
# Analysis of CK_2016-06-29_02 data
# Use Analysis block 1
###############################################################################################################

DATA=29
PANEL=2
data_dir="CK_2016-06-29_02"

file_panel="panel2.xlsx"
file_metadata="metadata_29_02.xlsx"

pca_score_cutoff=1
rand_seed_consensus=123
nmetaclusts=20

prefix_data="29_"
prefix_panel="02_"
prefix_pca="pca1_"
prefix_clust="cl20_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts}

# --------------------------------------------------
# Analysis of CK_2016-06-29_02 merging
# Use Analysis block 2
# --------------------------------------------------

file_merging="${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_merging.xlsx"
prefix_merging="merging_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging}

# --------------------------------------------------
# CK_2016-06-29_02 - CD4 and CD8 cluster extracting from merging
# Use Analysis block 3
# --------------------------------------------------

prefix_merging="merging_"

extract_cluster=('CD4' 'CD8')
extract_dir=('CK_2016-06-29_02_CD4_merging' 'CK_2016-06-29_02_CD8_merging')

for i in 0 1
do
  ./Analysis_block_3_cluster_extracting.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_extracting ${cluster_extracting} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging} --extract_cluster ${extract_cluster[$i]} --extract_dir ${extract_dir[$i]}
done

# --------------------------------------------------
# Analysis of CK_2016-06-29_02_CD4_merging using panel2CD4.xlsx
# and CK_2016-06-29_02_CD8_merging using panel2CD8.xlsx
# Use Analysis block 1
# --------------------------------------------------

DATA=29
PANEL=2

file_metadata="metadata_29_02.xlsx"

rand_seed_consensus=1234
nmetaclusts=20

prefix_pca="pca1_"
prefix_clust="cl20_"

data_dir=('CK_2016-06-29_02_CD4_merging' 'CK_2016-06-29_02_CD8_merging')
prefix_data=('29CD4_' '29CD8_')
file_panel=('panel2CD4.xlsx' 'panel2CD8.xlsx')
prefix_panel=('02CD4_' '02CD8_')
pca_score_cutoff=(1.2 1)

for i in 0 1
do
  ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff[$i]} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts}
done


# --------------------------------------------------
# Analysis of CK_2016-06-29_02_CD4_merging2 using panel2CD4.xlsx of cluster_merging2
# and CK_2016-06-29_02_CD8_merging2 using panel2CD8.xlsx of cluster_merging2
# Use Analysis block 2
# --------------------------------------------------

file_merging=("${prefix_data[0]}${prefix_panel[0]}${prefix_pca}${prefix_clust}cluster_merging2.xlsx" "${prefix_data[1]}${prefix_panel[1]}${prefix_pca}${prefix_clust}cluster_merging2.xlsx")
prefix_merging=('merging2_' 'merging2_')

for i in 0 1
do
  ./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging[$i]} --file_merging ${file_merging[$i]}
done

###############################################################################################################
# Analysis of CK_2016-06-29_03 data
# Use Analysis block 1
###############################################################################################################

DATA=29
PANEL=3
data_dir="CK_2016-06-29_03"

file_panel="panel3.xlsx"
file_metadata="metadata_29_03.xlsx"

pca_score_cutoff=0.9
rand_seed_consensus=1234
nmetaclusts=20

prefix_data="29_"
prefix_panel="03_"
prefix_pca="pca1_"
prefix_clust="cl20_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts}

# --------------------------------------------------
# Analysis of CK_2016-06-29_03 merging
# Use Analysis block 2
# --------------------------------------------------

file_merging="${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_merging4.xlsx"
prefix_merging="merging4_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging}

# --------------------------------------------------
# Analysis of CK_2016-06-29_03 data - using observables selected in CK_2016-06-23_03 with pca1 approach
# Use Analysis block 1
# --------------------------------------------------

DATA=29
PANEL=3
data_dir="CK_2016-06-29_03"

file_panel="panel3.xlsx"
file_metadata="metadata_29_03.xlsx"

pca_score_cutoff=0.9 # cutoff does not have to be specified bcs I use the clustering_observables.xls
rand_seed_consensus=1234
nmetaclusts=20

prefix_data="29_"
prefix_panel="03_"
prefix_pca="pca1v23_"
prefix_clust="cl20_"

# Use observables selected in 23_03_pca1

cp $RWD_MAIN/CK_2016-06-23_03/030_heatmaps/23_03_pca1_clustering_observables.xls $RWD_MAIN/CK_2016-06-29_03/030_heatmaps/29_03_pca1v23_clustering_observables.xls

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin}

# --------------------------------------------------
# Analysis of CK_2016-06-29_03 data - using a lower PCA score cutoff
# Use Analysis block 1
# --------------------------------------------------

DATA=29
PANEL=3
data_dir="CK_2016-06-29_03"

file_panel="panel3.xlsx"
file_metadata="metadata_29_03.xlsx"

pca_score_cutoff=0.3 # lower cutoff than for pca1
rand_seed_consensus=1234
nmetaclusts=20

prefix_data="29_"
prefix_panel="03_"
prefix_pca="pca1b_"
prefix_clust="cl20_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization false --pcascores false --select_observables ${select_observables} --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin}

































#
