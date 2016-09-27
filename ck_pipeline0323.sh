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
plottsne_expr=false
frequencies=false
expression=false
cluster_merging=false
cluster_extracting=false


## global parameters
tsne_pmin=3000

###############################################################################################################
# Analysis of CK_2016-03-23_01 data
# Use Analysis block 1
###############################################################################################################

DATA=0323
PANEL=1
data_dir="CK_2016-03-23_01"

file_panel="panel_0323_1_Tcells.xlsx"
file_metadata="metadata_0323_01.xlsx"

pca_score_cutoff=2.5
rand_seed_consensus=1234
nmetaclusts=20

prefix_data="0323_"
prefix_panel="01_"
prefix_pca="pca1_"
prefix_clust="cl20_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --plottsne_expr ${plottsne_expr} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin}

rand_seed_consensus=1234
nmetaclusts=6
prefix_clust="cl6_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --plottsne_expr false --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin}


rand_seed_consensus=1234
nmetaclusts=7
prefix_clust="cl7_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --plottsne_expr false --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin}

rand_seed_consensus=1234
nmetaclusts=8
prefix_clust="cl8_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --plottsne_expr false --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin}



# --------------------------------------------------
# Analysis of CK_2016-03-23_01 cluster_merging
# Use Analysis block 2
# --------------------------------------------------

prefix_clust="cl20_"

file_merging="${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_mergingNEW.xlsx"
prefix_merging="mergingNEW_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging}


# --------------------------------------------------
# CK_2016-03-23_01 - CD4 and CD8 cluster extracting from mergingNEW2
# Use Analysis block 3
# --------------------------------------------------

prefix_merging="mergingNEW2_"

extract_cluster=('CD4' 'CD8')
extract_dir=('CK_2016-03-23_01_CD4_mergingNEW2' 'CK_2016-03-23_01_CD8_mergingNEW2')

for i in 0 1
do
  ./Analysis_block_3_cluster_extracting.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_extracting ${cluster_extracting} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging} --extract_cluster ${extract_cluster[$i]} --extract_dir ${extract_dir[$i]}
done


# --------------------------------------------------
# Analysis of CK_2016-03-23_01_CD4_mergingNEW2 using panel1CD4.xlsx
# and CK_2016-03-23_01_CD8_mergingNEW2 using panel1CD8.xlsx
# Use Analysis block 1
# --------------------------------------------------

DATA=0323
PANEL=1

file_metadata="metadata_0323_01.xlsx"

rand_seed_consensus=1234
nmetaclusts=20

prefix_pca="pca1_"
prefix_clust="cl20_"

data_dir=('CK_2016-03-23_01_CD4_mergingNEW2' 'CK_2016-03-23_01_CD8_mergingNEW2')
prefix_data=('23CD4_' '23CD8_')
file_panel=('panel1CD4.xlsx' 'panel1CD8.xlsx')
prefix_panel=('01CD4_' '01CD8_')
pca_score_cutoff=(2 2)

for i in 0 1
do
  ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --plottsne_expr ${plottsne_expr} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff[$i]} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin}
done


rand_seed_consensus=12344
nmetaclusts=(6 7)
prefix_clust=("cl6_" "cl7_")

for i in 0 1
do
  ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --plottsne_expr false --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust[$i]} --pca_score_cutoff ${pca_score_cutoff[$i]} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts[$i]}
done


# --------------------------------------------------
# Analysis of CK_2016-03-23_01_CD4_mergingNEW2 using panel1CD4.xlsx for cluster_merging
# and CK_2016-03-23_01_CD8_mergingNEW2 using panel1CD8.xlsx for cluster_merging
# Use Analysis block 2
# --------------------------------------------------

prefix_clust="cl20_"

file_merging=("${prefix_data[0]}${prefix_panel[0]}${prefix_pca}${prefix_clust}cluster_merging.xlsx" "${prefix_data[1]}${prefix_panel[1]}${prefix_pca}${prefix_clust}cluster_merging.xlsx")
prefix_merging=('merging_' 'merging_')

for i in 0 1
do
  ./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging[$i]} --file_merging ${file_merging[$i]}
done



###############################################################################################################
# Analysis of CK_2016-03-23_03 data
# Use Analysis block 1
###############################################################################################################

DATA=0323
PANEL=3
data_dir="CK_2016-03-23_03"

file_panel="panel_0323_3_myeloid.xlsx"
file_metadata="metadata_0323_03.xlsx"

pca_score_cutoff=0.97
rand_seed_consensus=1234
nmetaclusts=20

prefix_data="0323_"
prefix_panel="03_"
prefix_pca="pca1_"
prefix_clust="cl20_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --plottsne_expr ${plottsne_expr} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin}


rand_seed_consensus=1234
nmetaclusts=5
prefix_clust="cl5_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --plottsne_expr false --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin}


rand_seed_consensus=1234
nmetaclusts=6
prefix_clust="cl6_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --plottsne_expr false --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin}


rand_seed_consensus=1234
nmetaclusts=9
prefix_clust="cl9_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --plottsne_expr false --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin}



# --------------------------------------------------
# Analysis of CK_2016-03-23_03 merging2
# Use Analysis block 2
# --------------------------------------------------

prefix_clust="cl20_"

file_merging="${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_merging2.xlsx"
prefix_merging="merging2_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging}










#
