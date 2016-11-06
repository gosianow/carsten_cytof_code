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
flowsom=true
flowsom_validation=false
heatmaps=true
runtsne=false
plottsne=true
plottsne_expr=false
frequencies=true
expression=true
cluster_merging=false
cluster_extracting=false
fcs_saving=false
cytokines_bimatrix_main=false
pd1_bimatrix=false
cytokines_expression=false
pd1_expression=false

## global parameters
tsne_pmin=1500


# ------------------------
# Use all markers - pca_score_cutoff = 0
# ------------------------

DATA=23
PANEL=1

file_metadata="metadata_23_01.xlsx"

rand_seed_consensus=1234
nmetaclusts=20

prefix_pca="pca0_"
prefix_clust="cl20_"

data_dir=('CK_2016-06-23_01_CD4_mergingNEW2' 'CK_2016-06-23_01_CD8_mergingNEW2')
prefix_data=('23CD4_' '23CD8_')
file_panel=('panel1CD4.xlsx' 'panel1CD8.xlsx')
prefix_panel=('01CD4_' '01CD8_')
pca_score_cutoff=(0 0)

# for i in 0 1
# do
#   ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --flowsom_validation ${flowsom_validation} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --plottsne_expr ${plottsne_expr} --frequencies ${frequencies} --expression false --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff[$i]} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin}
# done


# rand_seed_consensus=1234
# nmetaclusts=(5 5)
# prefix_clust=("cl5_" "cl5_")
#
# for i in 0 1
# do
#   ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --flowsom_validation false --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --plottsne_expr false --frequencies ${frequencies} --expression false --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust[$i]} --pca_score_cutoff ${pca_score_cutoff[$i]} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts[$i]}
# done





# ------------------------
# Use all markers - pca_score_cutoff = 0
# ------------------------


DATA=29
PANEL=1

file_metadata="metadata_29_01.xlsx"

rand_seed_consensus=1234
nmetaclusts=20

prefix_pca="pca0_"
prefix_clust="cl20_"

data_dir=('CK_2016-06-29_01_CD4_merging2' 'CK_2016-06-29_01_CD8_merging2')
prefix_data=('29CD4_' '29CD8_')
file_panel=('panel1CD4.xlsx' 'panel1CD8.xlsx')
prefix_panel=('01CD4_' '01CD8_')
pca_score_cutoff=(0 0)

# for i in 0 1
# do
#   ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --flowsom_validation ${flowsom_validation} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --plottsne_expr ${plottsne_expr} --frequencies ${frequencies} --expression false --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff[$i]} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin}
# done


rand_seed_consensus=1234
nmetaclusts=(5 5)
prefix_clust=("cl5_" "cl5_")

# for i in 0 1
# do
#   ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --flowsom_validation false --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --plottsne_expr false --frequencies ${frequencies} --expression false --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust[$i]} --pca_score_cutoff ${pca_score_cutoff[$i]} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts[$i]}
# done


rand_seed_consensus=1234
nmetaclusts=(6 5)
prefix_clust=("cl6_" "cl5_")

for i in 0 # run only for CD4
do
  ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --flowsom_validation false --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --plottsne_expr false --frequencies ${frequencies} --expression false --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust[$i]} --pca_score_cutoff ${pca_score_cutoff[$i]} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts[$i]}
done




#
