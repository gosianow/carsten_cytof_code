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
runtsne=true
plottsne=true
plottsne_expr=false
frequencies=true
expression=false
cluster_merging=false
cluster_extracting=false
fcs_saving=false
cytokines_bimatrix=false
pd1_bimatrix=false
cytokines_expression=false
pd1_expression=false

## global parameters
tsne_pmin=1500

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

# ./Analysis_block_1_main_norm_data.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --flowsom_validation ${flowsom_validation} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --plottsne_expr ${plottsne_expr} --frequencies ${frequencies} --expression false --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin}


rand_seed_consensus=1234

for i in 4 5
do

  echo ">>> reclustering into ${i} groups"
  nmetaclusts=$i
  prefix_clust="cl${i}_"

  ./Analysis_block_1_main_norm_data.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --flowsom_validation false --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --plottsne_expr false --frequencies ${frequencies} --expression false --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin}

done







#
