#!/bin/bash

###############################################################################################################
## Define paths to software and reference files

# RWD_MAIN=/Users/gosia/Dropbox/UZH/carsten_cytof
# RCODE=/Users/gosia/Dropbox/UZH/carsten_cytof_code
RWD_MAIN=/home/Shared/data/cytof/carsten_cytof
RCODE=/home/gosia/R/carsten_cytof_code

METADATA=$RWD_MAIN/FACSvalidation2_metadata
PANELS=$RWD_MAIN/FACSvalidation2_panels


## Define which analysis to re-run
data_normalization=false
pcascores=false
select_observables=true
flowsom=true
flowsom_validation=false
heatmaps=true
runtsne=true
plottsne=true
plottsne_expr=false
frequencies=false
expression=false
cluster_merging=false
cluster_extracting=false

## global parameters
tsne_pmin=3000

###############################################################################################################
# Analysis of FACSvalidation2 data
# Use Analysis block 1
###############################################################################################################


# ------------------------
# Use all markers - pca_score_cutoff = 0
# ------------------------


DATA=23
PANEL=1
data_dir="FACSvalidation2"

file_panel="panel_facs_validation2.xlsx"
file_metadata="metadata_facs_validation2.xlsx"

pca_score_cutoff=0
rand_seed_consensus=1234
nmetaclusts=20

prefix_data="fv_"
prefix_panel="fv_"
prefix_pca="pca0_"
prefix_clust="cl20_"


./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --flowsom_validation ${flowsom_validation} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --plottsne_expr ${plottsne_expr} --frequencies ${frequencies} --expression false --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin}


for i in 4 6
do
  nmetaclusts=$i
  prefix_clust="cl${i}_"

  ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --flowsom_validation false --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --plottsne_expr false --frequencies ${frequencies} --expression false --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin}

done

# ------------------------
# pca_score_cutoff = 6
# ------------------------

DATA=23
PANEL=1
data_dir="FACSvalidation2"

file_panel="panel_facs_validation2.xlsx"
file_metadata="metadata_facs_validation2.xlsx"

pca_score_cutoff=6
rand_seed_consensus=1234
nmetaclusts=20

prefix_data="fv_"
prefix_panel="fv_"
prefix_pca="pca1_"
prefix_clust="cl20_"


./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --flowsom_validation ${flowsom_validation} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --plottsne_expr ${plottsne_expr} --frequencies ${frequencies} --expression false --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin}








#
