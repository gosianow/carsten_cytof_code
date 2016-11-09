#!/bin/bash

###############################################################################################################
## Define paths to software and reference files

RWD_MAIN=/Users/gosia/Dropbox/UZH/carsten_cytof
RCODE=/Users/gosia/Dropbox/UZH/carsten_cytof_code
METADATA=$RWD_MAIN/CK_metadata
PANELS=$RWD_MAIN/CK_panels


## Define which analysis to re-run

# data_normalization=true
# pcascores=true
# select_observables=false
# flowsom=false
# flowsom_validation=false
# heatmaps=false
# runtsne=false
# plottsne=false
# plottsne_expr=false
# frequencies=false
# expression=false
# cluster_merging=false
# cluster_extracting=false
# fcs_saving=false
# cytokines_bimatrix=false
# pd1_bimatrix=false
# cytokines_expression=false
# pd1_expression=false


data_normalization=false
pcascores=false
select_observables=true
flowsom=true
flowsom_validation=true
heatmaps=true
runtsne=true
plottsne=true
plottsne_expr=true
frequencies=true
expression=true
cluster_merging=true
cluster_extracting=true
fcs_saving=false
cytokines_bimatrix=false
pd1_bimatrix=false
cytokines_expression=false
pd1_expression=false




## global parameters
tsne_pmin=1500


###############################################################################################################
# Analysis of CK_2016-06-29_03all2 data
# Use Analysis block 1
###############################################################################################################

# --------------------------------------------------
# Analysis of CK_2016-06-29_03all2_myeloid_merging3 data using panel3.xlsx
# Use Analysis block 1
# --------------------------------------------------

DATA=29
PANEL=3
data_dir="CK_2016-06-29_03all2_myeloid_merging3"

file_panel="panel3.xlsx"
file_metadata="metadata_29_03all2.xlsx"

pca_score_cutoff=1.6
nmetaclusts=20
rand_seed_consensus=1234

prefix_data="29mye_"
prefix_panel="03_"
prefix_pca="pca1_"
prefix_clust="cl20_"

rand_seed_consensus=1234

for i in 5
do
  nmetaclusts=$i
  prefix_clust="cl${i}_"

  ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --flowsom_validation false --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --plottsne_expr false --frequencies ${frequencies} --expression false --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin}
done

















#
