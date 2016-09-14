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


## global parameters
tsne_pmin=1500


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


# --------------------------------------------------
# Analysis of CK_2016-06-23_03 merging2
# Use Analysis block 2
# --------------------------------------------------

prefix_clust="cl20_"

file_merging="${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_merging4.xlsx"
prefix_merging="merging4_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging}





























#
