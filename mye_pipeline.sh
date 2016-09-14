#!/bin/bash

###############################################################################################################
## Define paths to software and reference files

RWD_MAIN=/Users/gosia/Dropbox/UZH/carsten_cytof
RCODE=/Users/gosia/Dropbox/UZH/carsten_cytof_code
METADATA=$RWD_MAIN/MyeEUNITER_metadata
PANELS=$RWD_MAIN/MyeEUNITER_panels


## Define which analysis to re-run
data_normalization=false
pcascores=false
select_observables=false
flowsom=false
heatmaps=false
runtsne=false
plottsne=false
frequencies=false # there is no frequency analysis here!
expression=false # there is no expression analysis here!
cluster_merging=false
cluster_extracting=false

## global parameters
tsne_pmin=5000 # In the CK analysis, I use 1500 per sample.

###############################################################################################################
# Analysis of MyeEUNITER data
# Use Analysis block 1
###############################################################################################################

data_dir="MyeEUNITER"

file_panel="panel_MyeEUNITER.xlsx"
file_metadata="metadata_MyeEUNITER.xlsx"

pca_score_cutoff=0 # We keep all the markers!
rand_seed_consensus=1234
nmetaclusts=20

prefix_data="mye_"
prefix_panel="mye_"
prefix_pca="pca1_"
prefix_clust="cl20_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin}


# --------------------------------------------------
# Analysis of MyeEUNITER cluster_merging
# Use Analysis block 2
# --------------------------------------------------

file_merging="${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_merging2.xlsx"
prefix_merging="merging2_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging}


file_merging="${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_merging3.xlsx"
prefix_merging="merging3_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging}





# --------------------------------------------------
# MyeEUNITER - CD66b+ cluster extracting from merging
# Use Analysis block 3
# --------------------------------------------------

prefix_merging="merging2_"

extract_cluster="CD66b+"
extract_dir="MyeEUNITER_CD66b_merging2"


./Analysis_block_3_cluster_extracting.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_extracting ${cluster_extracting} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging} --extract_cluster ${extract_cluster} --extract_dir ${extract_dir}


# --------------------------------------------------
# Analysis of MyeEUNITER_CD66b_merging2
# Use Analysis block 1
# --------------------------------------------------


data_dir="MyeEUNITER_CD66b_merging2"

file_panel="panel_MyeEUNITER.xlsx"
file_metadata="metadata_MyeEUNITER.xlsx"

pca_score_cutoff=0 # We keep all the markers!
rand_seed_consensus=1234
nmetaclusts=20

prefix_data="myeCD66b_"
prefix_panel="mye_"
prefix_pca="pca1_"
prefix_clust="cl20_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin}


















#
