#!/bin/bash

###############################################################################################################
## Define paths to software and reference files

RWD_MAIN=/Users/gosia/Dropbox/UZH/carsten_cytof
RCODE=/Users/gosia/Dropbox/UZH/carsten_cytof_code
METADATA=$RWD_MAIN/MyeEUNITER_metadata
PANELS=$RWD_MAIN/MyeEUNITER_panels


## Define which analysis to re-run
data_normalization=true
pcascores=true
select_observables=true
flowsom=true
heatmaps=true
runtsne=true
plottsne=true
frequencies=false # there is no frequency analysis here!
expression=false # there is no expression analysis here!
cluster_merging=true
cluster_extracting=true

## global parameters
tsne_pmin=5000 # In the CK analysis, I use 1500 per sample.

###############################################################################################################
# Analysis of MyeEUNITER data
# Use Analysis block 1
###############################################################################################################

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

# --------------------------------------------------
# Analysis of MyeEUNITER_CD66b_merging2 cluster_merging
# Use Analysis block 2
# --------------------------------------------------

file_merging="${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_merging.xlsx"
prefix_merging="merging_"

# ./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging}


# --------------------------------------------------
# MyeEUNITER_CD66b_merging2 - neutrophils cluster extracting from merging
# Use Analysis block 3
# --------------------------------------------------

prefix_merging="merging_"

extract_cluster="neutrophils"
extract_dir="MyeEUNITER_neutrophils_merging"


# ./Analysis_block_3_cluster_extracting.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_extracting ${cluster_extracting} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging} --extract_cluster ${extract_cluster} --extract_dir ${extract_dir}



# --------------------------------------------------
# Analysis of MyeEUNITER_CD66b_merging2
# Use Analysis block 1
# --------------------------------------------------


data_dir="MyeEUNITER_neutrophils_merging"

file_panel="panel_MyeEUNITER.xlsx"
file_metadata="metadata_MyeEUNITER.xlsx"

pca_score_cutoff=0 # We keep all the markers!
rand_seed_consensus=1234
nmetaclusts=20

prefix_data="myeNEUTROP_"
prefix_panel="mye_"
prefix_pca="pca1_"
prefix_clust="cl20_"

# ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin}

nmetaclusts=10
prefix_clust="cl10_"

# ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin}


nmetaclusts=5
prefix_clust="cl5_"

# ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin}




# --------------------------------------------------
# Analysis of MyeEUNITER_neutrophils_merging cluster_merging
# Use Analysis block 2
# --------------------------------------------------

prefix_clust="cl10_"
file_merging="${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_merging.xlsx"
prefix_merging="merging_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging}











#
