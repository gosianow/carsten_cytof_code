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


prefix_clust=('cl5_' 'cl8_')

file_merging=("${prefix_data[0]}${prefix_panel[0]}${prefix_pca}${prefix_clust[0]}cluster_merging5.xlsx" "${prefix_data[1]}${prefix_panel[1]}${prefix_pca}${prefix_clust[1]}cluster_merging4.xlsx")
prefix_merging=('merging5_' 'merging4_')

# for i in 1
# do
#   ./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust[$i]} --prefix_merging ${prefix_merging[$i]} --file_merging ${file_merging[$i]}
# done

prefix_clust=('cl5_' 'cl8_')

file_merging=("${prefix_data[0]}${prefix_panel[0]}${prefix_pca}${prefix_clust[0]}cluster_merging5.xlsx" "${prefix_data[1]}${prefix_panel[1]}${prefix_pca}${prefix_clust[1]}cluster_merging5.xlsx")
prefix_merging=('merging5_' 'merging5_')

# for i in  1 # run only for CD8
# do
#   ./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust[$i]} --prefix_merging ${prefix_merging[$i]} --file_merging ${file_merging[$i]}
# done


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


prefix_clust=('cl5_' 'cl8_')

file_merging=("${prefix_data[0]}${prefix_panel[0]}${prefix_pca}${prefix_clust[0]}cluster_merging5.xlsx" "${prefix_data[1]}${prefix_panel[1]}${prefix_pca}${prefix_clust[1]}cluster_merging4.xlsx")
prefix_merging=('merging5_' 'merging4_')

for i in 1
do
  ./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust[$i]} --prefix_merging ${prefix_merging[$i]} --file_merging ${file_merging[$i]}
done















#
