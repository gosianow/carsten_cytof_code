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
cluster_merging=true
cluster_extracting=false
fcs_saving=false
cytokines_bimatrix=true
pd1_bimatrix=false
cytokines_expression=false
pd1_expression=false

## global parameters
tsne_pmin=1500

# --------------------------------------------------
# Analysis of cytokines based on bimatrix
# in CK_2016-06-23_02_CD4_merging2 and CK_2016-06-23_02_CD8_merging2
# Use Analysis block 4
# --------------------------------------------------

# --------------------------------------------------
# Analysis of CK_2016-03-23_01_CD4_merging using panel1CD4.xlsx
# and CK_2016-03-23_01_CD8_merging using panel1CD8.xlsx
# Use Analysis block 1
# --------------------------------------------------

DATA=0323
PANEL=1

file_metadata="metadata_0323_01.xlsx"

rand_seed_consensus=1234
nmetaclusts=20

prefix_pca="pca1_"
prefix_clust="cl20_"

data_dir=('CK_2016-03-23_01_CD4_merging' 'CK_2016-03-23_01_CD8_merging')
prefix_data=('0323CD4_' '0323CD8_')
file_panel=('panel1CD4.xlsx' 'panel1CD8.xlsx')
prefix_panel=('01CD4_' '01CD8_')
pca_score_cutoff=(1.1 1.28)


nmetaclusts=(8 5)
prefix_clust=("cl8_" "cl5_")

for i in 0 # run only for CD4
do
  ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --flowsom_validation false --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --plottsne_expr false --frequencies false --expression false --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust[$i]} --pca_score_cutoff ${pca_score_cutoff[$i]} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts[$i]}
done




















#
