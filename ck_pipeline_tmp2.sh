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
flowsom_validation=false
heatmaps=false
runtsne=false
plottsne=false
plottsne_expr=false
frequencies=false
expression=false
cluster_merging=false
cluster_extracting=false
fcs_saving=false
cytokines_bimatrix_main=true
cytokines_bimatrix_cluster_merging=true
pd1_bimatrix=false
cytokines_expression=false
pd1_expression=false

## global parameters
tsne_pmin=1500

# -----------------------------
### Tmem cluster
# -----------------------------

DATA=29
PANEL=2
file_metadata="metadata_29_02.xlsx"

prefix_pca="pca1_"
prefix_clust="cl20_"

data_dir=('CK_2016-06-29_02_CD4_merging' 'CK_2016-06-29_02_CD8_merging')
prefix_data=('29CD4_' '29CD8_')
prefix_panel=('02CD4_' '02CD8_')

prefix_merging=('merging3_' 'merging3_') # name of merging from which the Tmem clusters are extracted
clsubset=("c('CM','EM')" "c('CM','EM','TE')")
prefix_clsubset=('Tmem_' 'Tmem_')
file_cytokines_cutoffs=('panel2CD4_29_cytokines_CM.xlsx' 'panel2CD8_29_cytokines_CM.xlsx')
prefix_cytokines_cutoffs=('cytCM_' 'cytCM_')


nmetaclusts=(40 20)

## Merging is done only when these files exist
file_merging_cyt=("${prefix_data[0]}${prefix_panel[0]}${prefix_pca}${prefix_merging[0]}${prefix_clsubset[0]}${prefix_cytokines_cutoffs[0]}raw2_cl${nmetaclusts[0]}_cytokine_merging2.xlsx" "${prefix_data[1]}${prefix_panel[1]}${prefix_pca}${prefix_merging[1]}${prefix_clsubset[1]}${prefix_cytokines_cutoffs[1]}raw2_cl${nmetaclusts[1]}_cytokine_merging2.xlsx")
prefix_merging_cyt=("cytmerging2_" "cytmerging2_")

for i in 0 1
do
  ./Analysis_block_4_cytokines_bimatrix_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cytokines_bimatrix_cluster_merging ${cytokines_bimatrix_cluster_merging} --PANELS ${PANELS} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging[$i]} --prefix_clsubset ${prefix_clsubset[$i]} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs[$i]} --nmetaclusts ${nmetaclusts[$i]} --file_merging_cyt ${file_merging_cyt[$i]} --prefix_merging_cyt ${prefix_merging_cyt[$i]}
done



















#
