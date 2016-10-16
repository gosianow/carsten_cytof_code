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
heatmaps=true
runtsne=false
plottsne=true
plottsne_expr=false
frequencies=true
expression=true
cluster_merging=true
cluster_extracting=false
fcs_saving=false
cytokines_bimatrix=false
pd1_bimatrix=false
cytokines_expression=true
pd1_expression=false

## global parameters
tsne_pmin=1500

# --------------------------------------------------
# Analysis of cytokines based on expression
# in CK_2016-06-29_02_CD4_merging and CK_2016-06-23_02_CD8_merging
# Use Analysis block 4
# --------------------------------------------------

DATA=29
PANEL=2
file_metadata="metadata_29_02.xlsx"

prefix_pca="pca1_"

data_dir=('CK_2016-06-29_02_CD4_merging' 'CK_2016-06-29_02_CD8_merging')
prefix_data=('29CD4_' '29CD8_')
prefix_panel=('02CD4_' '02CD8_')

prefix_merging=('merging3_' 'merging3_') # name of merging from which the Tmem clusters are extracted
clsubset=("c('CM','EM')" "c('CM','EM','TE')")
prefix_clsubset=('Tmem_' 'Tmem_')
file_cytokines_cutoffs=('panel2CD4_29_cytokines_CM.xlsx' 'panel2CD8_29_cytokines_CM.xlsx')
prefix_cytokines_cutoffs=('cytCM_' 'cytCM_')


for i in 0 1
do
  ./Analysis_block_4_cytokines_expression.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cytokines ${cytokines_expression} --PANELS ${PANELS} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging[$i]} --prefix_clsubset ${prefix_clsubset[$i]} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs[$i]} --file_cytokines_cutoffs ${file_cytokines_cutoffs[$i]} --clsubset ${clsubset[$i]}
done





#
