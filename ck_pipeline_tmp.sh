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
plottsne=true
frequencies=false
expression=false
cluster_merging=false
cluster_extracting=false
cytokines=true
fcs_saving=false
pd1=false


# --------------------------------------------------
# Analysis of CK_2016-06-23_02_CD4_merging2 using panel2CD4.xlsx
# and CK_2016-06-23_02_CD8_merging2 using panel2CD8.xlsx
# Use Analysis block 1
# --------------------------------------------------

DATA=23
PANEL=2
file_metadata="metadata_23_02.xlsx"

rand_seed_consensus=1234
nmetaclusts=20

prefix_pca="pca1_"
prefix_clust="cl20_"

data_dir=('CK_2016-06-23_02_CD4_merging2' 'CK_2016-06-23_02_CD8_merging2')
prefix_data=('23CD4_' '23CD8_')
file_panel=('panel2CD4.xlsx' 'panel2CD8.xlsx')
prefix_panel=('02CD4_' '02CD8_')
pca_score_cutoff=(0.96 1)

# --------------------------------------------------
# Analysis of cytokines
# in CK_2016-06-23_02_CD4_merging2 and CK_2016-06-23_02_CD8_merging2
# Use Analysis block 4
# --------------------------------------------------

prefix_merging=('merging_' 'merging_')
clsubset=("c('CM','EM')" "c('CM','EM','TE')")
prefix_clsubset=('Tmem_' 'Tmem_')
file_cytokines_cutoffs=('panel2CD4_cytokines_CM.xlsx' 'panel2CD8_cytokines_CM.xlsx')
prefix_cytokines_cutoffs=('cytCM_' 'cytCM_')
nmetaclusts=(40 20)

file_merging_cyt=("${prefix_data[0]}${prefix_panel[0]}${prefix_pca}${prefix_merging[0]}${prefix_clsubset[0]}${prefix_cytokines_cutoffs[0]}raw2_cl${nmetaclusts[0]}_cytokine_merging.xlsx" "${prefix_data[1]}${prefix_panel[1]}${prefix_pca}${prefix_merging[1]}${prefix_clsubset[1]}${prefix_cytokines_cutoffs[1]}raw2_cl${nmetaclusts[1]}_cytokine_merging.xlsx")
prefix_merging_cyt=("cytmerging_" "cytmerging_")


for i in 0 1
do
  ./Analysis_block_4_cytokines.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cytokines ${cytokines} --PANELS ${PANELS} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging[$i]} --prefix_clsubset ${prefix_clsubset[$i]} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs[$i]} --file_cytokines_cutoffs ${file_cytokines_cutoffs[$i]} --clsubset ${clsubset[$i]} --nmetaclusts ${nmetaclusts[$i]} --file_merging_cyt ${file_merging_cyt[$i]} --prefix_merging_cyt ${prefix_merging_cyt[$i]}
done









#
