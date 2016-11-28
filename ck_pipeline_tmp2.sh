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
cytokines_bimatrix=false
cytokines_bimatrix_main=false
cytokines_bimatrix_cluster_merging=false
pd1_bimatrix=false
pd1_bimatrix_main=true
pd1_bimatrix_cluster_merging=false
cytokines_expression=false
pd1_expression=false

## global parameters
tsne_pmin=1500

# --------------------------------------------------
# Analysis of PD1+ and PD1- cells based on bimatrix
# in CK_2016-06-23_02_CD4_merging2 and CK_2016-06-23_02_CD8_merging2
# Use Analysis block 5
# --------------------------------------------------

### Analysis of cytokines for PD1+

DATA=23
PANEL=2
file_metadata="metadata_23_02.xlsx"

prefix_pca="pca1_"
prefix_clust="cl20_"

data_dir=('CK_2016-06-23_02_CD4_merging2' 'CK_2016-06-23_02_CD8_merging2')
prefix_data=('23CD4_' '23CD8_')
file_panel=('panel2CD4.xlsx' 'panel2CD8.xlsx')
prefix_panel=('02CD4_' '02CD8_')

prefix_merging=('merging2_' 'merging2_')
clsubset=("c('CM','EM','TE','TM')" "c('CM','EM','TE','TM')")
prefix_clsubset=('Tmem_' 'Tmem_')
file_cytokines_cutoffs=('panel2CD4_23_cytokines_CM.xlsx' 'panel2CD8_23_cytokines_CM_new.xlsx')
prefix_cytokines_cutoffs=('cytCM_raw2_pd1_' 'cytCMnew_raw2_pd1_')

outdir="070_pd1_bimatrix"

som_dim=(5 5)
nmetaclusts=(25 25)

for i in 0 1
do
  ./Analysis_block_5_pd1_bimatrix_0_bimatrix.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --pd1_bimatrix ${pd1_bimatrix} --PANELS ${PANELS} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging[$i]} --prefix_clsubset ${prefix_clsubset[$i]} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs[$i]} --file_cytokines_cutoffs ${file_cytokines_cutoffs[$i]} --clsubset ${clsubset[$i]} --outdir ${outdir}

  ./Analysis_block_6_cytokines_bimatrix_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cytokines_bimatrix_main ${pd1_bimatrix_main} --PANELS ${PANELS} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging[$i]} --prefix_clsubset ${prefix_clsubset[$i]} --prefix_cytokines_cutoffs "${prefix_cytokines_cutoffs[$i]}positive_"  --som_dim ${som_dim[$i]} --nmetaclusts ${nmetaclusts[$i]} --outdir ${outdir}

done




#
