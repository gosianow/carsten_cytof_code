#!/bin/bash
## Define paths to software and reference files

RWD_MAIN=/Users/gosia/Dropbox/UZH/carsten_cytof
RCODE=$RWD_MAIN

## Define which analysis to re-run
pcascores=false
select_observables=false
flowsom=false
cluster_merging=false
heatmaps=false
runtsne=false
plottsne=false
frequencies=false
cluster_merging=false
cluster_extracting=false
cytokines=true
fcs_saving=false


##########################################
# Analysis of CK_2016-06-23_02_CD4_merging2 and CK_2016-06-23_02_CD8_merging2
# for cluster_merging_CD4.xlsx and cluster_merging_CD8.xlsx
##########################################

data='CK_2016-06-23_02'
merging='merging2'
extr_name=('CD4' 'CD8')
extr_dir=('CD4' 'CD8')

clusters2analyse="c('CM','EM','TE')"
cluster_name="Tmem"
nmetaclusts=40

if ${cytokines}; then
for indx in 0
do

  RWD=$RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"
  
  pca_prefix="pnl${extr_dir[$indx]}_pca1_" # 'pnlCD4_pca1_'
  merging_prefix="merging_${extr_dir[$indx]}" # 'merging_CD4'
  
  ### Analysis of positive-negative (cytokine) markers
  
  ## based on cytokines_CM_RAW.xlsx
 R CMD BATCH --no-save --no-restore "--args rwd='$RWD' cytokines_prefix='${pca_prefix}${merging_prefix}_cyt${cluster_name}_' path_panel='panel_${extr_dir[$indx]}.xlsx' path_cytokines_cutoffs='panel_${extr_dir[$indx]}_cytokines_CM_RAW.xlsx' path_clustering='${pca_prefix}${merging_prefix}_clustering.xls' path_clustering_labels='${pca_prefix}${merging_prefix}_clustering_labels.xls' clusters2analyse=${clusters2analyse} cutoff_colname='positive_cutoff_raw' data2analyse='raw' cytokines_suffix='_${nmetaclusts}cl_raw' nmetaclusts=${nmetaclusts}" $RCODE/06_cytokines.R $ROUT/06_cytokines.Rout
 tail $ROUT/06_cytokines.Rout
  

done
fi




if ${cytokines}; then
for indx in 1
do

  RWD=$RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"
  
  pca_prefix="pnl${extr_dir[$indx]}_pca1_" # 'pnlCD4_pca1_'
  merging_prefix="merging_${extr_dir[$indx]}_2" # 'merging_CD4_2'
  
  ### Analysis of positive-negative (cytokine) markers
  
  ## based on cytokines_CM_RAW.xlsx
 R CMD BATCH --no-save --no-restore "--args rwd='$RWD' cytokines_prefix='${pca_prefix}${merging_prefix}_cyt${cluster_name}_' path_panel='panel_${extr_dir[$indx]}.xlsx' path_cytokines_cutoffs='panel_${extr_dir[$indx]}_cytokines_CM_RAW.xlsx' path_clustering='${pca_prefix}${merging_prefix}_clustering.xls' path_clustering_labels='${pca_prefix}${merging_prefix}_clustering_labels.xls' clusters2analyse=${clusters2analyse} cutoff_colname='positive_cutoff_raw' data2analyse='raw' cytokines_suffix='_${nmetaclusts}cl_raw' nmetaclusts=${nmetaclusts}" $RCODE/06_cytokines.R $ROUT/06_cytokines.Rout
 tail $ROUT/06_cytokines.Rout
  

done
fi




