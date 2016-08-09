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
# Analysis of CK_2016-06-23_02_CD4_merging2 using panel_CD4.xlsx
# and CK_2016-06-23_02_CD8_merging2 using panel_CD8.xlsx
##########################################

data='CK_2016-06-23_02'
merging='merging2'
extr_name=('CD4' 'CD8')
extr_dir=('CD4' 'CD8')


if ${cytokines}; then
for indx in 0 1
do

  RWD=$RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"
  
  pca_prefix="pnl${extr_dir[$indx]}_pca1_"
  
  ### Analysis of positive-negative (cytokine) markers
  
  ## based on cytokines.xlsx
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' cytokines_prefix='${pca_prefix}cyt_' path_panel='panel_${extr_dir[$indx]}.xlsx' path_cytokines_cutoffs='panel_${extr_dir[$indx]}_cytokines.xlsx'" $RCODE/06_cytokines.R $ROUT/06_cytokines.Rout
  tail $ROUT/06_cytokines.Rout
  
  ## based on cytokines2.xlsx
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' cytokines_prefix='${pca_prefix}cyt2_' path_panel='panel_${extr_dir[$indx]}.xlsx' path_cytokines_cutoffs='panel_${extr_dir[$indx]}_cytokines2.xlsx'" $RCODE/06_cytokines.R $ROUT/06_cytokines.Rout
  tail $ROUT/06_cytokines.Rout
  

done
fi






