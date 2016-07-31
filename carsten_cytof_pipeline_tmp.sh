#!/bin/bash
## Define paths to software and reference files

RWD_MAIN=/Users/gosia/Dropbox/UZH/carsten_cytof
RCODE=$RWD_MAIN


##########################################
# Analysis of CK_2016-06-23_02_CD4 data from merging2
# and CK_2016-06-23_02_CD8 data from merging2
##########################################

data='CK_2016-06-23_02'
merging='merging2'
extr_name=('CD4' 'CD8')
extr_dir=('CD4' 'CD8')

for indx in 0 1
do
  
  echo "${extr_name[$indx]}"
  
  RWD=$RWD_MAIN/$data
  ROUT=$RWD/Rout
  mkdir -p $ROUT

  ### Cluster extracting
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_clustering='pca1_${merging}_clustering.xls' path_clustering_labels='pca1_${merging}_clustering_labels.xls' extract_cluster='${extr_name[$indx]}' extract_dir='${RWD_MAIN}/${data}_${extr_dir[$indx]}_${merging}/010_cleanfcs'" $RCODE/02_cluster_extracting.R $ROUT/02_cluster_extracting.Rout
  tail $ROUT/02_cluster_extracting.Rout

  ### Copy metadata and panel 
  cp ${RWD}/metadata.xlsx $RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}/
  cp ${RWD}/panel.xlsx $RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}/


  RWD=$RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}
  ROUT=$RWD/Rout
  mkdir -p $ROUT

  ### PCA scores
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD'" $RCODE/01_pcascores.R $ROUT/01_pcascores.Rout
  tail $ROUT/01_pcascores.Rout

done


extr_dir=('CD4' 'CD8')
pca_cutoff=(1 1)

for indx in 0 1
do
  
  echo "${data}_${extr_dir[$indx]}_${merging}"
  
  RWD=$RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  
  ### Select observables for clustering
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' observ_prefix='pca1_' pca_score_cutoff=${pca_cutoff[$indx]} pca_skip_top=0" $RCODE/02_select_observables.R $ROUT/02_select_observables.Rout
  tail $ROUT/02_select_observables.Rout


  ### FlowSOM clustering
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' flowsom_prefix='pca1_cl20_' path_clustering_observables='pca1_clustering_observables.xls'" $RCODE/02_flowsom.R $ROUT/02_flowsom.Rout
  tail $ROUT/02_flowsom.Rout


  ### Heatmaps
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='pca1_cl20_' path_clustering_observables='pca1_clustering_observables.xls' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls' path_pcascore='princompscore_by_sample.xls' xspace=2" $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
  tail $ROUT/02_heatmaps.Rout


  ### Run tSNE
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsne_prefix='pca1_cl20_'  path_clustering_observables='pca1_clustering_observables.xls' tsne_pmin=1500" $RCODE/03_runtsne.R $ROUT/03_runtsne.Rout
  tail $ROUT/03_runtsne.Rout


  ### Plot tSNE
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='pca1_cl20_'  path_rtsne_out='pca1_cl20_rtsne_out.rda' path_rtsne_data='pca1_cl20_rtsne_data.xls' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
  tail $ROUT/03_plottsne.Rout


  ### Get cluster frequencies
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='pca1_cl20_' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
  tail $ROUT/04_frequencies.Rout

done



##########################################
# Analysis of CK_2016-06-29_01_CD4 data from merging2
# and CK_2016-06-29_01_CD8 data from merging2
##########################################

data='CK_2016-06-29_01'
merging='merging2'
extr_name=('CD4' 'CD8')
extr_dir=('CD4' 'CD8')

for indx in 0 1
do
  
  echo "${extr_name[$indx]}"
  
  RWD=$RWD_MAIN/$data
  ROUT=$RWD/Rout
  mkdir -p $ROUT

  ### Cluster extracting
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_clustering='pca1_${merging}_clustering.xls' path_clustering_labels='pca1_${merging}_clustering_labels.xls' extract_cluster='${extr_name[$indx]}' extract_dir='${RWD_MAIN}/${data}_${extr_dir[$indx]}_${merging}/010_cleanfcs'" $RCODE/02_cluster_extracting.R $ROUT/02_cluster_extracting.Rout
  tail $ROUT/02_cluster_extracting.Rout

  ### Copy metadata and panel 
  cp ${RWD}/metadata.xlsx $RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}/
  cp ${RWD}/panel.xlsx $RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}/


  RWD=$RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}
  ROUT=$RWD/Rout
  mkdir -p $ROUT

  ### PCA scores
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD'" $RCODE/01_pcascores.R $ROUT/01_pcascores.Rout
  tail $ROUT/01_pcascores.Rout

done


extr_dir=('CD4' 'CD8')
pca_cutoff=(1.84 2)

for indx in 0 1
do
  
  echo "${data}_${extr_dir[$indx]}_${merging}"
  
  RWD=$RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  
  ### Select observables for clustering
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' observ_prefix='pca1_' pca_score_cutoff=${pca_cutoff[$indx]} pca_skip_top=0" $RCODE/02_select_observables.R $ROUT/02_select_observables.Rout
  tail $ROUT/02_select_observables.Rout


  ### FlowSOM clustering
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' flowsom_prefix='pca1_cl20_' path_clustering_observables='pca1_clustering_observables.xls'" $RCODE/02_flowsom.R $ROUT/02_flowsom.Rout
  tail $ROUT/02_flowsom.Rout


  ### Heatmaps
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='pca1_cl20_' path_clustering_observables='pca1_clustering_observables.xls' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls' path_pcascore='princompscore_by_sample.xls' xspace=2" $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
  tail $ROUT/02_heatmaps.Rout


  ### Run tSNE
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsne_prefix='pca1_cl20_'  path_clustering_observables='pca1_clustering_observables.xls' tsne_pmin=1500" $RCODE/03_runtsne.R $ROUT/03_runtsne.Rout
  tail $ROUT/03_runtsne.Rout


  ### Plot tSNE
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='pca1_cl20_'  path_rtsne_out='pca1_cl20_rtsne_out.rda' path_rtsne_data='pca1_cl20_rtsne_data.xls' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
  tail $ROUT/03_plottsne.Rout


  ### Get cluster frequencies
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='pca1_cl20_' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
  tail $ROUT/04_frequencies.Rout

done




