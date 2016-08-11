#!/bin/bash
## Define paths to software and reference files

RWD_MAIN=/Users/gosia/Dropbox/UZH/carsten_cytof
RCODE=$RWD_MAIN

## Define which analysis to re-run
pcascores=false
select_observables=false
flowsom=false
heatmaps=false
runtsne=false
plottsne=false
frequencies=false
cluster_merging=false
cluster_extracting=false
cytokines=false
fcs_saving=false
pd1=true

##############################################################################
# Analysis of CK_2016-06-23_01 data
##############################################################################

RWD=$RWD_MAIN/CK_2016-06-23_01
ROUT=$RWD/Rout
mkdir -p $ROUT
echo "$RWD"

### PCA scores
if ${pcascores}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' pcas_prefix='' path_panel='panel.xlsx'" $RCODE/01_pcascores.R $ROUT/01_pcascores.Rout
tail $ROUT/01_pcascores.Rout
fi

### Select observables for clustering
if ${select_observables}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' observ_prefix='pca1_' pca_score_cutoff=3 pca_skip_top=0" $RCODE/02_select_observables.R $ROUT/02_select_observables.Rout
tail $ROUT/02_select_observables.Rout
fi

### FlowSOM clustering
if ${flowsom}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' flowsom_prefix='pca1_cl20_' path_clustering_observables='pca1_clustering_observables.xls' nmetaclusts=20 rand_seed_consensus=123" $RCODE/02_flowsom.R $ROUT/02_flowsom.Rout
tail $ROUT/02_flowsom.Rout
fi

### Heatmaps
if ${heatmaps}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='pca1_cl20_' path_panel='panel.xlsx' path_clustering_observables='pca1_clustering_observables.xls' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls' path_pcascore='princompscore_by_sample.xls'" $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
tail $ROUT/02_heatmaps.Rout
fi

### Run tSNE
if ${runtsne}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsne_prefix='pca1_cl20_'  path_clustering_observables='pca1_clustering_observables.xls' tsne_pmin=1500" $RCODE/03_runtsne.R $ROUT/03_runtsne.Rout
tail $ROUT/03_runtsne.Rout
fi

### Plot tSNE
if ${plottsne}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='pca1_cl20_'  path_rtsne_out='pca1_cl20_rtsne_out.rda' path_rtsne_data='pca1_cl20_rtsne_data.xls' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
tail $ROUT/03_plottsne.Rout
fi


### Get cluster frequencies
if ${frequencies}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='pca1_cl20_' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
tail $ROUT/04_frequencies.Rout
fi


##########################################
# Analysis of CK_2016-06-23_01 mergingNEW
##########################################

### Cluster merging
if ${cluster_merging}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' merging_prefix='pca1_mergingNEW_' path_cluster_merging='cluster_mergingNEW.xlsx' path_clustering='pca1_cl20_clustering.xls'" $RCODE/02_cluster_merging.R $ROUT/02_cluster_merging.Rout
tail $ROUT/02_cluster_merging.Rout
fi

### Heatmaps
if ${heatmaps}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='pca1_mergingNEW_' path_panel='panel.xlsx' path_clustering_observables='pca1_clustering_observables.xls' path_clustering='pca1_mergingNEW_clustering.xls' path_clustering_labels='pca1_mergingNEW_clustering_labels.xls' path_pcascore='princompscore_by_sample.xls' " $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
tail $ROUT/02_heatmaps.Rout
fi

### Plot tSNE
if ${plottsne}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='pca1_mergingNEW_'  path_rtsne_out='pca1_cl20_rtsne_out.rda' path_rtsne_data='pca1_cl20_rtsne_data.xls' path_clustering='pca1_mergingNEW_clustering.xls' path_clustering_labels='pca1_mergingNEW_clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
tail $ROUT/03_plottsne.Rout
fi

### Get cluster frequencies
if ${frequencies}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='pca1_mergingNEW_' path_clustering='pca1_mergingNEW_clustering.xls' path_clustering_labels='pca1_mergingNEW_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
tail $ROUT/04_frequencies.Rout
fi

##########################################
# Analysis of CK_2016-06-23_01 mergingNEW2
##########################################

### Cluster merging
if ${cluster_merging}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' merging_prefix='pca1_mergingNEW2_' path_cluster_merging='cluster_mergingNEW2.xlsx' path_clustering='pca1_cl20_clustering.xls'" $RCODE/02_cluster_merging.R $ROUT/02_cluster_merging.Rout
tail $ROUT/02_cluster_merging.Rout
fi

### Heatmaps
if ${heatmaps}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='pca1_mergingNEW2_' path_panel='panel.xlsx' path_clustering_observables='pca1_clustering_observables.xls' path_clustering='pca1_mergingNEW2_clustering.xls' path_clustering_labels='pca1_mergingNEW2_clustering_labels.xls' path_pcascore='princompscore_by_sample.xls' " $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
tail $ROUT/02_heatmaps.Rout
fi

### Plot tSNE
if ${plottsne}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='pca1_mergingNEW2_'  path_rtsne_out='pca1_cl20_rtsne_out.rda' path_rtsne_data='pca1_cl20_rtsne_data.xls' path_clustering='pca1_mergingNEW2_clustering.xls' path_clustering_labels='pca1_mergingNEW2_clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
tail $ROUT/03_plottsne.Rout
fi

### Get cluster frequencies
if ${frequencies}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='pca1_mergingNEW2_' path_clustering='pca1_mergingNEW2_clustering.xls' path_clustering_labels='pca1_mergingNEW2_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
tail $ROUT/04_frequencies.Rout
fi


##########################################
# Analysis of CK_2016-06-23_01_CD4_mergingNEW2 and CK_2016-06-23_01_CD8_mergingNEW2
# using panel.xlsx
##########################################

data='CK_2016-06-23_01'
merging='mergingNEW2'
extr_name=('CD4' 'CD8')
extr_dir=('CD4' 'CD8')

if ${cluster_extracting}; then
for indx in 0 1
do
  
  echo "${data}_${extr_dir[$indx]}_${merging}"
  
  RWD=$RWD_MAIN/$data
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"

  ### Cluster extracting
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_clustering='pca1_${merging}_clustering.xls' path_clustering_labels='pca1_${merging}_clustering_labels.xls' extract_cluster='${extr_name[$indx]}' extract_dir='${RWD_MAIN}/${data}_${extr_dir[$indx]}_${merging}/010_cleanfcs'" $RCODE/02_cluster_extracting.R $ROUT/02_cluster_extracting.Rout
  tail $ROUT/02_cluster_extracting.Rout

done
fi

if ${pcascores}; then
for indx in 0 1
do
  
  echo "${data}_${extr_dir[$indx]}_${merging}"
  
  RWD=$RWD_MAIN/$data
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"
  
  ### Copy metadata and panel 
  cp ${RWD}/metadata.xlsx $RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}/
  cp ${RWD}/panel.xlsx $RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}/


  RWD=$RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"

  ### PCA scores
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' pcas_prefix='' path_panel='panel.xlsx'" $RCODE/01_pcascores.R $ROUT/01_pcascores.Rout
  tail $ROUT/01_pcascores.Rout

done
fi

##########################################
# Analysis of CK_2016-06-23_01_CD4_mergingNEW2 using panel_CD4.xlsx
# and CK_2016-06-23_01_CD8_mergingNEW2 using panel_CD8.xlsx
##########################################

data='CK_2016-06-23_01'
merging='mergingNEW2'
extr_name=('CD4' 'CD8')
extr_dir=('CD4' 'CD8')

if ${pcascores}; then
for indx in 0 1
do
  
  echo "${data}_${extr_dir[$indx]}_${merging}"
  
  RWD=$RWD_MAIN/$data
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"
  
  ### Copy metadata and panel 
  cp ${RWD}/metadata.xlsx $RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}/
  cp ${RWD}/panel_${extr_dir[$indx]}.xlsx $RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}/


  RWD=$RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"

  ### PCA scores
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' pcas_prefix='pnl${extr_dir[$indx]}_' path_panel='panel_${extr_dir[$indx]}.xlsx'" $RCODE/01_pcascores.R $ROUT/01_pcascores.Rout
  tail $ROUT/01_pcascores.Rout

done
fi


pca_cutoff=(2 2)

for indx in 0 1
do
  
  echo "${data}_${extr_dir[$indx]}_${merging}"
  
  RWD=$RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"
  
  pca_prefix="pnl${extr_dir[$indx]}_pca1_"
  
  ### Select observables for clustering
  if ${select_observables}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' observ_prefix='${pca_prefix}' pca_score_cutoff=${pca_cutoff[$indx]} pca_skip_top=0" $RCODE/02_select_observables.R $ROUT/02_select_observables.Rout
  tail $ROUT/02_select_observables.Rout
  fi

  ### FlowSOM clustering
  if ${flowsom}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' flowsom_prefix='${pca_prefix}cl20_' path_clustering_observables='${pca_prefix}clustering_observables.xls' nmetaclusts=20 rand_seed_consensus=1234" $RCODE/02_flowsom.R $ROUT/02_flowsom.Rout
  tail $ROUT/02_flowsom.Rout
  fi

  ### Heatmaps
  if ${heatmaps}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='${pca_prefix}cl20_' path_panel='panel_${extr_dir[$indx]}.xlsx' path_clustering_observables='${pca_prefix}clustering_observables.xls' path_clustering='${pca_prefix}cl20_clustering.xls' path_clustering_labels='${pca_prefix}cl20_clustering_labels.xls' path_pcascore='pnl${extr_dir[$indx]}_princompscore_by_sample.xls' " $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
  tail $ROUT/02_heatmaps.Rout
  fi

  ### Run tSNE
  if ${runtsne}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsne_prefix='${pca_prefix}cl20_'  path_clustering_observables='${pca_prefix}clustering_observables.xls' tsne_pmin=1500" $RCODE/03_runtsne.R $ROUT/03_runtsne.Rout
  tail $ROUT/03_runtsne.Rout
  fi

  ### Plot tSNE
  if ${plottsne}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='${pca_prefix}cl20_'  path_rtsne_out='${pca_prefix}cl20_rtsne_out.rda' path_rtsne_data='${pca_prefix}cl20_rtsne_data.xls' path_clustering='${pca_prefix}cl20_clustering.xls' path_clustering_labels='${pca_prefix}cl20_clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
  tail $ROUT/03_plottsne.Rout
  fi

  ### Get cluster frequencies
  if ${frequencies}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${pca_prefix}cl20_' path_clustering='${pca_prefix}cl20_clustering.xls' path_clustering_labels='${pca_prefix}cl20_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
  tail $ROUT/04_frequencies.Rout
  fi
  
done


##########################################
# Analysis of CK_2016-06-23_01_CD4_mergingNEW2 and CK_2016-06-23_01_CD8_mergingNEW2
# for cluster_merging_CD4_2.xlsx and cluster_merging_CD8_2.xlsx
##########################################

data='CK_2016-06-23_01'
merging='mergingNEW2'
extr_name=('CD4' 'CD8')
extr_dir=('CD4' 'CD8')


for indx in 0 1
do
  
  echo "${data}_${extr_dir[$indx]}_${merging}"
  
  RWD=$RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"
  
  pca_prefix="pnl${extr_dir[$indx]}_pca1_"
  merging_prefix="merging_${extr_dir[$indx]}_2"
  
  ### Cluster merging
  if ${cluster_merging}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' merging_prefix='${pca_prefix}${merging_prefix}_' path_cluster_merging='cluster_${merging_prefix}.xlsx' path_clustering='${pca_prefix}cl20_clustering.xls'" $RCODE/02_cluster_merging.R $ROUT/02_cluster_merging.Rout
  tail $ROUT/02_cluster_merging.Rout
  fi

  ### Heatmaps
  if ${heatmaps}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='${pca_prefix}${merging_prefix}_' path_panel='panel_${extr_dir[$indx]}.xlsx' path_clustering_observables='${pca_prefix}clustering_observables.xls' path_clustering='${pca_prefix}${merging_prefix}_clustering.xls' path_clustering_labels='${pca_prefix}${merging_prefix}_clustering_labels.xls' path_pcascore='pnl${extr_dir[$indx]}_princompscore_by_sample.xls' " $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
  tail $ROUT/02_heatmaps.Rout
  fi

  ### Plot tSNE
  if ${plottsne}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='${pca_prefix}${merging_prefix}_'  path_rtsne_out='${pca_prefix}cl20_rtsne_out.rda' path_rtsne_data='${pca_prefix}cl20_rtsne_data.xls' path_clustering='${pca_prefix}${merging_prefix}_clustering.xls' path_clustering_labels='${pca_prefix}${merging_prefix}_clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
  tail $ROUT/03_plottsne.Rout
  fi

  ### Get cluster frequencies
  if ${frequencies}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${pca_prefix}${merging_prefix}_' path_clustering='${pca_prefix}${merging_prefix}_clustering.xls' path_clustering_labels='${pca_prefix}${merging_prefix}_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
  tail $ROUT/04_frequencies.Rout
  fi

done





##########################################
# Analysis of CK_2016-06-23_01 data using panel2
##########################################

RWD=$RWD_MAIN/CK_2016-06-23_01
ROUT=$RWD/Rout
mkdir -p $ROUT
echo "$RWD"

### PCA scores
if ${pcascores}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' pcas_prefix='pnl2_' path_panel='panel2.xlsx'" $RCODE/01_pcascores.R $ROUT/01_pcascores.Rout
tail $ROUT/01_pcascores.Rout
fi

##############################################################################
# Analysis of CK_2016-06-23_02 data
##############################################################################


RWD=$RWD_MAIN/CK_2016-06-23_02
ROUT=$RWD/Rout
mkdir -p $ROUT
echo "$RWD"

### PCA scores
if ${pcascores}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' pcas_prefix='' path_panel='panel.xlsx'" $RCODE/01_pcascores.R $ROUT/01_pcascores.Rout
tail $ROUT/01_pcascores.Rout
fi

### Select observables for clustering
if ${select_observables}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' observ_prefix='pca1_' pca_score_cutoff=1 pca_skip_top=0" $RCODE/02_select_observables.R $ROUT/02_select_observables.Rout
tail $ROUT/02_select_observables.Rout
fi

### FlowSOM clustering
if ${flowsom}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' flowsom_prefix='pca1_cl20_' path_clustering_observables='pca1_clustering_observables.xls' nmetaclusts=20 rand_seed_consensus=1234" $RCODE/02_flowsom.R $ROUT/02_flowsom.Rout
tail $ROUT/02_flowsom.Rout
fi

### Heatmaps
if ${heatmaps}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='pca1_cl20_' path_panel='panel.xlsx' path_clustering_observables='pca1_clustering_observables.xls' path_clustering='pca1_cl20_clustering.xls'  path_clustering_labels='pca1_cl20_clustering_labels.xls' path_pcascore='princompscore_by_sample.xls' " $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
tail $ROUT/02_heatmaps.Rout
fi


### Run tSNE
if ${runtsne}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsne_prefix='pca1_cl20_'  path_clustering_observables='pca1_clustering_observables.xls' tsne_pmin=1500" $RCODE/03_runtsne.R $ROUT/03_runtsne.Rout
tail $ROUT/03_runtsne.Rout
fi

### Plot tSNE
if ${plottsne}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='pca1_cl20_'  path_rtsne_out='pca1_cl20_rtsne_out.rda' path_rtsne_data='pca1_cl20_rtsne_data.xls' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
tail $ROUT/03_plottsne.Rout
fi

### Get cluster frequencies
if ${frequencies}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='pca1_cl20_' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
tail $ROUT/04_frequencies.Rout
fi


##########################################
# Analysis of CK_2016-06-23_02 merging2
##########################################

### Cluster merging
if ${cluster_merging}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' merging_prefix='pca1_merging2_' path_cluster_merging='cluster_merging2.xlsx' path_clustering='pca1_cl20_clustering.xls'" $RCODE/02_cluster_merging.R $ROUT/02_cluster_merging.Rout
tail $ROUT/02_cluster_merging.Rout
fi

### Heatmaps
if ${heatmaps}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='pca1_merging2_' path_panel='panel.xlsx' path_clustering_observables='pca1_clustering_observables.xls' path_clustering='pca1_merging2_clustering.xls' path_clustering_labels='pca1_merging2_clustering_labels.xls' path_pcascore='princompscore_by_sample.xls' " $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
tail $ROUT/02_heatmaps.Rout
fi

### Plot tSNE
if ${plottsne}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='pca1_merging2_'  path_rtsne_out='pca1_cl20_rtsne_out.rda' path_rtsne_data='pca1_cl20_rtsne_data.xls' path_clustering='pca1_merging2_clustering.xls' path_clustering_labels='pca1_merging2_clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
tail $ROUT/03_plottsne.Rout
fi

### Get cluster frequencies
if ${frequencies}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='pca1_merging2_' path_clustering='pca1_merging2_clustering.xls' path_clustering_labels='pca1_merging2_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
tail $ROUT/04_frequencies.Rout
fi

##########################################
# Analysis of CK_2016-06-23_02_CD4_merging2 and CK_2016-06-23_02_CD8_merging2
# using panel.xlsx
##########################################

data='CK_2016-06-23_02'
merging='merging2'
extr_name=('CD4' 'CD8')
extr_dir=('CD4' 'CD8')

if ${cluster_extracting}; then
for indx in 0 1
do
  
  echo "${data}_${extr_dir[$indx]}_${merging}"
  
  RWD=$RWD_MAIN/$data
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"

  ### Cluster extracting
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_clustering='pca1_${merging}_clustering.xls' path_clustering_labels='pca1_${merging}_clustering_labels.xls' extract_cluster='${extr_name[$indx]}' extract_dir='${RWD_MAIN}/${data}_${extr_dir[$indx]}_${merging}/010_cleanfcs'" $RCODE/02_cluster_extracting.R $ROUT/02_cluster_extracting.Rout
  tail $ROUT/02_cluster_extracting.Rout

done
fi


if ${pcascores}; then
for indx in 0 1
do
  
  echo "${data}_${extr_dir[$indx]}_${merging}"
  
  RWD=$RWD_MAIN/$data
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"
  
  ### Copy metadata and panel 
  cp ${RWD}/metadata.xlsx $RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}/
  cp ${RWD}/panel.xlsx $RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}/


  RWD=$RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"

  ### PCA scores
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' pcas_prefix='' path_panel='panel.xlsx'" $RCODE/01_pcascores.R $ROUT/01_pcascores.Rout
  tail $ROUT/01_pcascores.Rout

done
fi

##########################################
# Analysis of CK_2016-06-23_02_CD4_merging2 using panel_CD4.xlsx
# and CK_2016-06-23_02_CD8_merging2 using panel_CD8.xlsx
##########################################

data='CK_2016-06-23_02'
merging='merging2'
extr_name=('CD4' 'CD8')
extr_dir=('CD4' 'CD8')

if ${pcascores}; then
for indx in 0 1
do
  
  echo "${data}_${extr_dir[$indx]}_${merging}"
  
  RWD=$RWD_MAIN/$data
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"
  
  ### Copy metadata and panel 
  cp ${RWD}/metadata.xlsx $RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}/
  cp ${RWD}/panel_${extr_dir[$indx]}.xlsx $RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}/


  RWD=$RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"

  ### PCA scores
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' pcas_prefix='pnl${extr_dir[$indx]}_' path_panel='panel_${extr_dir[$indx]}.xlsx'" $RCODE/01_pcascores.R $ROUT/01_pcascores.Rout
  tail $ROUT/01_pcascores.Rout

done
fi


if ${fcs_saving}; then
for indx in 0 1
do
  
  RWD=$RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"

  ### Save fcs files with arcsineh and arcsideh+01 normalized expression into 070_dumpfcs/
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' save_prefix='pnl${extr_dir[$indx]}_' path_panel='panel_${extr_dir[$indx]}.xlsx'" $RCODE/02_fcs_saving.R $ROUT/02_fcs_saving.Rout
  tail $ROUT/02_fcs_saving.Rout

done
fi



extr_dir=('CD4' 'CD8')
pca_cutoff=(0.96 1)

for indx in 0 1
do
  
  echo "${data}_${extr_dir[$indx]}_${merging}"
  
  RWD=$RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"
  
  pca_prefix="pnl${extr_dir[$indx]}_pca1_"
  
  ### Select observables for clustering
  if ${select_observables}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' observ_prefix='${pca_prefix}' pca_score_cutoff=${pca_cutoff[$indx]} pca_skip_top=0" $RCODE/02_select_observables.R $ROUT/02_select_observables.Rout
  tail $ROUT/02_select_observables.Rout
  fi

  ### FlowSOM clustering
  if ${flowsom}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' flowsom_prefix='${pca_prefix}cl20_' path_clustering_observables='${pca_prefix}clustering_observables.xls' nmetaclusts=20 rand_seed_consensus=1234" $RCODE/02_flowsom.R $ROUT/02_flowsom.Rout
  tail $ROUT/02_flowsom.Rout
  fi

  ### Heatmaps
  if ${heatmaps}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='${pca_prefix}cl20_' path_panel='panel_${extr_dir[$indx]}.xlsx' path_clustering_observables='${pca_prefix}clustering_observables.xls' path_clustering='${pca_prefix}cl20_clustering.xls' path_clustering_labels='${pca_prefix}cl20_clustering_labels.xls' path_pcascore='pnl${extr_dir[$indx]}_princompscore_by_sample.xls' " $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
  tail $ROUT/02_heatmaps.Rout
  fi

  ### Run tSNE
  if ${runtsne}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsne_prefix='${pca_prefix}cl20_'  path_clustering_observables='${pca_prefix}clustering_observables.xls' tsne_pmin=1500" $RCODE/03_runtsne.R $ROUT/03_runtsne.Rout
  tail $ROUT/03_runtsne.Rout
  fi

  ### Plot tSNE
  if ${plottsne}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='${pca_prefix}cl20_'  path_rtsne_out='${pca_prefix}cl20_rtsne_out.rda' path_rtsne_data='${pca_prefix}cl20_rtsne_data.xls' path_clustering='${pca_prefix}cl20_clustering.xls' path_clustering_labels='${pca_prefix}cl20_clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
  tail $ROUT/03_plottsne.Rout
  fi

  ### Get cluster frequencies
  if ${frequencies}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${pca_prefix}cl20_' path_clustering='${pca_prefix}cl20_clustering.xls' path_clustering_labels='${pca_prefix}cl20_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
  tail $ROUT/04_frequencies.Rout
  fi
  
done


##########################################
# Analysis of CK_2016-06-23_02_CD4_merging2 and CK_2016-06-23_02_CD8_merging2
# for cluster_merging_CD4.xlsx and cluster_merging_CD8.xlsx
##########################################

data='CK_2016-06-23_02'
merging='merging2'
extr_name=('CD4' 'CD8')
extr_dir=('CD4' 'CD8')


for indx in 0 1
do
  
  echo "${data}_${extr_dir[$indx]}_${merging}"
  
  RWD=$RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"
  
  pca_prefix="pnl${extr_dir[$indx]}_pca1_"
  merging_prefix="merging_${extr_dir[$indx]}"
  
  ### Cluster merging
  if ${cluster_merging}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' merging_prefix='${pca_prefix}${merging_prefix}_' path_cluster_merging='cluster_${merging_prefix}.xlsx' path_clustering='${pca_prefix}cl20_clustering.xls'" $RCODE/02_cluster_merging.R $ROUT/02_cluster_merging.Rout
  tail $ROUT/02_cluster_merging.Rout
  fi

  ### Heatmaps
  if ${heatmaps}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='${pca_prefix}${merging_prefix}_' path_panel='panel_${extr_dir[$indx]}.xlsx' path_clustering_observables='${pca_prefix}clustering_observables.xls' path_clustering='${pca_prefix}${merging_prefix}_clustering.xls' path_clustering_labels='${pca_prefix}${merging_prefix}_clustering_labels.xls' path_pcascore='pnl${extr_dir[$indx]}_princompscore_by_sample.xls' " $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
  tail $ROUT/02_heatmaps.Rout
  fi

  ### Plot tSNE
  if ${plottsne}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='${pca_prefix}${merging_prefix}_'  path_rtsne_out='${pca_prefix}cl20_rtsne_out.rda' path_rtsne_data='${pca_prefix}cl20_rtsne_data.xls' path_clustering='${pca_prefix}${merging_prefix}_clustering.xls' path_clustering_labels='${pca_prefix}${merging_prefix}_clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
  tail $ROUT/03_plottsne.Rout
  fi

  ### Get cluster frequencies
  if ${frequencies}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${pca_prefix}${merging_prefix}_' path_clustering='${pca_prefix}${merging_prefix}_clustering.xls' path_clustering_labels='${pca_prefix}${merging_prefix}_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
  tail $ROUT/04_frequencies.Rout
  fi


done


clusters2analyse="c('CM','EM','TE')"
cluster_name="Tmem"

### Analysis of positive-negative (cytokine) markers
nmetaclusts=40

## Based on cluster_merging_CD4.xlsx
if ${cytokines}; then
for indx in 0
do

  RWD=$RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"
  
  pca_prefix="pnl${extr_dir[$indx]}_pca1_" # 'pnlCD4_pca1_'
  merging_prefix="merging_${extr_dir[$indx]}" # 'merging_CD4'

  ## based on cytokines_CM_RAW.xlsx
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' cytokines_prefix='${pca_prefix}${merging_prefix}_cyt${cluster_name}_' path_panel='panel_${extr_dir[$indx]}.xlsx' path_cytokines_cutoffs='panel_${extr_dir[$indx]}_cytokines_CM_RAW.xlsx' path_clustering='${pca_prefix}${merging_prefix}_clustering.xls' path_clustering_labels='${pca_prefix}${merging_prefix}_clustering_labels.xls' clusters2analyse=${clusters2analyse} cutoff_colname='positive_cutoff_raw' data2analyse='raw' cytokines_suffix='_${nmetaclusts}cl_raw' nmetaclusts=${nmetaclusts}" $RCODE/06_cytokines.R $ROUT/06_cytokines.Rout
 tail $ROUT/06_cytokines.Rout
  

done
fi


nmetaclusts=20

## Based on cluster_merging_CD8_2.xlsx
if ${cytokines}; then
for indx in 1
do

  RWD=$RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"
  
  pca_prefix="pnl${extr_dir[$indx]}_pca1_" # 'pnlCD8_pca1_'
  merging_prefix="merging_${extr_dir[$indx]}_2" # 'merging_CD8_2'
  
  ## based on cytokines_CM_RAW.xlsx
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' cytokines_prefix='${pca_prefix}${merging_prefix}_cyt${cluster_name}_' path_panel='panel_${extr_dir[$indx]}.xlsx' path_cytokines_cutoffs='panel_${extr_dir[$indx]}_cytokines_CM_RAW.xlsx' path_clustering='${pca_prefix}${merging_prefix}_clustering.xls' path_clustering_labels='${pca_prefix}${merging_prefix}_clustering_labels.xls' clusters2analyse=${clusters2analyse} cutoff_colname='positive_cutoff_raw' data2analyse='raw' cytokines_suffix='_${nmetaclusts}cl_raw' nmetaclusts=${nmetaclusts}" $RCODE/06_cytokines.R $ROUT/06_cytokines.Rout
 tail $ROUT/06_cytokines.Rout
  

done
fi



### Analysis of PD-1
nmetaclusts=20

## Based on cluster_merging_CD4.xlsx
if ${pd1}; then
for indx in 0
do

  RWD=$RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"
  
  pca_prefix="pnl${extr_dir[$indx]}_pca1_" # 'pnlCD4_pca1_'
  merging_prefix="merging_${extr_dir[$indx]}" # 'merging_CD4'

  ## based on cytokines_CM_RAW.xlsx
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' cytokines_prefix='${pca_prefix}${merging_prefix}_cyt${cluster_name}_' path_panel='panel_${extr_dir[$indx]}.xlsx' path_cytokines_cutoffs='panel_${extr_dir[$indx]}_cytokines_CM_RAW.xlsx' path_clustering='${pca_prefix}${merging_prefix}_clustering.xls' path_clustering_labels='${pca_prefix}${merging_prefix}_clustering_labels.xls' clusters2analyse=${clusters2analyse} cutoff_colname='positive_cutoff_raw' data2analyse='raw' cytokines_suffix='_${nmetaclusts}cl_raw' nmetaclusts=${nmetaclusts}" $RCODE/06_cytokines.R $ROUT/06_cytokines.Rout
 tail $ROUT/06_cytokines.Rout
  

done
fi


## Based on cluster_merging_CD8_2.xlsx
if ${pd1}; then
for indx in 1
do

  RWD=$RWD_MAIN/${data}_${extr_dir[$indx]}_${merging}
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"
  
  pca_prefix="pnl${extr_dir[$indx]}_pca1_" # 'pnlCD8_pca1_'
  merging_prefix="merging_${extr_dir[$indx]}_2" # 'merging_CD8_2'
  
  ## based on cytokines_CM_RAW.xlsx
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' cytokines_prefix='${pca_prefix}${merging_prefix}_cyt${cluster_name}_' path_panel='panel_${extr_dir[$indx]}.xlsx' path_cytokines_cutoffs='panel_${extr_dir[$indx]}_cytokines_CM_RAW.xlsx' path_clustering='${pca_prefix}${merging_prefix}_clustering.xls' path_clustering_labels='${pca_prefix}${merging_prefix}_clustering_labels.xls' clusters2analyse=${clusters2analyse} cutoff_colname='positive_cutoff_raw' data2analyse='raw' cytokines_suffix='_${nmetaclusts}cl_raw' nmetaclusts=${nmetaclusts}" $RCODE/06_cytokines.R $ROUT/06_cytokines.Rout
 tail $ROUT/06_cytokines.Rout
  

done
fi

##############################################################################
# Analysis of CK_2016-06-23_03 data
##############################################################################

RWD=$RWD_MAIN/CK_2016-06-23_03
ROUT=$RWD/Rout
mkdir -p $ROUT
echo "$RWD"

### PCA scores
if ${pcascores}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' pcas_prefix='' path_panel='panel.xlsx'" $RCODE/01_pcascores.R $ROUT/01_pcascores.Rout
tail $ROUT/01_pcascores.Rout
fi

### Select observables for clustering
if ${select_observables}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' observ_prefix='pca1_' pca_score_cutoff=0.9 pca_skip_top=0" $RCODE/02_select_observables.R $ROUT/02_select_observables.Rout
tail $ROUT/02_select_observables.Rout
fi

### FlowSOM clustering
if ${flowsom}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' flowsom_prefix='pca1_cl20_' path_clustering_observables='pca1_clustering_observables.xls' nmetaclusts=20 rand_seed_consensus=1234" $RCODE/02_flowsom.R $ROUT/02_flowsom.Rout
tail $ROUT/02_flowsom.Rout
fi

### Heatmaps
if ${heatmaps}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='pca1_cl20_' path_panel='panel.xlsx' path_clustering_observables='pca1_clustering_observables.xls' path_clustering='pca1_cl20_clustering.xls'  path_clustering_labels='pca1_cl20_clustering_labels.xls' path_pcascore='princompscore_by_sample.xls' " $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
tail $ROUT/02_heatmaps.Rout
fi

### Run tSNE
if ${runtsne}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsne_prefix='pca1_cl20_'  path_clustering_observables='pca1_clustering_observables.xls' tsne_pmin=1500" $RCODE/03_runtsne.R $ROUT/03_runtsne.Rout
tail $ROUT/03_runtsne.Rout
fi

### Plot tSNE
if ${plottsne}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='pca1_cl20_'  path_rtsne_out='pca1_cl20_rtsne_out.rda' path_rtsne_data='pca1_cl20_rtsne_data.xls' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
tail $ROUT/03_plottsne.Rout
fi

### Get cluster frequencies
if ${frequencies}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='pca1_cl20_' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
tail $ROUT/04_frequencies.Rout
fi

####################################
### PCA2 approach 

### Select observables for clustering
if ${select_observables}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' observ_prefix='pca2_' pca_score_cutoff=0.9 pca_skip_top=2" $RCODE/02_select_observables.R $ROUT/02_select_observables.Rout
tail $ROUT/02_select_observables.Rout
fi

### FlowSOM clustering
if ${flowsom}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' flowsom_prefix='pca2_cl20_' path_clustering_observables='pca2_clustering_observables.xls' nmetaclusts=20 rand_seed_consensus=123" $RCODE/02_flowsom.R $ROUT/02_flowsom.Rout
tail $ROUT/02_flowsom.Rout
fi

### Heatmaps
if ${heatmaps}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='pca2_cl20_' path_panel='panel.xlsx' path_clustering_observables='pca2_clustering_observables.xls' path_clustering='pca2_cl20_clustering.xls'  path_clustering_labels='pca2_cl20_clustering_labels.xls' path_pcascore='princompscore_by_sample.xls' " $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
tail $ROUT/02_heatmaps.Rout
fi

### Run tSNE
if ${runtsne}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsne_prefix='pca2_cl20_'  path_clustering_observables='pca2_clustering_observables.xls' tsne_pmin=1500" $RCODE/03_runtsne.R $ROUT/03_runtsne.Rout
tail $ROUT/03_runtsne.Rout
fi

### Plot tSNE
if ${plottsne}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='pca2_cl20_'  path_rtsne_out='pca2_cl20_rtsne_out.rda' path_rtsne_data='pca2_cl20_rtsne_data.xls' path_clustering='pca2_cl20_clustering.xls' path_clustering_labels='pca2_cl20_clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
tail $ROUT/03_plottsne.Rout
fi



##########################################
# Analysis of CK_2016-06-23_03 merging2
##########################################

### Cluster merging
if ${cluster_merging}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' merging_prefix='pca1_merging2_' path_cluster_merging='cluster_merging2.xlsx' path_clustering='pca1_cl20_clustering.xls'" $RCODE/02_cluster_merging.R $ROUT/02_cluster_merging.Rout
tail $ROUT/02_cluster_merging.Rout
fi

### Heatmaps
if ${heatmaps}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='pca1_merging2_' path_panel='panel.xlsx' path_clustering_observables='pca1_clustering_observables.xls' path_clustering='pca1_merging2_clustering.xls' path_clustering_labels='pca1_merging2_clustering_labels.xls' path_pcascore='princompscore_by_sample.xls' " $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
tail $ROUT/02_heatmaps.Rout
fi

### Plot tSNE
if ${plottsne}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='pca1_merging2_'  path_rtsne_out='pca1_cl20_rtsne_out.rda' path_rtsne_data='pca1_cl20_rtsne_data.xls' path_clustering='pca1_merging2_clustering.xls' path_clustering_labels='pca1_merging2_clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
tail $ROUT/03_plottsne.Rout
fi

### Get cluster frequencies
if ${frequencies}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='pca1_merging2_' path_clustering='pca1_merging2_clustering.xls' path_clustering_labels='pca1_merging2_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
tail $ROUT/04_frequencies.Rout
fi

##########################################
# Analysis of CK_2016-06-23_03 merging3
##########################################

### Cluster merging
if ${cluster_merging}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' merging_prefix='pca1_merging3_' path_cluster_merging='cluster_merging3.xlsx' path_clustering='pca1_cl20_clustering.xls'" $RCODE/02_cluster_merging.R $ROUT/02_cluster_merging.Rout
tail $ROUT/02_cluster_merging.Rout
fi

### Heatmaps
if ${heatmaps}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='pca1_merging3_' path_panel='panel.xlsx' path_clustering_observables='pca1_clustering_observables.xls' path_clustering='pca1_merging3_clustering.xls' path_clustering_labels='pca1_merging3_clustering_labels.xls' path_pcascore='princompscore_by_sample.xls' " $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
tail $ROUT/02_heatmaps.Rout
fi

### Plot tSNE
if ${plottsne}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='pca1_merging3_'  path_rtsne_out='pca1_cl20_rtsne_out.rda' path_rtsne_data='pca1_cl20_rtsne_data.xls' path_clustering='pca1_merging3_clustering.xls' path_clustering_labels='pca1_merging3_clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
tail $ROUT/03_plottsne.Rout
fi

### Get cluster frequencies
if ${frequencies}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='pca1_merging3_' path_clustering='pca1_merging3_clustering.xls' path_clustering_labels='pca1_merging3_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
tail $ROUT/04_frequencies.Rout
fi

##############################################################################
# Analysis of CK_2016-06-29_01 data
##############################################################################

RWD=$RWD_MAIN/CK_2016-06-29_01
ROUT=$RWD/Rout
mkdir -p $ROUT
echo "$RWD"

### PCA scores
if ${pcascores}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' pcas_prefix='' path_panel='panel.xlsx'" $RCODE/01_pcascores.R $ROUT/01_pcascores.Rout
tail $ROUT/01_pcascores.Rout
fi

### Select observables for clustering
if ${select_observables}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' observ_prefix='pca1_' pca_score_cutoff=3.5 pca_skip_top=0" $RCODE/02_select_observables.R $ROUT/02_select_observables.Rout
tail $ROUT/02_select_observables.Rout
fi

### FlowSOM clustering
if ${flowsom}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' flowsom_prefix='pca1_cl20_' path_clustering_observables='pca1_clustering_observables.xls' nmetaclusts=20 rand_seed_consensus=1234" $RCODE/02_flowsom.R $ROUT/02_flowsom.Rout
tail $ROUT/02_flowsom.Rout
fi

### Heatmaps
if ${heatmaps}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='pca1_cl20_' path_panel='panel.xlsx' path_clustering_observables='pca1_clustering_observables.xls' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls' path_pcascore='princompscore_by_sample.xls' " $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
tail $ROUT/02_heatmaps.Rout
fi

### Run tSNE
if ${runtsne}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsne_prefix='pca1_cl20_'  path_clustering_observables='pca1_clustering_observables.xls' tsne_pmin=1500" $RCODE/03_runtsne.R $ROUT/03_runtsne.Rout
tail $ROUT/03_runtsne.Rout
fi

### Plot tSNE
if ${plottsne}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='pca1_cl20_'  path_rtsne_out='pca1_cl20_rtsne_out.rda' path_rtsne_data='pca1_cl20_rtsne_data.xls' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
tail $ROUT/03_plottsne.Rout
fi

### Get cluster frequencies
if ${frequencies}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='pca1_cl20_' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
tail $ROUT/04_frequencies.Rout
fi

##########################################
# Analysis of CK_2016-06-29_01 merging2
##########################################

### Cluster merging
if ${cluster_merging}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' merging_prefix='pca1_merging2_' path_cluster_merging='cluster_merging2.xlsx' path_clustering='pca1_cl20_clustering.xls'" $RCODE/02_cluster_merging.R $ROUT/02_cluster_merging.Rout
tail $ROUT/02_cluster_merging.Rout
fi

### Heatmaps
if ${heatmaps}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='pca1_merging2_' path_panel='panel.xlsx' path_clustering_observables='pca1_clustering_observables.xls' path_clustering='pca1_merging2_clustering.xls' path_clustering_labels='pca1_merging2_clustering_labels.xls' path_pcascore='princompscore_by_sample.xls' " $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
tail $ROUT/02_heatmaps.Rout
fi

### Plot tSNE
if ${plottsne}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='pca1_merging2_'  path_rtsne_out='pca1_cl20_rtsne_out.rda' path_rtsne_data='pca1_cl20_rtsne_data.xls' path_clustering='pca1_merging2_clustering.xls' path_clustering_labels='pca1_merging2_clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
tail $ROUT/03_plottsne.Rout
fi

### Get cluster frequencies
if ${frequencies}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='pca1_merging2_' path_clustering='pca1_merging2_clustering.xls' path_clustering_labels='pca1_merging2_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
tail $ROUT/04_frequencies.Rout
fi

##########################################
# Analysis of CK_2016-06-29_01_CD4_merging2 and CK_2016-06-29_01_CD8_merging2
# using panel.xlsx
##########################################

data='CK_2016-06-29_01'
merging='merging2'
extr_name=('CD4' 'CD8')
extr_dir=('CD4' 'CD8')

##########################################
# Analysis of CK_2016-06-29_01_CD4_merging2 using panel_CD4.xlsx
# and CK_2016-06-29_01_CD8_merging2 using panel_CD8.xlsx
##########################################


##############################################################################
# Analysis of CK_2016-06-29_02 data
##############################################################################

RWD=$RWD_MAIN/CK_2016-06-29_02
ROUT=$RWD/Rout
mkdir -p $ROUT
echo "$RWD"

### PCA scores
if ${pcascores}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' pcas_prefix='' path_panel='panel.xlsx'" $RCODE/01_pcascores.R $ROUT/01_pcascores.Rout
tail $ROUT/01_pcascores.Rout
fi

### Select observables for clustering
if ${select_observables}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' observ_prefix='pca1_' pca_score_cutoff=1 pca_skip_top=0" $RCODE/02_select_observables.R $ROUT/02_select_observables.Rout
tail $ROUT/02_select_observables.Rout
fi

### FlowSOM clustering
if ${flowsom}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' flowsom_prefix='pca1_cl20_' path_clustering_observables='pca1_clustering_observables.xls' nmetaclusts=20 rand_seed_consensus=123" $RCODE/02_flowsom.R $ROUT/02_flowsom.Rout
tail $ROUT/02_flowsom.Rout
fi

### Heatmaps
if ${heatmaps}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='pca1_cl20_' path_panel='panel.xlsx' path_clustering_observables='pca1_clustering_observables.xls' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls' path_pcascore='princompscore_by_sample.xls' " $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
tail $ROUT/02_heatmaps.Rout
fi

### Run tSNE
if ${runtsne}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsne_prefix='pca1_cl20_'  path_clustering_observables='pca1_clustering_observables.xls' tsne_pmin=1500" $RCODE/03_runtsne.R $ROUT/03_runtsne.Rout
tail $ROUT/03_runtsne.Rout
fi

### Plot tSNE
if ${plottsne}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='pca1_cl20_'  path_rtsne_out='pca1_cl20_rtsne_out.rda' path_rtsne_data='pca1_cl20_rtsne_data.xls' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
tail $ROUT/03_plottsne.Rout
fi

### Get cluster frequencies
if ${frequencies}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='pca1_cl20_' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
tail $ROUT/04_frequencies.Rout
fi

##########################################
# Analysis of CK_2016-06-29_02 merging
##########################################

### Cluster merging
if ${cluster_merging}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' merging_prefix='pca1_merging_' path_cluster_merging='cluster_merging.xlsx' path_clustering='pca1_cl20_clustering.xls'" $RCODE/02_cluster_merging.R $ROUT/02_cluster_merging.Rout
tail $ROUT/02_cluster_merging.Rout
fi

### Heatmaps
if ${heatmaps}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='pca1_merging_' path_panel='panel.xlsx' path_clustering_observables='pca1_clustering_observables.xls' path_clustering='pca1_merging_clustering.xls' path_clustering_labels='pca1_merging_clustering_labels.xls' path_pcascore='princompscore_by_sample.xls' " $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
tail $ROUT/02_heatmaps.Rout
fi

### Plot tSNE
if ${plottsne}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='pca1_merging_'  path_rtsne_out='pca1_cl20_rtsne_out.rda' path_rtsne_data='pca1_cl20_rtsne_data.xls' path_clustering='pca1_merging_clustering.xls' path_clustering_labels='pca1_merging_clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
tail $ROUT/03_plottsne.Rout
fi

### Get cluster frequencies
if ${frequencies}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='pca1_merging_' path_clustering='pca1_merging_clustering.xls' path_clustering_labels='pca1_merging_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
tail $ROUT/04_frequencies.Rout
fi

##########################################
# Analysis of CK_2016-06-29_02_CD4_merging and CK_2016-06-29_02_CD8_merging
# using panel.xlsx
##########################################

data='CK_2016-06-29_02'
merging='merging'
extr_name=('CD4' 'CD8')
extr_dir=('CD4' 'CD8')


##########################################
# Analysis of CK_2016-06-29_02_CD4_merging using using panel_CD4.xlsx
# and CK_2016-06-29_02_CD8_merging using panel_CD8.xlsx
##########################################



##############################################################################
# Analysis of CK_2016-06-29_03 data
##############################################################################


RWD=$RWD_MAIN/CK_2016-06-29_03
ROUT=$RWD/Rout
mkdir -p $ROUT
echo "$RWD"

### PCA scores
if ${pcascores}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' pcas_prefix='' path_panel='panel.xlsx'" $RCODE/01_pcascores.R $ROUT/01_pcascores.Rout
tail $ROUT/01_pcascores.Rout
fi

### Select observables for clustering
if ${select_observables}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' observ_prefix='pca1_' pca_score_cutoff=0.9 pca_skip_top=0" $RCODE/02_select_observables.R $ROUT/02_select_observables.Rout
tail $ROUT/02_select_observables.Rout
fi

### FlowSOM clustering
if ${flowsom}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' flowsom_prefix='pca1_cl20_' path_clustering_observables='pca1_clustering_observables.xls' nmetaclusts=20 rand_seed_consensus=1234" $RCODE/02_flowsom.R $ROUT/02_flowsom.Rout
tail $ROUT/02_flowsom.Rout
fi

### Heatmaps
if ${heatmaps}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='pca1_cl20_' path_panel='panel.xlsx' path_clustering_observables='pca1_clustering_observables.xls' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls' path_pcascore='princompscore_by_sample.xls' " $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
tail $ROUT/02_heatmaps.Rout
fi

### Run tSNE
if ${runtsne}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsne_prefix='pca1_cl20_'  path_clustering_observables='pca1_clustering_observables.xls' tsne_pmin=1500" $RCODE/03_runtsne.R $ROUT/03_runtsne.Rout
tail $ROUT/03_runtsne.Rout
fi

### Plot tSNE
if ${plottsne}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='pca1_cl20_'  path_rtsne_out='pca1_cl20_rtsne_out.rda' path_rtsne_data='pca1_cl20_rtsne_data.xls' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
tail $ROUT/03_plottsne.Rout
fi

### Get cluster frequencies
if ${frequencies}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='pca1_cl20_' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
tail $ROUT/04_frequencies.Rout
fi

#############################################
# Analysis of CK_2016-06-29_03 data 
# using observables selected in pca1 CK_2016-06-23_03
#############################################


### FlowSOM clustering
if ${flowsom}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' flowsom_prefix='pca1v23_cl20_' path_clustering_observables='${RWD_MAIN}/CK_2016-06-23_03/030_heatmaps/pca1_clustering_observables.xls' nmetaclusts=20 rand_seed_consensus=1234" $RCODE/02_flowsom.R $ROUT/02_flowsom.Rout
tail $ROUT/02_flowsom.Rout
fi

### Heatmaps
if ${heatmaps}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='pca1v23_cl20_' path_panel='panel.xlsx' path_clustering_observables='${RWD_MAIN}/CK_2016-06-23_03/030_heatmaps/pca1_clustering_observables.xls' path_clustering='pca1v23_cl20_clustering.xls' path_clustering_labels='pca1v23_cl20_clustering_labels.xls' path_pcascore='${RWD_MAIN}/CK_2016-06-23_03/020_pcascores/princompscore_by_sample.xls' " $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
tail ${ROUT}/02_heatmaps.Rout
fi

### Run tSNE
if ${runtsne}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsne_prefix='pca1v23_cl20_'  path_clustering_observables='${RWD_MAIN}/CK_2016-06-23_03/030_heatmaps/pca1_clustering_observables.xls' tsne_pmin=1500" $RCODE/03_runtsne.R $ROUT/03_runtsne.Rout
tail $ROUT/03_runtsne.Rout
fi

### Plot tSNE
if ${plottsne}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='pca1v23_cl20_'  path_rtsne_out='pca1v23_cl20_rtsne_out.rda' path_rtsne_data='pca1v23_cl20_rtsne_data.xls' path_clustering='pca1v23_cl20_clustering.xls' path_clustering_labels='pca1v23_cl20_clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
tail $ROUT/03_plottsne.Rout
fi

### Get cluster frequencies
if ${frequencies}; then
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='pca1v23_cl20_' path_clustering='pca1v23_cl20_clustering.xls' path_clustering_labels='pca1v23_cl20_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
tail $ROUT/04_frequencies.Rout
fi





























