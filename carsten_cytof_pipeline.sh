#!/bin/bash
## Define paths to software and reference files

RWD_MAIN=/Users/gosia/Dropbox/UZH/carsten_cytof
RCODE=$RWD_MAIN

##############################################################################
# Analysis of CK_2016-06-23_01 data
##############################################################################

RWD=$RWD_MAIN/CK_2016-06-23_01
ROUT=$RWD/Rout
mkdir -p $ROUT

### PCA scores
R CMD BATCH --no-save --no-restore "--args rwd='$RWD'" $RCODE/01_pcascores.R $ROUT/01_pcascores.Rout
tail $ROUT/01_pcascores.Rout


### FlowSOM clustering
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' flowsom_prefix='pca1_cl20_' pca_score_cutoff=3 pca_skip_top=0" $RCODE/02_flowsom.R $ROUT/02_flowsom.Rout
tail $ROUT/02_flowsom.Rout


### Heatmaps
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='pca1_cl20_' path_clustering='pca1_cl20_clustering.xls' path_clustering_observables='pca1_cl20_clustering_observables.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls' xspace=2" $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
tail $ROUT/02_heatmaps.Rout


### Run tSNE
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsne_prefix='pca1_cl20_'  path_clustering_observables='pca1_cl20_clustering_observables.xls' tsne_pmin=1500" $RCODE/03_runtsne.R $ROUT/03_runtsne.Rout
tail $ROUT/03_runtsne.Rout


### Plot tSNE
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='pca1_cl20_'  path_rtsne_out='pca1_cl20_rtsne_out.rda' path_rtsne_data='pca1_cl20_rtsne_data.xls' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls'  tsne_cmin=1000" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
tail $ROUT/03_plottsne.Rout


### Cluster frequencies
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='pca1_cl20_' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
tail $ROUT/04_frequencies.Rout



##########################################
# Analysis of CK_2016-06-23_01 mergingNEW
##########################################

### Cluster merging
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' merging_prefix='pca1_mergingNEW_' path_cluster_merging='cluster_mergingNEW.xlsx' path_clustering='pca1_cl20_clustering.xls'" $RCODE/02_cluster_merging.R $ROUT/02_cluster_merging.Rout
tail $ROUT/02_cluster_merging.Rout


### Heatmaps
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='pca1_mergingNEW_' path_clustering='pca1_mergingNEW_clustering.xls' path_clustering_observables='pca1_cl20_clustering_observables.xls' path_clustering_labels='pca1_mergingNEW_clustering_labels.xls' xspace=7" $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
tail $ROUT/02_heatmaps.Rout


### Plot tSNE
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='pca1_mergingNEW_'  path_rtsne_out='pca1_cl20_rtsne_out.rda' path_rtsne_data='pca1_cl20_rtsne_data.xls' path_clustering='pca1_mergingNEW_clustering.xls' path_clustering_labels='pca1_mergingNEW_clustering_labels.xls'  tsne_cmin=1000" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
tail $ROUT/03_plottsne.Rout


### Cluster frequencies
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='pca1_mergingNEW_' path_clustering='pca1_mergingNEW_clustering.xls' path_clustering_labels='pca1_mergingNEW_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
tail $ROUT/04_frequencies.Rout


##########################################
# Analysis of CK_2016-06-23_01 mergingNEW2
##########################################

### Cluster merging
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' merging_prefix='pca1_mergingNEW2_' path_cluster_merging='cluster_mergingNEW2.xlsx' path_clustering='pca1_cl20_clustering.xls'" $RCODE/02_cluster_merging.R $ROUT/02_cluster_merging.Rout
tail $ROUT/02_cluster_merging.Rout


### Heatmaps
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='pca1_mergingNEW2_' path_clustering='pca1_mergingNEW2_clustering.xls' path_clustering_observables='pca1_cl20_clustering_observables.xls' path_clustering_labels='pca1_mergingNEW2_clustering_labels.xls' xspace=7" $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
tail $ROUT/02_heatmaps.Rout


### Plot tSNE
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='pca1_mergingNEW2_'  path_rtsne_out='pca1_cl20_rtsne_out.rda' path_rtsne_data='pca1_cl20_rtsne_data.xls' path_clustering='pca1_mergingNEW2_clustering.xls' path_clustering_labels='pca1_mergingNEW2_clustering_labels.xls'  tsne_cmin=1000" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
tail $ROUT/03_plottsne.Rout


### Cluster frequencies
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='pca1_mergingNEW2_' path_clustering='pca1_mergingNEW2_clustering.xls' path_clustering_labels='pca1_mergingNEW2_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
tail $ROUT/04_frequencies.Rout



##########################################
# Analysis of CK_2016-06-23_01_CD4 data
##########################################






##############################################################################
# Analysis of CK_2016-06-23_02 data
##############################################################################


RWD=$RWD_MAIN/CK_2016-06-23_02
ROUT=$RWD/Rout
mkdir -p $ROUT

### PCA scores
R CMD BATCH --no-save --no-restore "--args rwd='$RWD'" $RCODE/01_pcascores.R $ROUT/01_pcascores.Rout
tail $ROUT/01_pcascores.Rout


### FlowSOM clustering
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' flowsom_prefix='pca1_cl20_' pca_score_cutoff=1" $RCODE/02_flowsom.R $ROUT/02_flowsom.Rout
tail $ROUT/02_flowsom.Rout


### Heatmaps
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='pca1_cl20_' path_clustering='pca1_cl20_clustering.xls' path_clustering_observables='pca1_cl20_clustering_observables.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls' xspace=2" $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
tail $ROUT/02_heatmaps.Rout



### Run tSNE
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsne_prefix='pca1_cl20_'  path_clustering_observables='pca1_cl20_clustering_observables.xls' tsne_pmin=1500" $RCODE/03_runtsne.R $ROUT/03_runtsne.Rout
tail $ROUT/03_runtsne.Rout


### Plot tSNE
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='pca1_cl20_'  path_rtsne_out='pca1_cl20_rtsne_out.rda' path_rtsne_data='pca1_cl20_rtsne_data.xls' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls'  tsne_cmin=1000" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
tail $ROUT/03_plottsne.Rout


### Cluster frequencies
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='pca1_cl20_' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
tail $ROUT/04_frequencies.Rout


##############################################################################
# Analysis of CK_2016-06-23_03 data
##############################################################################

RWD=$RWD_MAIN/CK_2016-06-23_03
ROUT=$RWD/Rout
mkdir -p $ROUT

### PCA scores
R CMD BATCH --no-save --no-restore "--args rwd='$RWD'" $RCODE/01_pcascores.R $ROUT/01_pcascores.Rout
tail $ROUT/01_pcascores.Rout


####################################
### PCA1 approach 

### FlowSOM clustering
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' flowsom_prefix='pca1_cl20_' pca_score_cutoff=0.9 pca_skip_top=0" $RCODE/02_flowsom.R $ROUT/02_flowsom.Rout
tail $ROUT/02_flowsom.Rout


### Heatmaps
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='pca1_cl20_' path_clustering='pca1_cl20_clustering.xls' path_clustering_observables='pca1_cl20_clustering_observables.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls' xspace=2" $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
tail $ROUT/02_heatmaps.Rout


### Run tSNE
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsne_prefix='pca1_cl20_'  path_clustering_observables='pca1_cl20_clustering_observables.xls' tsne_pmin=1500" $RCODE/03_runtsne.R $ROUT/03_runtsne.Rout
tail $ROUT/03_runtsne.Rout


### Plot tSNE
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='pca1_cl20_'  path_rtsne_out='pca1_cl20_rtsne_out.rda' path_rtsne_data='pca1_cl20_rtsne_data.xls' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls'  tsne_cmin=1000" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
tail $ROUT/03_plottsne.Rout


### Cluster frequencies
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='pca1_cl20_' path_clustering='pca1_cl20_clustering.xls' path_clustering_labels='pca1_cl20_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
tail $ROUT/04_frequencies.Rout


####################################
### PCA2 approach 

### FlowSOM clustering
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' flowsom_prefix='pca2_cl20_' pca_score_cutoff=0.9 pca_skip_top=2" $RCODE/02_flowsom.R $ROUT/02_flowsom.Rout
tail $ROUT/02_flowsom.Rout


### Heatmaps
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='pca2_cl20_' path_clustering='pca2_cl20_clustering.xls' path_clustering_observables='pca2_cl20_clustering_observables.xls' path_clustering_labels='pca2_cl20_clustering_labels.xls' xspace=2" $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
tail $ROUT/02_heatmaps.Rout


### Run tSNE
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsne_prefix='pca2_cl20_'  path_clustering_observables='pca2_cl20_clustering_observables.xls' tsne_pmin=1500" $RCODE/03_runtsne.R $ROUT/03_runtsne.Rout
tail $ROUT/03_runtsne.Rout


### Plot tSNE
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='pca2_cl20_'  path_rtsne_out='pca2_cl20_rtsne_out.rda' path_rtsne_data='pca2_cl20_rtsne_data.xls' path_clustering='pca2_cl20_clustering.xls' path_clustering_labels='pca2_cl20_clustering_labels.xls'  tsne_cmin=1000" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
tail $ROUT/03_plottsne.Rout




##########################################
# Analysis of CK_2016-06-23_03 merging2
##########################################

### Cluster merging
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' merging_prefix='pca1_merging2_' path_cluster_merging='cluster_merging2.xlsx' path_clustering='pca1_cl20_clustering.xls'" $RCODE/02_cluster_merging.R $ROUT/02_cluster_merging.Rout
tail $ROUT/02_cluster_merging.Rout


### Heatmaps
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='pca1_merging2_' path_clustering='pca1_merging2_clustering.xls' path_clustering_observables='pca1_cl20_clustering_observables.xls' path_clustering_labels='pca1_merging2_clustering_labels.xls' xspace=7" $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
tail $ROUT/02_heatmaps.Rout


### Plot tSNE
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='pca1_merging2_'  path_rtsne_out='pca1_cl20_rtsne_out.rda' path_rtsne_data='pca1_cl20_rtsne_data.xls' path_clustering='pca1_merging2_clustering.xls' path_clustering_labels='pca1_merging2_clustering_labels.xls'  tsne_cmin=1000" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
tail $ROUT/03_plottsne.Rout


### Cluster frequencies
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='pca1_merging2_' path_clustering='pca1_merging2_clustering.xls' path_clustering_labels='pca1_merging2_clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
tail $ROUT/04_frequencies.Rout




##############################################################################
# Analysis of CK_2016-06-29_01 data
##############################################################################





##############################################################################
# Analysis of CK_2016-06-29_02 data
##############################################################################



















