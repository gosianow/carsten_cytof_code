#!/bin/bash

###############################################################################################################
## Define paths to software and reference files

RWD_MAIN=/Users/gosia/Dropbox/UZH/carsten_cytof
RCODE=$RWD_MAIN
METADATA=$RWD_MAIN/MyeEUNITER
PANELS=$RWD_MAIN/MyeEUNITER


## Define which analysis to re-run
pcascores=true
select_observables=true
flowsom=true
heatmaps=true
runtsne=true
plottsne=true
frequencies=true
cluster_merging=true
cluster_extracting=true
cytokines=false
fcs_saving=false
pd1=false

###############################################################################################################
# Analysis of FACS_data data
# Use Analysis block 1
###############################################################################################################

RWD=$RWD_MAIN/MyeEUNITER
ROUT=$RWD/Rout
mkdir -p $ROUT
echo "$RWD"

file_panel="panel_mye.xlsx"
file_metadata="metadata_mye.xlsx"

pca_score_cutoff=0 ## Update!!!
rand_seed_consensus=1234
nmetaclusts=20

prefix_data="mye_"
prefix_panel=""
prefix_pca="pca1_"
prefix_clust="cl20_"

### PCA scores
if ${pcascores}; then
  echo "01_pcascores"
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' pcas_prefix='${prefix_data}${prefix_panel}' path_panel='${PANELS}/${file_panel}'" $RCODE/01_pcascores.R $ROUT/01_pcascores.Rout
tail $ROUT/01_pcascores.Rout
fi

### Select observables for clustering
if ${select_observables}; then
  echo "02_select_observables"
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' observ_prefix='${prefix_data}${prefix_panel}${prefix_pca}' path_pca_score='${prefix_data}${prefix_panel}princompscore_by_sample.xls' pca_score_cutoff=${pca_score_cutoff} pca_skip_top=0" $RCODE/02_select_observables.R $ROUT/02_select_observables.Rout
tail $ROUT/02_select_observables.Rout
fi

### FlowSOM clustering
if ${flowsom}; then
  echo "02_flowsom"
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' flowsom_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}' path_clustering_observables='${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' nmetaclusts=${nmetaclusts} rand_seed_consensus=${rand_seed_consensus}" $RCODE/02_flowsom.R $ROUT/02_flowsom.Rout
tail $ROUT/02_flowsom.Rout
fi

### Heatmaps
if ${heatmaps}; then
  echo "02_heatmaps"
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' heatmap_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}'  path_clustering_observables='${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' path_clustering='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering.xls'  path_clustering_labels='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering_labels.xls'" $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
tail $ROUT/02_heatmaps.Rout
fi

### Run tSNE
if ${runtsne}; then
  echo "03_runtsne"
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' tsne_prefix='${prefix_data}${prefix_panel}${prefix_pca}' path_clustering_observables='${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' tsne_pmin=1500" $RCODE/03_runtsne.R $ROUT/03_runtsne.Rout
tail $ROUT/03_runtsne.Rout
fi

### Plot tSNE - norm data
if ${plottsne}; then
  echo "03_plottsne"
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' tsnep_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}' path_rtsne_out='${prefix_data}${prefix_panel}${prefix_pca}rtsne_out_norm.rda' path_rtsne_data='${prefix_data}${prefix_panel}${prefix_pca}rtsne_data_norm.xls' path_clustering='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering.xls' path_clustering_labels='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10 tsnep_suffix='_norm'" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
tail $ROUT/03_plottsne.Rout
fi

### Plot tSNE - raw data
if ${plottsne}; then
  echo "03_plottsne"
R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' tsnep_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}' path_rtsne_out='${prefix_data}${prefix_panel}${prefix_pca}rtsne_out_raw.rda' path_rtsne_data='${prefix_data}${prefix_panel}${prefix_pca}rtsne_data_raw.xls' path_clustering='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering.xls' path_clustering_labels='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10 tsnep_suffix='_raw'" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
tail $ROUT/03_plottsne.Rout
fi











#
