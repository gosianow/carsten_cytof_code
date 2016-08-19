#!/bin/bash
## Define functions

Analysis_block_1_main () {

  ### PCA scores
  if ${pcascores}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' pcas_prefix='${prefix_data}${prefix_panel}' path_panel='${PANELS}/${file_panel}'" $RCODE/01_pcascores.R $ROUT/01_pcascores.Rout
  tail $ROUT/01_pcascores.Rout
  fi

  ### Select observables for clustering
  if ${select_observables}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' observ_prefix='${prefix_data}${prefix_panel}${prefix_pca}' path_pca_score='${prefix_data}${prefix_panel}princompscore_by_sample.xls' pca_score_cutoff=${pca_score_cutoff} pca_skip_top=0" $RCODE/02_select_observables.R $ROUT/02_select_observables.Rout
  tail $ROUT/02_select_observables.Rout
  fi

  ### FlowSOM clustering
  if ${flowsom}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' flowsom_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}' path_clustering_observables='${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' nmetaclusts=${nmetaclusts} rand_seed_consensus=${rand_seed_consensus}" $RCODE/02_flowsom.R $ROUT/02_flowsom.Rout
  tail $ROUT/02_flowsom.Rout
  fi

  ### Heatmaps
  if ${heatmaps}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' heatmap_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}' path_panel='${PANELS}/${file_panel}' path_clustering_observables='${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' path_clustering='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering.xls'  path_clustering_labels='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering_labels.xls' path_pca_score='${prefix_data}${prefix_panel}princompscore_by_sample.xls' " $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
  tail $ROUT/02_heatmaps.Rout
  fi

  ### Run tSNE
  if ${runtsne}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' tsne_prefix='${prefix_data}${prefix_panel}${prefix_pca}' path_clustering_observables='${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' tsne_pmin=1500" $RCODE/03_runtsne.R $ROUT/03_runtsne.Rout
  tail $ROUT/03_runtsne.Rout
  fi

  ### Plot tSNE based on normalized data
  if ${plottsne}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' tsnep_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}' path_rtsne_out='${prefix_data}${prefix_panel}${prefix_pca}rtsne_out_norm.rda' path_rtsne_data='${prefix_data}${prefix_panel}${prefix_pca}rtsne_data_norm.xls' path_clustering='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering.xls' path_clustering_labels='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10 tsnep_suffix='_norm'" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
  tail $ROUT/03_plottsne.Rout
  fi

  ### Plot tSNE based on raw data
  if ${plottsne}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' tsnep_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}' path_rtsne_out='${prefix_data}${prefix_panel}${prefix_pca}rtsne_out_raw.rda' path_rtsne_data='${prefix_data}${prefix_panel}${prefix_pca}rtsne_data_raw.xls' path_clustering='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering.xls' path_clustering_labels='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10 tsnep_suffix='_raw'" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
  tail $ROUT/03_plottsne.Rout
  fi

  ### Get cluster frequencies
  if ${frequencies}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' freq_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}' path_clustering='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering.xls' path_clustering_labels='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
  tail $ROUT/04_frequencies.Rout
  fi


}


Analysis_block_2_cluster_merging () {

  ### Cluster merging
  if ${cluster_merging}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' merging_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}' path_cluster_merging='${file_merging}' path_clustering='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering.xls'" $RCODE/02_cluster_merging.R $ROUT/02_cluster_merging.Rout
  tail $ROUT/02_cluster_merging.Rout
  fi

  if [ ! -e "$RWD/030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls" ]; then
    echo "File '$RWD/030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls' does NOT exist!"
    exit
  fi

  ### Heatmaps
  if ${heatmaps}; then
      R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' heatmap_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}' path_panel='${PANELS}/${file_panel}' path_clustering_observables='${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' path_clustering='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls'  path_clustering_labels='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering_labels.xls' path_pca_score='${prefix_data}${prefix_panel}princompscore_by_sample.xls' " $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
      tail $ROUT/02_heatmaps.Rout
  fi

  ### Plot tSNE
  if ${plottsne}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' tsnep_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}' path_rtsne_out='${prefix_data}${prefix_panel}${prefix_pca}rtsne_out_norm.rda' path_rtsne_data='${prefix_data}${prefix_panel}${prefix_pca}rtsne_data_norm.xls' path_clustering='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls' path_clustering_labels='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10 tsnep_suffix='_norm'" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
  tail $ROUT/03_plottsne.Rout
  fi

  ### Get cluster frequencies
  if ${frequencies}; then
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' freq_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}' path_clustering='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls' path_clustering_labels='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
  tail $ROUT/04_frequencies.Rout
  fi


}


Analysis_block_3_cluster_extracting () {

  ### Cluster extracting
  if ${cluster_extracting}; then
  for i in 0 1
  do

    R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' path_clustering='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls' path_clustering_labels='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering_labels.xls' extract_cluster='${extract_cluster[$i]}' extract_dir='${RWD_MAIN}/${extract_dir[$i]}/010_cleanfcs'" $RCODE/02_cluster_extracting.R $ROUT/02_cluster_extracting.Rout
    tail $ROUT/02_cluster_extracting.Rout

  done
  fi

}


Analysis_block_4_main_CD4_CD8 () {

  ### PCA scores
  if ${pcascores}; then
  for i in 0 1
  do

    RWD=$RWD_MAIN/${data_dir[$i]}
    ROUT=$RWD/Rout
    mkdir -p $ROUT
    echo "$RWD"

    ### PCA scores
    R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' pcas_prefix='${prefix_data[$i]}${prefix_panel[$i]}' path_panel='${PANELS}/${file_panel[$i]}'" $RCODE/01_pcascores.R $ROUT/01_pcascores.Rout
    tail $ROUT/01_pcascores.Rout

  done
  fi


  for i in 0 1
  do

    RWD=$RWD_MAIN/${data_dir[$i]}
    ROUT=$RWD/Rout
    mkdir -p $ROUT
    echo "$RWD"

    ### Select observables for clustering
    if ${select_observables}; then
    R CMD BATCH --no-save --no-restore "--args rwd='$RWD' observ_prefix='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}' path_pca_score='${prefix_data[$i]}${prefix_panel[$i]}princompscore_by_sample.xls' pca_score_cutoff=${pca_score_cutoff[$i]} pca_skip_top=0" $RCODE/02_select_observables.R $ROUT/02_select_observables.Rout
    tail $ROUT/02_select_observables.Rout
    fi

    ### FlowSOM clustering
    if ${flowsom}; then
    R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' flowsom_prefix='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}${prefix_clust}' path_clustering_observables='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}clustering_observables.xls' nmetaclusts=${nmetaclusts} rand_seed_consensus=${rand_seed_consensus}" $RCODE/02_flowsom.R $ROUT/02_flowsom.Rout
    tail $ROUT/02_flowsom.Rout
    fi

    ### Heatmaps
    if ${heatmaps}; then
    R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' heatmap_prefix='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}${prefix_clust}' path_panel='${PANELS}/${file_panel[$i]}' path_clustering_observables='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}clustering_observables.xls' path_clustering='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}${prefix_clust}clustering.xls'  path_clustering_labels='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}${prefix_clust}clustering_labels.xls' path_pca_score='${prefix_data[$i]}${prefix_panel[$i]}princompscore_by_sample.xls' " $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
    tail $ROUT/02_heatmaps.Rout
    fi

    ### Run tSNE
    if ${runtsne}; then
    R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' tsne_prefix='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}' path_clustering_observables='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}clustering_observables.xls' tsne_pmin=1500" $RCODE/03_runtsne.R $ROUT/03_runtsne.Rout
    tail $ROUT/03_runtsne.Rout
    fi

    ### Plot tSNE
    if ${plottsne}; then
    R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' tsnep_prefix='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}${prefix_clust}' path_rtsne_out='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}rtsne_out_norm.rda' path_rtsne_data='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}rtsne_data_norm.xls' path_clustering='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}${prefix_clust}clustering.xls' path_clustering_labels='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}${prefix_clust}clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10 tsnep_suffix='_norm'" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
    tail $ROUT/03_plottsne.Rout
    fi

    ### Get cluster frequencies
    if ${frequencies}; then
    R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' freq_prefix='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}${prefix_clust}' path_clustering='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}${prefix_clust}clustering.xls' path_clustering_labels='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}${prefix_clust}clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
    tail $ROUT/04_frequencies.Rout
    fi


  done


}

Analysis_block_5_cluster_merging_CD4_CD8 () {

  for i in 0 1
  do

    RWD=$RWD_MAIN/${data_dir[$i]}
    ROUT=$RWD/Rout
    mkdir -p $ROUT
    echo "$RWD"

    ### Cluster merging
    if ${cluster_merging}; then
    R CMD BATCH --no-save --no-restore "--args rwd='$RWD' merging_prefix='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}${prefix_merging[$i]}' path_cluster_merging='${file_merging[$i]}' path_clustering='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}${prefix_clust}clustering.xls'" $RCODE/02_cluster_merging.R $ROUT/02_cluster_merging.Rout
    tail $ROUT/02_cluster_merging.Rout
    fi

    if [ ! -e "$RWD/030_heatmaps/${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}${prefix_merging[$i]}clustering.xls" ]; then
      echo "File '$RWD/030_heatmaps/${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}${prefix_merging[$i]}clustering.xls' does NOT exist!"
      exit
    fi

    ### Heatmaps
    if ${heatmaps}; then
    R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' heatmap_prefix='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}${prefix_merging[$i]}' path_panel='${PANELS}/${file_panel[$i]}' path_clustering_observables='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}clustering_observables.xls' path_clustering='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}${prefix_merging[$i]}clustering.xls'  path_clustering_labels='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}${prefix_merging[$i]}clustering_labels.xls' path_pca_score='${prefix_data[$i]}${prefix_panel[$i]}princompscore_by_sample.xls' " $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
    tail $ROUT/02_heatmaps.Rout
    fi

    ### Plot tSNE
    if ${plottsne}; then
    R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' tsnep_prefix='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}${prefix_merging[$i]}' path_rtsne_out='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}rtsne_out_norm.rda' path_rtsne_data='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}rtsne_data_norm.xls' path_clustering='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}${prefix_merging[$i]}clustering.xls' path_clustering_labels='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}${prefix_merging[$i]}clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10 tsnep_suffix='_norm'" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
    tail $ROUT/03_plottsne.Rout
    fi

    ### Get cluster frequencies
    if ${frequencies}; then
    R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' freq_prefix='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}${prefix_merging[$i]}' path_clustering='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}${prefix_merging[$i]}clustering.xls' path_clustering_labels='${prefix_data[$i]}${prefix_panel[$i]}${prefix_pca}${prefix_merging[$i]}clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
    tail $ROUT/04_frequencies.Rout
    fi

  done


}


## Define paths to software and reference files

RWD_MAIN=/Users/gosia/Dropbox/UZH/carsten_cytof
RCODE=$RWD_MAIN
METADATA=$RWD_MAIN/CK_metadata
PANELS=$RWD_MAIN/CK_panels


## Define which analysis to re-run
pcascores=false
select_observables=false
flowsom=false
heatmaps=false
runtsne=true
plottsne=true
frequencies=false
cluster_merging=false
cluster_extracting=false
cytokines=false
fcs_saving=false
pd1=false



###############################################################################################################

###############################################################################################################
# Analysis of CK_2016-06-29_03 data
# Use Analysis block 1
###############################################################################################################

DATA=29
PANEL=3
RWD=$RWD_MAIN/CK_2016-06-29_03
ROUT=$RWD/Rout
mkdir -p $ROUT
echo "$RWD"

file_panel="panel3.xlsx"
file_metadata="metadata_29_03.xlsx"

pca_score_cutoff=0.9
rand_seed_consensus=1234
nmetaclusts=20

prefix_data="29_"
prefix_panel="03_"
prefix_pca="pca1_"
prefix_clust="cl20_"

Analysis_block_1_main

















#
