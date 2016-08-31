#!/bin/bash

###############################################################################################################
## Define paths to software and reference files

RWD_MAIN=/Users/gosia/Dropbox/UZH/carsten_cytof
RCODE=/Users/gosia/Dropbox/UZH/carsten_cytof_code
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


###############################################################################################################
# Analysis of FACS_data data
# Use Analysis block 1
###############################################################################################################

RWD=$RWD_MAIN/MyeEUNITER
ROUT=$RWD/Rout
mkdir -p $ROUT
echo "$RWD"

file_panel="panel_MyeEUNITER.xlsx"
file_metadata="metadata_MyeEUNITER.xlsx"

pca_score_cutoff=0 # We keep all the markers!
rand_seed_consensus=1234

tsne_pmin=5000 # In the CK analysis, I use 1500 per sample.

prefix_data="mye_"
prefix_panel=""
prefix_pca="pca1_"



# ### Data normalization
# if ${data_normalization}; then
#   echo "01_data_normalization"
#   R CMD BATCH --no-save --no-restore "--args rwd='$RWD' data_prefix='${prefix_data}${prefix_panel}' data_outdir='010_data' path_metadata='${METADATA}/${file_metadata}' path_panel='${PANELS}/${file_panel}'" $RCODE/01_data_normalization.R $ROUT/01_data_normalization.Rout
#   tail $ROUT/01_data_normalization.Rout
# fi
#
# ### PCA scores based on raw data
# if ${pcascores}; then
#   echo "01_pcascores"
#   R CMD BATCH --no-save --no-restore "--args rwd='$RWD' pcas_prefix='${prefix_data}${prefix_panel}' pcas_outdir='020_pcascores' path_data='010_data/${prefix_data}${prefix_panel}expr_raw.rds' path_metadata='${METADATA}/${file_metadata}' path_panel='${PANELS}/${file_panel}'" $RCODE/01_pcascores.R $ROUT/01_pcascores.Rout
#   tail $ROUT/01_pcascores.Rout
# fi
#
# ### Select observables for clustering
# if ${select_observables}; then
#   echo "02_select_observables"
#   R CMD BATCH --no-save --no-restore "--args rwd='$RWD' observ_prefix='${prefix_data}${prefix_panel}${prefix_pca}' observ_outdir='030_heatmaps' path_pca_score='020_pcascores/${prefix_data}${prefix_panel}princompscore_by_sample.xls' pca_score_cutoff=${pca_score_cutoff} pca_skip_top=0" $RCODE/02_select_observables.R $ROUT/02_select_observables.Rout
#   tail $ROUT/02_select_observables.Rout
# fi
#
#
# ### Run tSNE
# if ${runtsne}; then
#   echo "03_runtsne"
#
#   # on raw data
#   R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsne_prefix='${prefix_data}${prefix_panel}${prefix_pca}raw_' tsne_outdir='040_tsnemaps' path_data='010_data/${prefix_data}${prefix_panel}expr_raw.rds' path_metadata='${METADATA}/${file_metadata}'  path_clustering_observables='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' tsne_pmin=${tsne_pmin}" $RCODE/03_runtsne.R $ROUT/03_runtsne.Rout
#   tail $ROUT/03_runtsne.Rout
#
# fi
#
#
# ### Plot tSNE with marker expression as a heat
# if ${plottsne}; then
#   echo "03_plottsne"
#
#   # raw tSNE and raw marker expression
#   R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='${prefix_data}${prefix_panel}${prefix_pca}raw_mark_raw_' tsnep_outdir='040_tsnemaps_expr' path_data='010_data/${prefix_data}${prefix_panel}expr_raw.rds' path_metadata='${METADATA}/${file_metadata}'  path_rtsne_out='040_tsnemaps/${prefix_data}${prefix_panel}${prefix_pca}raw_rtsne_out.rda' path_rtsne_data='040_tsnemaps/${prefix_data}${prefix_panel}${prefix_pca}raw_rtsne_data.xls' path_clustering_observables='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' pdf_width=15 pdf_height=10" $RCODE/03_plottsne_expr.R $ROUT/03_plottsne_expr.Rout
#   tail $ROUT/03_plottsne_expr.Rout
#
#   # raw tSNE and 01 normalized marker expression
#   R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='${prefix_data}${prefix_panel}${prefix_pca}raw_mark_norm_' tsnep_outdir='040_tsnemaps_expr' path_data='010_data/${prefix_data}${prefix_panel}expr_norm.rds' path_metadata='${METADATA}/${file_metadata}'  path_rtsne_out='040_tsnemaps/${prefix_data}${prefix_panel}${prefix_pca}raw_rtsne_out.rda' path_rtsne_data='040_tsnemaps/${prefix_data}${prefix_panel}${prefix_pca}raw_rtsne_data.xls' path_clustering_observables='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' pdf_width=15 pdf_height=10" $RCODE/03_plottsne_expr.R $ROUT/03_plottsne_expr.Rout
#   tail $ROUT/03_plottsne_expr.Rout
#
# fi




nmetaclusts=(20 4)
prefix_clust=("cl20_" "cl4_")

for i in "${nmetaclusts[@]}"
do

  ### FlowSOM clustering based on raw data
  if ${flowsom}; then
    echo "02_flowsom"
    R CMD BATCH --no-save --no-restore "--args rwd='$RWD' flowsom_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust[$i]}' flowsom_outdir='030_heatmaps' path_data='010_data/${prefix_data}${prefix_panel}expr_raw.rds' path_clustering_observables='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' nmetaclusts=${nmetaclusts[$i]} rand_seed_consensus=${rand_seed_consensus}" $RCODE/02_flowsom.R $ROUT/02_flowsom.Rout
    tail $ROUT/02_flowsom.Rout
  fi

  ### Heatmaps
  if ${heatmaps}; then
    echo "02_heatmaps"

    # based on raw expression
    R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust[$i]}raw_' heatmap_outdir='030_heatmaps' path_data='010_data/${prefix_data}${prefix_panel}expr_raw.rds' path_metadata='${METADATA}/${file_metadata}'   path_clustering_observables='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' path_clustering='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust[$i]}clustering.xls'  path_clustering_labels='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust[$i]}clustering_labels.xls' path_marker_selection='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust[$i]}marker_selection.txt' aggregate_fun='median' pheatmap_palette='YlGnBu' pheatmap_palette_rev=FALSE pheatmap_scale=TRUE" $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
    tail $ROUT/02_heatmaps.Rout

    # based on 01 normalized expression
    R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust[$i]}norm_' heatmap_outdir='030_heatmaps' path_data='010_data/${prefix_data}${prefix_panel}expr_norm.rds' path_metadata='${METADATA}/${file_metadata}'   path_clustering_observables='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' path_clustering='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust[$i]}clustering.xls'  path_clustering_labels='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust[$i]}clustering_labels.xls' path_marker_selection='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust[$i]}marker_selection.txt' aggregate_fun='median' pheatmap_palette='RdYlBu' pheatmap_palette_rev=TRUE pheatmap_scale=TRUE" $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
    tail $ROUT/02_heatmaps.Rout

  fi


  ### Plot tSNE
  if ${plottsne}; then
    echo "03_plottsne"

    ### Based on raw data
    R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust[$i]}raw_' tsnep_outdir='040_tsnemaps' path_metadata='${METADATA}/${file_metadata}'  path_rtsne_out='040_tsnemaps/${prefix_data}${prefix_panel}${prefix_pca}raw_rtsne_out.rda' path_rtsne_data='040_tsnemaps/${prefix_data}${prefix_panel}${prefix_pca}raw_rtsne_data.xls' path_clustering='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust[$i]}clustering.xls' path_clustering_labels='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust[$i]}clustering_labels.xls'  pdf_width=15 pdf_height=10" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
    tail $ROUT/03_plottsne.Rout

  fi

done

























#
