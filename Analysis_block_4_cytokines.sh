#!/bin/bash

# -----------------------------------------------------
# argument parcing
# -----------------------------------------------------

while [[ ${1} ]]; do
  case "${1}" in
    --RCODE)
    RCODE=${2}
    shift
    ;;
    --RWD_MAIN)
    RWD_MAIN=${2}
    shift
    ;;
    --data_dir)
    data_dir=${2}
    shift
    ;;
    --cytokines)
    cytokines=${2}
    shift
    ;;
    --PANELS)
    PANELS=${2}
    shift
    ;;
    --METADATA)
    METADATA=${2}
    shift
    ;;
    --file_metadata)
    file_metadata=${2}
    shift
    ;;
    --prefix_data)
    prefix_data=${2}
    shift
    ;;
    --prefix_panel)
    prefix_panel=${2}
    shift
    ;;
    --prefix_pca)
    prefix_pca=${2}
    shift
    ;;
    --prefix_merging)
    prefix_merging=${2}
    shift
    ;;
    --prefix_clsubset)
    prefix_clsubset=${2}
    shift
    ;;
    --file_cytokines_cutoffs)
    file_cytokines_cutoffs=${2}
    shift
    ;;
    --prefix_cytokines_cutoffs)
    prefix_cytokines_cutoffs=${2}
    shift
    ;;
    --clsubset)
    clsubset=${2}
    shift
    ;;
    --nmetaclusts)
    nmetaclusts=${2}
    shift
    ;;
    --file_merging_cyt)
    file_merging_cyt=${2}
    shift
    ;;
    --prefix_merging_cyt)
    prefix_merging_cyt=${2}
    shift
    ;;

    *)
    echo "Unknown parameter: ${1}" >&2
  esac

  if ! shift; then
    echo 'Missing parameter argument.' >&2
  fi
done

# -----------------------------------------------------
# function
# -----------------------------------------------------

### Cytokine analysis
if ${cytokines}; then

  RWD=$RWD_MAIN/${data_dir}
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"

  ### Analysis based on different cutoffs for base and tx
  prefix_cytokines="${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}${prefix_clsubset}${prefix_cytokines_cutoffs}raw2_"
  prefix_clust="cl${nmetaclusts}_"

  ### Create the bimatrix
  echo "06_cytokines_bimatrix"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' cytokines_prefix='${prefix_cytokines}' cytokines_outdir='060_cytokines' path_data='010_data/${prefix_data}${prefix_panel}expr_raw.rds' path_metadata='${METADATA}/${file_metadata}' path_cytokines_cutoffs='${PANELS}/${file_cytokines_cutoffs}' path_clustering='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls' path_clustering_labels='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering_labels.xls' clsubset=${clsubset} cutoff_colname=c('positive_cutoff_raw_base','positive_cutoff_raw_tx')" $RCODE/06_cytokines_bimatrix.R $ROUT/06_cytokines_bimatrix.Rout
  tail $ROUT/06_cytokines_bimatrix.Rout


  ### FlowSOM clustering of bimatrix
  echo "02_flowsom"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' flowsom_prefix='${prefix_cytokines}${prefix_clust}' flowsom_outdir='060_cytokines' path_data='060_cytokines/${prefix_cytokines}bimatrix.txt' path_clustering_observables='060_cytokines/${prefix_cytokines}clustering_observables.xls' nmetaclusts=${nmetaclusts} rand_seed_consensus=1234" $RCODE/02_flowsom.R $ROUT/02_flowsom.Rout
  tail $ROUT/02_flowsom.Rout


  ### Heatmaps
  echo "02_heatmaps"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='${prefix_cytokines}${prefix_clust}' heatmap_outdir='060_cytokines' path_data='060_cytokines/${prefix_cytokines}bimatrix.txt' path_metadata='${METADATA}/${file_metadata}'   path_clustering_observables='060_cytokines/${prefix_cytokines}clustering_observables.xls' path_clustering='060_cytokines/${prefix_cytokines}${prefix_clust}clustering.xls'  path_clustering_labels='060_cytokines/${prefix_cytokines}${prefix_clust}clustering_labels.xls' path_marker_selection='${prefix_cytokines}marker_selection.txt' aggregate_fun='mean' pheatmap_palette='RdYlBu' pheatmap_palette_rev=TRUE pheatmap_scale=FALSE" $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
  tail $ROUT/02_heatmaps.Rout


  ### tSNE plot with bimatrix clusters (based on raw data)
  echo "03_plottsne"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='${prefix_cytokines}${prefix_clust}raw_' tsnep_outdir='060_cytokines' path_metadata='${METADATA}/${file_metadata}'  path_rtsne_out='040_tsnemaps/${prefix_data}${prefix_panel}${prefix_pca}raw_rtsne_out.rda' path_rtsne_data='040_tsnemaps/${prefix_data}${prefix_panel}${prefix_pca}raw_rtsne_data.xls' path_clustering='060_cytokines/${prefix_cytokines}${prefix_clust}clustering.xls' path_clustering_labels='060_cytokines/${prefix_cytokines}${prefix_clust}clustering_labels.xls' pdf_width=15 pdf_height=10" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
  tail $ROUT/03_plottsne.Rout


  ### Get cluster frequencies
  echo "04_frequencies"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${prefix_cytokines}${prefix_clust}' freq_outdir='060_cytokines' path_metadata='${METADATA}/${file_metadata}'  path_clustering='060_cytokines/${prefix_cytokines}${prefix_clust}clustering.xls' path_clustering_labels='060_cytokines/${prefix_cytokines}${prefix_clust}clustering_labels.xls' path_fun_models='$RCODE/00_models.R'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
  tail $ROUT/04_frequencies.Rout

  #############################################################################
  ### Cluster merging analysis
  #############################################################################

  if [ ! -e "$RWD/${file_merging_cyt}" ]; then
    echo "File '$RWD/${file_merging_cyt}' does NOT exist! No merging analysis!"
    exit
  fi


  ### Cluster merging
  echo "02_cluster_merging"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' merging_prefix='${prefix_cytokines}${prefix_merging_cyt}' merging_outdir='060_cytokines' path_cluster_merging='${file_merging_cyt}' path_clustering='060_cytokines/${prefix_cytokines}${prefix_clust}clustering.xls'" $RCODE/02_cluster_merging.R $ROUT/02_cluster_merging.Rout
  tail $ROUT/02_cluster_merging.Rout


  if [ ! -e "$RWD/060_cytokines/${prefix_cytokines}${prefix_merging_cyt}clustering.xls" ]; then
    echo "File '$RWD/060_cytokines/${prefix_cytokines}${prefix_merging_cyt}clustering.xls' does NOT exist!"
    exit
  fi


  ### Heatmaps of merged clusters
  echo "02_heatmaps"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='${prefix_cytokines}${prefix_merging_cyt}' heatmap_outdir='060_cytokines' path_data='060_cytokines/${prefix_cytokines}bimatrix.txt' path_metadata='${METADATA}/${file_metadata}'   path_clustering_observables='060_cytokines/${prefix_cytokines}clustering_observables.xls' path_clustering='060_cytokines/${prefix_cytokines}${prefix_merging_cyt}clustering.xls'  path_clustering_labels='060_cytokines/${prefix_cytokines}${prefix_merging_cyt}clustering_labels.xls' path_marker_selection='${prefix_cytokines}marker_selection.txt' aggregate_fun='mean' pheatmap_palette='RdYlBu' pheatmap_palette_rev=TRUE pheatmap_scale=FALSE" $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
  tail $ROUT/02_heatmaps.Rout


  ### tSNE plot with bimatrix clusters (based on raw data)
  echo "03_plottsne"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='${prefix_cytokines}${prefix_merging_cyt}raw_' tsnep_outdir='060_cytokines' path_metadata='${METADATA}/${file_metadata}'  path_rtsne_out='040_tsnemaps/${prefix_data}${prefix_panel}${prefix_pca}raw_rtsne_out.rda' path_rtsne_data='040_tsnemaps/${prefix_data}${prefix_panel}${prefix_pca}raw_rtsne_data.xls' path_clustering='060_cytokines/${prefix_cytokines}${prefix_merging_cyt}clustering.xls' path_clustering_labels='060_cytokines/${prefix_cytokines}${prefix_merging_cyt}clustering_labels.xls' pdf_width=15 pdf_height=10" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
  tail $ROUT/03_plottsne.Rout


  ### Get cluster frequencies
  echo "04_frequencies"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${prefix_cytokines}${prefix_merging_cyt}' freq_outdir='060_cytokines' path_metadata='${METADATA}/${file_metadata}'  path_clustering='060_cytokines/${prefix_cytokines}${prefix_merging_cyt}clustering.xls' path_clustering_labels='060_cytokines/${prefix_cytokines}${prefix_merging_cyt}clustering_labels.xls' path_fun_models='$RCODE/00_models.R'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
  tail $ROUT/04_frequencies.Rout


fi
































#
