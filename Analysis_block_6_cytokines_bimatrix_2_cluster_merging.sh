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
    --cytokines_bimatrix_cluster_merging)
    cytokines_bimatrix_cluster_merging=${2}
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
    --prefix_cytokines_cutoffs)
    prefix_cytokines_cutoffs=${2}
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
    --outdir)
    outdir=${2}
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

#############################################################################
### Cytokine analysis
### Cluster merging analysis
#############################################################################


if ${cytokines_bimatrix_cluster_merging}; then

  RWD=$RWD_MAIN/${data_dir}
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"

  prefix_cytokines="${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}${prefix_clsubset}${prefix_cytokines_cutoffs}"
  prefix_clust="cl${nmetaclusts}_"


  if [ ! -e "$RWD/${file_merging_cyt}" ]; then
    echo "File '$RWD/${file_merging_cyt}' does NOT exist! No merging analysis!"
    exit
  fi


  ### Cluster merging
  echo ">>> 02_cluster_merging"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' merging_prefix='${prefix_cytokines}${prefix_merging_cyt}' merging_outdir='${outdir}/01_clustering' path_cluster_merging='${file_merging_cyt}' path_clustering='${outdir}/01_clustering/${prefix_cytokines}${prefix_clust}clustering.xls'" $RCODE/02_cluster_merging.R $ROUT/02_cluster_merging.Rout
  tail $ROUT/02_cluster_merging.Rout


  if [ ! -e "$RWD/${outdir}/01_clustering/${prefix_cytokines}${prefix_merging_cyt}clustering.xls" ]; then
    echo "File '$RWD/${outdir}/01_clustering/${prefix_cytokines}${prefix_merging_cyt}clustering.xls' does NOT exist!"
    exit
  fi


  ### Heatmaps of merged clusters
  echo ">>> 02_heatmaps"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='${prefix_cytokines}${prefix_merging_cyt}' heatmap_outdir='${outdir}/01_clustering' path_data='${outdir}/01_clustering/${prefix_cytokines}bimatrix.txt' path_metadata='${METADATA}/${file_metadata}'   path_clustering_observables='${outdir}/01_clustering/${prefix_cytokines}clustering_observables.xls' path_clustering='${outdir}/01_clustering/${prefix_cytokines}${prefix_merging_cyt}clustering.xls'  path_clustering_labels='${outdir}/01_clustering/${prefix_cytokines}${prefix_merging_cyt}clustering_labels.xls' path_marker_selection='${prefix_cytokines}marker_selection.txt' aggregate_fun='mean' pheatmap_palette='RdYlBu' pheatmap_palette_rev=TRUE pheatmap_scale=FALSE" $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
  tail $ROUT/02_heatmaps.Rout


  ### Get cluster frequencies
  echo ">>> 04_frequencies"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${prefix_cytokines}${prefix_merging_cyt}' freq_outdir='${outdir}/03_frequencies' path_metadata='${METADATA}/${file_metadata}'  path_clustering='${outdir}/01_clustering/${prefix_cytokines}${prefix_merging_cyt}clustering.xls' path_clustering_labels='${outdir}/01_clustering/${prefix_cytokines}${prefix_merging_cyt}clustering_labels.xls' path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_1dataset_3responses.R'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
  tail $ROUT/04_frequencies.Rout


fi
































#
