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
    --pd1_bimatrix_main)
    pd1_bimatrix_main=${2}
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
    --som_dim)
    som_dim=${2}
    shift
    ;;
    --nmetaclusts)
    nmetaclusts=${2}
    shift
    ;;


    *)
    echo "Unknown parameter: ${1}" >&2
    return 1
  esac

  if ! shift; then
    echo 'Missing parameter argument.' >&2
    return 1
  fi
done

# -----------------------------------------------------
# function
# -----------------------------------------------------

if ${pd1_bimatrix_main}; then

  RWD=$RWD_MAIN/${data_dir}
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"

  ### Use different cutoffs for base and tx and raw data
  prefix_pd1main="${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}${prefix_clsubset}${prefix_cytokines_cutoffs}"
  prefix_clust="cl${nmetaclusts}_"

  ### Create the bimatrix
  echo ">>> 07_pd1_bimatrix"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' pd1_prefix='${prefix_pd1main}' pd1_outdir='070_pd1_bimatrix/01_clustering' path_data='010_data/${prefix_data}${prefix_panel}expr_raw.rds' path_metadata='${METADATA}/${file_metadata}' path_cytokines_cutoffs='${PANELS}/${file_cytokines_cutoffs}' path_clustering='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls' path_clustering_labels='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering_labels.xls' clsubset=${clsubset} cutoff_colname=c('positive_cutoff_raw_base','positive_cutoff_raw_tx')" $RCODE/07_pd1_bimatrix.R $ROUT/07_pd1_bimatrix.Rout
  tail $ROUT/07_pd1_bimatrix.Rout


  ### Frequency analysis of PD1+ and PD1-
  echo ">>> 04_frequencies"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${prefix_pd1main}' freq_outdir='070_pd1_bimatrix/03_frequencies_auto' path_metadata='${METADATA}/${file_metadata}'  path_clustering='070_pd1_bimatrix/01_clustering/${prefix_pd1main}clustering.xls' path_clustering_labels='070_pd1_bimatrix/01_clustering/${prefix_pd1main}clustering_labels.xls' path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_1dataset_3responses.R'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
  tail $ROUT/04_frequencies.Rout


  ### Analysis of cytokines for PD1+
  prefix_pd1type="positive_"


  ### Clustering of bimatrix is based on SOM only when som_dim^2 = nmetaclusts, otherwise cluster consesnsus is applied additionally
  echo ">>> 06_cytokines_bimatrix_clustering"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' clust_prefix='${prefix_pd1main}${prefix_pd1type}${prefix_clust}' clust_outdir='070_pd1_bimatrix/01_clustering' path_data='070_pd1_bimatrix/01_clustering/${prefix_pd1main}${prefix_pd1type}bimatrix.txt' path_clustering_observables='070_pd1_bimatrix/01_clustering/${prefix_pd1main}${prefix_pd1type}clustering_observables.xls' som_dim=${som_dim} nmetaclusts=${nmetaclusts}" $RCODE/06_cytokines_bimatrix_clustering.R $ROUT/06_cytokines_bimatrix_clustering.Rout
  tail $ROUT/06_cytokines_bimatrix_clustering.Rout


  ### Heatmaps
  echo ">>> 02_heatmaps"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='${prefix_pd1main}${prefix_pd1type}${prefix_clust}' heatmap_outdir='070_pd1_bimatrix/01_clustering' path_data='070_pd1_bimatrix/01_clustering/${prefix_pd1main}${prefix_pd1type}bimatrix.txt' path_metadata='${METADATA}/${file_metadata}'   path_clustering_observables='070_pd1_bimatrix/01_clustering/${prefix_pd1main}clustering_observables.xls' path_clustering='070_pd1_bimatrix/01_clustering/${prefix_pd1main}${prefix_pd1type}${prefix_clust}clustering.xls'  path_clustering_labels='070_pd1_bimatrix/01_clustering/${prefix_pd1main}${prefix_pd1type}${prefix_clust}clustering_labels.xls' path_marker_selection='${prefix_pd1main}marker_selection.txt' aggregate_fun='mean' pheatmap_palette='RdYlBu' pheatmap_palette_rev=TRUE pheatmap_scale=FALSE" $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
  tail $ROUT/02_heatmaps.Rout


  ### Frequency analysis of cytokines
  echo ">>> 04_frequencies"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${prefix_pd1main}${prefix_pd1type}${prefix_clust}' freq_outdir='070_pd1_bimatrix/03_frequencies_auto' path_metadata='${METADATA}/${file_metadata}'  path_clustering='070_pd1_bimatrix/01_clustering/${prefix_pd1main}${prefix_pd1type}${prefix_clust}clustering.xls' path_clustering_labels='070_pd1_bimatrix/01_clustering/${prefix_pd1main}${prefix_pd1type}${prefix_clust}clustering_labels.xls' path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_1dataset_3responses.R' pdf_hight=8" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
  tail $ROUT/04_frequencies.Rout




fi











#
