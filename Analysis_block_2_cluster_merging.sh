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
    --cluster_merging)
    cluster_merging=${2}
    shift
    ;;
    --heatmaps)
    heatmaps=${2}
    shift
    ;;
    --plottsne)
    plottsne=${2}
    shift
    ;;
    --frequencies)
    frequencies=${2}
    shift
    ;;
    --expression)
    expression=${2}
    shift
    ;;
    --METADATA)
    METADATA=${2}
    shift
    ;;
    --PANELS)
    PANELS=${2}
    shift
    ;;
    --file_metadata)
    file_metadata=${2}
    shift
    ;;
    --file_panel)
    file_panel=${2}
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
    --prefix_clust)
    prefix_clust=${2}
    shift
    ;;
    --prefix_merging)
    prefix_merging=${2}
    shift
    ;;
    --file_merging)
    file_merging=${2}
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

RWD=$RWD_MAIN/${data_dir}
ROUT=$RWD/Rout
mkdir -p $ROUT

if ${cluster_merging} || ${heatmaps} || ${plottsne} || ${frequencies} || ${expression}
then
  echo "$RWD"
fi

### Cluster merging
if ${cluster_merging}; then
  echo ">>> 02_cluster_merging"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' merging_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}' merging_outdir='030_heatmaps' path_cluster_merging='${file_merging}' path_clustering='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering.xls'" $RCODE/02_cluster_merging.R $ROUT/02_cluster_merging.Rout
  tail $ROUT/02_cluster_merging.Rout

  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' codes_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}' codes_outdir='030_heatmaps' path_codes='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}codes.xls' pdf_width=22 pdf_height=7 path_cluster_merging='${file_merging}'" $RCODE/02_som_codes.R $ROUT/02_som_codes.Rout
  tail $ROUT/02_som_codes.Rout
fi

if [ ! -e "$RWD/030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls" ]; then
  echo "File '$RWD/030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls' does NOT exist!"
  exit
fi

### Heatmaps
if ${heatmaps}; then
  echo ">>> 02_heatmaps"

  # based on raw data
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}raw_' heatmap_outdir='030_heatmaps' path_data='010_data/${prefix_data}${prefix_panel}expr_raw.rds' path_metadata='${METADATA}/${file_metadata}'   path_clustering_observables='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' path_clustering='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls'  path_clustering_labels='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering_labels.xls' path_marker_selection='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}marker_selection.txt' aggregate_fun='median' pheatmap_palette='YlGnBu' pheatmap_palette_rev=FALSE pheatmap_scale=TRUE extra_path_data='010_data/${prefix_data}${prefix_panel}expr_norm.rds' extra_pheatmap_palette='RdYlBu' extra_pheatmap_palette_rev=TRUE extra_suffix='_norm'" $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
  tail $ROUT/02_heatmaps.Rout

fi

### Plot tSNE
if ${plottsne}; then
  echo ">>> 03_plottsne"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}raw_' tsnep_outdir='040_tsnemaps' path_metadata='${METADATA}/${file_metadata}'  path_rtsne_out='040_tsnemaps/${prefix_data}${prefix_panel}${prefix_pca}raw_rtsne_out.rda' path_rtsne_data='040_tsnemaps/${prefix_data}${prefix_panel}${prefix_pca}raw_rtsne_data.xls' path_clustering='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls' path_clustering_labels='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering_labels.xls'  pdf_width=22 pdf_height=7 tsne_distse=1 tsne_quantse=0.9" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
  tail $ROUT/03_plottsne.Rout

fi

### Plot tSNE with different perplexity

for perplexity in 50 100 200 1000
do

  if ${plottsne}; then
    echo ">>> 03_plottsne"
    R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}raw_perp${perplexity}_' tsnep_outdir='040_tsnemaps' path_metadata='${METADATA}/${file_metadata}'  path_rtsne_out='040_tsnemaps/${prefix_data}${prefix_panel}${prefix_pca}raw_perp${perplexity}_rtsne_out.rda' path_rtsne_data='040_tsnemaps/${prefix_data}${prefix_panel}${prefix_pca}raw_perp${perplexity}_rtsne_data.xls' path_clustering='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls' path_clustering_labels='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering_labels.xls'  pdf_width=22 pdf_height=7 tsne_distse=1 tsne_quantse=0.9" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
    tail $ROUT/03_plottsne.Rout

  fi

done


### Get cluster frequencies
if ${frequencies}; then
  echo ">>> 04_frequencies"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}' freq_outdir='050_frequencies' path_metadata='${METADATA}/${file_metadata}'  path_clustering='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls' path_clustering_labels='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering_labels.xls' path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_1dataset_3responses.R'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
  tail $ROUT/04_frequencies.Rout
fi


### Explore the expression of markers
if ${expression}; then
  echo ">>> 04_expression"

  # based on raw data
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' expr_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}raw_' expr_outdir='080_expression' path_data='010_data/${prefix_data}${prefix_panel}expr_raw.rds' path_metadata='${METADATA}/${file_metadata}'  path_clustering_observables='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' path_clustering='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls'  path_clustering_labels='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering_labels.xls'  path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_1dataset_3responses.R' analysis_type='all'" $RCODE/04_expression.R $ROUT/04_expression.Rout
  tail $ROUT/04_expression.Rout

  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' expr_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}raw_' expr_outdir='080_expression' path_data='010_data/${prefix_data}${prefix_panel}expr_raw.rds' path_metadata='${METADATA}/${file_metadata}'  path_clustering_observables='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' path_clustering='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls'  path_clustering_labels='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering_labels.xls'  path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_1dataset_3responses.R' analysis_type='clust'" $RCODE/04_expression.R $ROUT/04_expression.Rout
  tail $ROUT/04_expression.Rout

fi












#
