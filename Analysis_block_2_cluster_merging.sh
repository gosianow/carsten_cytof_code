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
echo "$RWD"

### Cluster merging
if ${cluster_merging}; then
  echo "02_cluster_merging"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' merging_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}' path_cluster_merging='${file_merging}' path_clustering='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering.xls'" $RCODE/02_cluster_merging.R $ROUT/02_cluster_merging.Rout
  tail $ROUT/02_cluster_merging.Rout
fi

if [ ! -e "$RWD/030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls" ]; then
  echo "File '$RWD/030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls' does NOT exist!"
  exit
fi

### Heatmaps
if ${heatmaps}; then
  echo "02_heatmaps"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' heatmap_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}'  path_clustering_observables='${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' path_clustering='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls'  path_clustering_labels='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering_labels.xls' path_marker_selection='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}marker_selection.txt'" $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
  tail $ROUT/02_heatmaps.Rout
fi

### Plot tSNE
if ${plottsne}; then
  echo "03_plottsne"

  ### Based on 0-1 normalized data
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' tsnep_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}norm_' path_rtsne_out='${prefix_data}${prefix_panel}${prefix_pca}rtsne_out_norm.rda' path_rtsne_data='${prefix_data}${prefix_panel}${prefix_pca}rtsne_data_norm.xls' path_clustering='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls' path_clustering_labels='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10 tsnep_suffix=''" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
  tail $ROUT/03_plottsne.Rout

  ### Based on raw data
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' tsnep_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}raw_' path_rtsne_out='${prefix_data}${prefix_panel}${prefix_pca}rtsne_out_raw.rda' path_rtsne_data='${prefix_data}${prefix_panel}${prefix_pca}rtsne_data_raw.xls' path_clustering='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls' path_clustering_labels='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering_labels.xls'  tsne_cmin=1000 pdf_width=15 pdf_height=10 tsnep_suffix=''" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
  tail $ROUT/03_plottsne.Rout

fi

### Get cluster frequencies
if ${frequencies}; then
  echo "04_frequencies"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}' freq_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}' path_clustering='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls' path_clustering_labels='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering_labels.xls'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
  tail $ROUT/04_frequencies.Rout
fi











#
