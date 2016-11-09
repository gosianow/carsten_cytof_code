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
    --cluster_extracting)
    cluster_extracting=${2}
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
    --extract_cluster)
    extract_cluster=${2}
    shift
    ;;
    --extract_dir)
    extract_dir=${2}
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

### Cluster extracting
if ${cluster_extracting}; then
  echo "02_cluster_extracting"

  RWD=$RWD_MAIN/${data_dir}
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"

  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' extract_outdir='${RWD_MAIN}/${extract_dir}/010_cleanfcs'  path_metadata='${METADATA}/${file_metadata}'  path_clustering='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls' path_clustering_labels='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering_labels.xls' extract_cluster=${extract_cluster}" $RCODE/02_cluster_extracting.R $ROUT/02_cluster_extracting.Rout
  tail $ROUT/02_cluster_extracting.Rout


fi











#
