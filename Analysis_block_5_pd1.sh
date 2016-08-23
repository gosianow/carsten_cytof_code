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
    --pd1)
    pd1=${2}
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
    --path_cytokines_cutoffs)
    path_cytokines_cutoffs=${2}
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

### PD-1 analysis
if ${pd1}; then
  
    RWD=$RWD_MAIN/${data_dir}
    ROUT=$RWD/Rout
    mkdir -p $ROUT
    echo "$RWD"

    echo "07_pd1"

    R CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_metadata='${METADATA}/${file_metadata}'  pd1_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}${prefix_clsubset}${prefix_cytokines_cutoffs}raw_'  path_cytokines_cutoffs='${PANELS}/${path_cytokines_cutoffs}' path_clustering='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls' path_clustering_labels='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering_labels.xls' clsubset=${clsubset} cutoff_colname='positive_cutoff_raw_base' data2analyse='raw' nmetaclusts=${nmetaclusts}  path_marker_selection='${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}${prefix_clsubset}marker_selection.txt' path_rtsne_out='${prefix_data}${prefix_panel}${prefix_pca}rtsne_out_norm.rda' path_rtsne_data='${prefix_data}${prefix_panel}${prefix_pca}rtsne_data_norm.xls' pdf_width=15 pdf_height=10 tsnep_suffix='_norm'" $RCODE/07_pd1.R $ROUT/07_pd1.Rout
    tail $ROUT/07_pd1.Rout

  done

fi











#
