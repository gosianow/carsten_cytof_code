#!/bin/bash

shopt -s expand_aliases
source ~/.bash_aliases

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
    --cytokines_bimatrix)
    cytokines_bimatrix=${2}
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

if ${cytokines_bimatrix}; then

  RWD=$RWD_MAIN/${data_dir}
  ROUT=$RWD/Rout
  mkdir -p $ROUT
  echo "$RWD"

  ### Analysis based on different cutoffs for base and tx
  prefix_cytokines="${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}${prefix_clsubset}${prefix_cytokines_cutoffs}"

  ### Create the bimatrix
  echo ">>> 06_cytokines_bimatrix"
  R33 CMD BATCH --no-save --no-restore "--args rwd='$RWD' cytokines_prefix='${prefix_cytokines}' cytokines_outdir='${outdir}/01_clustering' path_data='010_data/${prefix_data}${prefix_panel}expr_raw.rds' path_metadata='${METADATA}/${file_metadata}' path_cytokines_cutoffs='${PANELS}/${file_cytokines_cutoffs}' path_clustering='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls' path_clustering_labels='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering_labels.xls' clsubset=${clsubset} cutoff_colname=c('positive_cutoff_raw_base','positive_cutoff_raw_tx')" $RCODE/06_cytokines_bimatrix.R $ROUT/06_cytokines_bimatrix.Rout
  tail $ROUT/06_cytokines_bimatrix.Rout


fi



























#
