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

  prefix_cytokines="${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}${prefix_clsubset}${prefix_cytokines_cutoffs}"

  ### Get the expression of cytokines for the Tmem cluster
  echo "06_cytokines_expression"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' cytokines_prefix='${prefix_cytokines}' cytokines_outdir='060_cytokines_expression' path_data='010_data/${prefix_data}${prefix_panel}expr_raw.rds' path_cytokines_cutoffs='${PANELS}/${file_cytokines_cutoffs}' path_clustering='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering.xls' path_clustering_labels='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_merging}clustering_labels.xls' clsubset=${clsubset} cutoff_colname='positive_cutoff_raw_base'" $RCODE/06_cytokines_expression.R $ROUT/06_cytokines_expression.Rout
  tail $ROUT/06_cytokines_expression.Rout


  ### Expression analysis
  echo "04_expression"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' expr_prefix='${prefix_cytokines}' expr_outdir='060_cytokines_expression' path_data='060_cytokines_expression/${prefix_cytokines}expr_raw.rds' path_metadata='${METADATA}/${file_metadata}'  path_clustering_observables='060_cytokines_expression/${prefix_cytokines}clustering_observables.xls' path_clustering='060_cytokines_expression/${prefix_cytokines}clustering.xls'  path_clustering_labels='060_cytokines_expression/${prefix_cytokines}clustering_labels.xls'  path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_1dataset_3responses.R' analysis_type='all'" $RCODE/04_expression.R $ROUT/04_expression.Rout
  tail $ROUT/04_expression.Rout


fi
































#
