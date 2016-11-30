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
    --METADATA)
    METADATA=${2}
    shift
    ;;
    --PANELS)
    PANELS=${2}
    shift
    ;;
    --cytokines_fcs_saving)
    cytokines_fcs_saving=${2}
    shift
    ;;
    --new_data_dir)
    new_data_dir=${2}
    shift
    ;;
    --data_name)
    data_name=${2}
    shift
    ;;
    --data_dir)
    data_dir=${2}
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
    --prefix_data_merging)
    prefix_data_merging=${2}
    shift
    ;;
    --prefix_clust)
    prefix_clust=${2}
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

if ${cytokines_fcs_saving}; then

  RWD=$RWD_MAIN/${new_data_dir}
  ROUT=$RWD/Rout
  echo "$RWD"
  echo ">>> 06_cytokines_fcs_saving"

  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' save_prefix='${prefix_data_merging}${prefix_clust}${data_name}_' save_outdir='${outdir}/06_dumpfcs' path_panel='${PANELS}/${file_panel}' path_metadata='${METADATA}/${file_metadata}' path_fcs='$RWD_MAIN/${data_dir}/010_cleanfcs' path_clustering='${outdir}/01_clustering/${prefix_data_merging}${prefix_clust}${data_name}_clustering.xls' path_clustering_labels='${outdir}/01_clustering/${prefix_data_merging}${prefix_clust}${data_name}_clustering_labels.xls'" $RCODE/06_cytokines_fcs_saving.R $ROUT/06_cytokines_fcs_saving.Rout
  tail $ROUT/06_cytokines_fcs_saving.Rout

fi











#
