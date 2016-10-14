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
    --fcs_saving)
    fcs_saving=${2}
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

if ${fcs_saving}; then

  RWD=$RWD_MAIN/${data_dir}
  ROUT=$RWD/Rout
  echo "$RWD"
  echo "02_fcs_saving"

  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' save_prefix='${prefix_data}${prefix_panel}' save_outdir='060_dumpfcs' path_metadata='${METADATA}/${file_metadata}' path_panel='${PANELS}/${file_panel}'" $RCODE/02_fcs_saving.R $ROUT/02_fcs_saving.Rout
  tail $ROUT/02_fcs_saving.Rout

fi











#
