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
    --new_data_dir)
    new_data_dir=${2}
    shift
    ;;
    --cytokines_merging)
    cytokines_merging=${2}
    shift
    ;;
    --data_name1)
    data_name1=${2}
    shift
    ;;
    --data_name2)
    data_name2=${2}
    shift
    ;;
    --data_dir1)
    data_dir1=${2}
    shift
    ;;
    --data_dir2)
    data_dir2=${2}
    shift
    ;;
    --file_metadata1)
    file_metadata1=${2}
    shift
    ;;
    --file_metadata2)
    file_metadata2=${2}
    shift
    ;;
    --prefix_data1)
    prefix_data1=${2}
    shift
    ;;
    --prefix_data2)
    prefix_data2=${2}
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
    --prefix_merging1)
    prefix_merging1=${2}
    shift
    ;;
    --prefix_merging2)
    prefix_merging2=${2}
    shift
    ;;
    --prefix_data_merging)
    prefix_data_merging=${2}
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
    --prefix_merging_cyt1)
    prefix_merging_cyt1=${2}
    shift
    ;;
    --prefix_merging_cyt2)
    prefix_merging_cyt2=${2}
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

RWD=$RWD_MAIN/${new_data_dir}
ROUT=$RWD/Rout
mkdir -p $ROUT

if ${cytokines_merging}
then
  echo "$RWD"
fi

if ${cytokines_merging}; then

  ### Get cluster frequencies - 3 responses
  echo "08_frequencies_merged"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${prefix_data_merging}' freq_outdir='08_cytokines_merged' path_metadata=c('${METADATA}/${file_metadata1}','${METADATA}/${file_metadata2}')  path_counts=c('${RWD_MAIN}/${data_dir1}/060_cytokines_bimatrix/03_frequencies/${prefix_data1}${prefix_panel}${prefix_pca}${prefix_merging1}${prefix_clsubset}${prefix_cytokines_cutoffs}raw2_${prefix_merging_cyt1}counts.xls','${RWD_MAIN}/${data_dir2}/060_cytokines_bimatrix/03_frequencies/${prefix_data2}${prefix_panel}${prefix_pca}${prefix_merging2}${prefix_clsubset}${prefix_cytokines_cutoffs}raw2_${prefix_merging_cyt2}counts.xls') data_name=c('${data_name1}','${data_name2}') path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_2datasets_3responses.R'" $RCODE/08_frequencies_merged.R $ROUT/08_frequencies_merged.Rout
  tail $ROUT/08_frequencies_merged.Rout

  ### Get cluster frequencies - 2 responses
  echo "08_frequencies_merged"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${prefix_data_merging}' freq_outdir='08_cytokines_merged_2responses' path_metadata=c('${METADATA}/${file_metadata1}','${METADATA}/${file_metadata2}')  path_counts=c('${RWD_MAIN}/${data_dir1}/060_cytokines_bimatrix/03_frequencies/${prefix_data1}${prefix_panel}${prefix_pca}${prefix_merging1}${prefix_clsubset}${prefix_cytokines_cutoffs}raw2_${prefix_merging_cyt1}counts.xls','${RWD_MAIN}/${data_dir2}/060_cytokines_bimatrix/03_frequencies/${prefix_data2}${prefix_panel}${prefix_pca}${prefix_merging2}${prefix_clsubset}${prefix_cytokines_cutoffs}raw2_${prefix_merging_cyt2}counts.xls') data_name=c('${data_name1}','${data_name2}') path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_2datasets_2responses.R'" $RCODE/08_frequencies_merged.R $ROUT/08_frequencies_merged.Rout
  tail $ROUT/08_frequencies_merged.Rout

fi













#
