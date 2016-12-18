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
    --new_data_dir)
    new_data_dir=${2}
    shift
    ;;
    --frequencies_merged)
    frequencies_merged=${2}
    shift
    ;;
    --expression_merged)
    expression_merged=${2}
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
    --data_name3)
    data_name3=${2}
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
    --data_dir3)
    data_dir3=${2}
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
    --file_metadata3)
    file_metadata3=${2}
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
    --prefix_data3)
    prefix_data3=${2}
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
    --prefix_merging3)
    prefix_merging3=${2}
    shift
    ;;
    --prefix_data_merging)
    prefix_data_merging=${2}
    shift
    ;;
    --METADATA)
    METADATA=${2}
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

if ${frequencies_merged} || ${expression_merged} || ${runtsne_merged} || ${plottsne_merged}
then
  echo "$RWD"
fi


### Get cluster frequencies
if ${frequencies_merged}; then
  echo "08_frequencies_merged"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${prefix_data_merging}' freq_outdir='08_frequencies_merged' path_metadata=c('${METADATA}/${file_metadata1}','${METADATA}/${file_metadata2}','${METADATA}/${file_metadata3}')  path_counts=c('${RWD_MAIN}/${data_dir1}/050_frequencies/${prefix_data1}${prefix_panel}${prefix_pca}${prefix_merging1}counts.xls','${RWD_MAIN}/${data_dir2}/050_frequencies/${prefix_data2}${prefix_panel}${prefix_pca}${prefix_merging2}counts.xls','${RWD_MAIN}/${data_dir3}/050_frequencies/${prefix_data3}${prefix_panel}${prefix_pca}${prefix_merging3}counts.xls') data_name=c('${data_name1}','${data_name2}','${data_name3}') path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_3datasets_3responses.R' path_fun_plot_frequencies='$RCODE/00_plot_frequencies.R' path_fun_plot_heatmaps='$RCODE/00_plot_heatmaps_for_sign_freqs.R'" $RCODE/08_frequencies_merged.R $ROUT/08_frequencies_merged.Rout
  tail $ROUT/08_frequencies_merged.Rout
fi

### Marker expression analysis
if ${expression_merged}; then
  echo "08_expression_merged"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' expr_prefix='${prefix_data_merging}' expr_outdir='08_expression_merged' path_metadata=c('${METADATA}/${file_metadata1}','${METADATA}/${file_metadata2}','${METADATA}/${file_metadata3}')  path_expression=c('${RWD_MAIN}/${data_dir1}/080_expression/${prefix_data1}${prefix_panel}${prefix_pca}${prefix_merging1}raw_expr_all.xls','${RWD_MAIN}/${data_dir2}/080_expression/${prefix_data2}${prefix_panel}${prefix_pca}${prefix_merging2}raw_expr_all.xls','${RWD_MAIN}/${data_dir3}/080_expression/${prefix_data3}${prefix_panel}${prefix_pca}${prefix_merging3}raw_expr_all.xls') data_name=c('${data_name1}','${data_name2}','${data_name3}') path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_3datasets_3responses.R' path_fun_plot_heatmaps='$RCODE/00_plot_heatmaps_for_sign_expr.R' path_fun_plot_expression='$RCODE/00_plot_expression.R' analysis_type='all'" $RCODE/08_expression_merged.R $ROUT/08_expression_merged.Rout
  tail $ROUT/08_expression_merged.Rout
fi

















#
