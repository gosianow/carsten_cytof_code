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
    --runtsne_merged)
    runtsne_merged=${2}
    shift
    ;;
    --plottsne_merged)
    plottsne_merged=${2}
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


### Get cluster frequencies - models with NR, R and HD
if ${frequencies_merged}; then
  echo ">>> 08_frequencies_merged"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${prefix_data_merging}' freq_outdir='08_frequencies_merged' path_metadata=c('${METADATA}/${file_metadata1}','${METADATA}/${file_metadata2}')  path_counts=c('${RWD_MAIN}/${data_dir1}/050_frequencies/${prefix_data1}${prefix_panel}${prefix_pca}${prefix_merging1}counts.xls','${RWD_MAIN}/${data_dir2}/050_frequencies/${prefix_data2}${prefix_panel}${prefix_pca}${prefix_merging2}counts.xls') data_name=c('${data_name1}','${data_name2}') path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_2datasets_3responses.R'" $RCODE/08_frequencies_merged.R $ROUT/08_frequencies_merged.Rout
  tail $ROUT/08_frequencies_merged.Rout
fi

### Get cluster frequencies - models with NR, R
if ${frequencies_merged}; then
  echo ">>> 08_frequencies_merged"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${prefix_data_merging}' freq_outdir='08_frequencies_merged_2responses' path_metadata=c('${METADATA}/${file_metadata1}','${METADATA}/${file_metadata2}')  path_counts=c('${RWD_MAIN}/${data_dir1}/050_frequencies/${prefix_data1}${prefix_panel}${prefix_pca}${prefix_merging1}counts.xls','${RWD_MAIN}/${data_dir2}/050_frequencies/${prefix_data2}${prefix_panel}${prefix_pca}${prefix_merging2}counts.xls') data_name=c('${data_name1}','${data_name2}') path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_2datasets_2responses.R'" $RCODE/08_frequencies_merged.R $ROUT/08_frequencies_merged.Rout
  tail $ROUT/08_frequencies_merged.Rout
fi


### Marker expression analysis - models with NR, R and HD
if ${expression_merged}; then
  echo ">>> 08_expression_merged"
  # all
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' expr_prefix='${prefix_data_merging}' expr_outdir='08_expression_merged' path_metadata=c('${METADATA}/${file_metadata1}','${METADATA}/${file_metadata2}')  path_expression=c('${RWD_MAIN}/${data_dir1}/080_expression/${prefix_data1}${prefix_panel}${prefix_pca}${prefix_merging1}raw_expr_all.xls','${RWD_MAIN}/${data_dir2}/080_expression/${prefix_data2}${prefix_panel}${prefix_pca}${prefix_merging2}raw_expr_all.xls') data_name=c('${data_name1}','${data_name2}') path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_2datasets_3responses.R' analysis_type='all'" $RCODE/08_expression_merged.R $ROUT/08_expression_merged.Rout
  tail $ROUT/08_expression_merged.Rout
  # clust
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' expr_prefix='${prefix_data_merging}' expr_outdir='08_expression_merged' path_metadata=c('${METADATA}/${file_metadata1}','${METADATA}/${file_metadata2}')  path_expression=c('${RWD_MAIN}/${data_dir1}/080_expression/${prefix_data1}${prefix_panel}${prefix_pca}${prefix_merging1}raw_expr_clust.xls','${RWD_MAIN}/${data_dir2}/080_expression/${prefix_data2}${prefix_panel}${prefix_pca}${prefix_merging2}raw_expr_clust.xls') data_name=c('${data_name1}','${data_name2}') path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_2datasets_3responses.R' analysis_type='clust'" $RCODE/08_expression_merged.R $ROUT/08_expression_merged.Rout
  tail $ROUT/08_expression_merged.Rout
fi


### Marker expression analysis - models with NR, R
if ${expression_merged}; then
  echo ">>> 08_expression_merged"
  # all
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' expr_prefix='${prefix_data_merging}' expr_outdir='08_expression_merged_2responses' path_metadata=c('${METADATA}/${file_metadata1}','${METADATA}/${file_metadata2}')  path_expression=c('${RWD_MAIN}/${data_dir1}/080_expression/${prefix_data1}${prefix_panel}${prefix_pca}${prefix_merging1}raw_expr_all.xls','${RWD_MAIN}/${data_dir2}/080_expression/${prefix_data2}${prefix_panel}${prefix_pca}${prefix_merging2}raw_expr_all.xls') data_name=c('${data_name1}','${data_name2}') path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_2datasets_2responses.R' analysis_type='all'" $RCODE/08_expression_merged.R $ROUT/08_expression_merged.Rout
  tail $ROUT/08_expression_merged.Rout
  # clust
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' expr_prefix='${prefix_data_merging}' expr_outdir='08_expression_merged_2responses' path_metadata=c('${METADATA}/${file_metadata1}','${METADATA}/${file_metadata2}')  path_expression=c('${RWD_MAIN}/${data_dir1}/080_expression/${prefix_data1}${prefix_panel}${prefix_pca}${prefix_merging1}raw_expr_clust.xls','${RWD_MAIN}/${data_dir2}/080_expression/${prefix_data2}${prefix_panel}${prefix_pca}${prefix_merging2}raw_expr_clust.xls') data_name=c('${data_name1}','${data_name2}') path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_2datasets_2responses.R' analysis_type='clust'" $RCODE/08_expression_merged.R $ROUT/08_expression_merged.Rout
  tail $ROUT/08_expression_merged.Rout
fi

# if ${runtsne_merged}; then
#   echo "08_runtsne_merged"
#   R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsne_prefix='${prefix_data_merging}' tsne_outdir='08_tsnemaps_merged' path_metadata=c('${METADATA}/${file_metadata1}','${METADATA}/${file_metadata2}')  path_data=c('${RWD_MAIN}/${data_dir1}/010_data/${prefix_data1}${prefix_panel}expr_raw.rds','${RWD_MAIN}/${data_dir2}/010_data/${prefix_data2}${prefix_panel}expr_raw.rds') path_clustering_observables=c('${RWD_MAIN}/${data_dir1}/030_heatmaps/${prefix_data1}${prefix_panel}${prefix_pca}clustering_observables.xls','${RWD_MAIN}/${data_dir2}/030_heatmaps/${prefix_data2}${prefix_panel}${prefix_pca}clustering_observables.xls') data_name=c('${data_name1}','${data_name2}') tsne_pmin=1500" $RCODE/08_runtsne_merged.R $ROUT/08_runtsne_merged.Rout
#   tail $ROUT/08_runtsne_merged.Rout
# fi
#
#
# if ${plottsne_merged}; then
#   echo "08_plottsne_merged"
#   R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='${prefix_data_merging}' tsnep_outdir='08_tsnemaps_merged' path_rtsne_out='08_tsnemaps_merged/${prefix_data_merging}rtsne_out.rda' path_rtsne_data='08_tsnemaps_merged/${prefix_data_merging}rtsne_data.xls'  path_metadata=c('${METADATA}/${file_metadata1}','${METADATA}/${file_metadata2}')  path_clustering=c('${RWD_MAIN}/${data_dir1}/030_heatmaps/${prefix_data1}${prefix_panel}${prefix_pca}${prefix_merging1}clustering.xls','${RWD_MAIN}/${data_dir2}/030_heatmaps/${prefix_data2}${prefix_panel}${prefix_pca}${prefix_merging2}clustering.xls') path_clustering_labels=c('${RWD_MAIN}/${data_dir1}/030_heatmaps/${prefix_data1}${prefix_panel}${prefix_pca}${prefix_merging1}clustering_labels.xls','${RWD_MAIN}/${data_dir2}/030_heatmaps/${prefix_data2}${prefix_panel}${prefix_pca}${prefix_merging2}clustering_labels.xls') data_name=c('${data_name1}','${data_name2}') pdf_width=15 pdf_height=10" $RCODE/08_plottsne_merged.R $ROUT/08_plottsne_merged.Rout
#   tail $ROUT/08_plottsne_merged.Rout
# fi

















#
