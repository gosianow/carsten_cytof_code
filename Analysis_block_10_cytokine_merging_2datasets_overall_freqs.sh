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
    --METADATA)
    METADATA=${2}
    shift
    ;;
    --new_data_dir)
    new_data_dir=${2}
    shift
    ;;
    --cytokines_merging_v4)
    cytokines_merging_v4=${2}
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

RWD=$RWD_MAIN/${new_data_dir}
ROUT=$RWD/Rout
mkdir -p $ROUT

if ${cytokines_merging_v4}
then
  echo "$RWD"
fi

if ${cytokines_merging_v4}; then

  FDR_cutoff=(0.05 0.1)
  suffix=('_top005' '_top01')

  for j in 0 1
  do

    ### Analysis of individual cytokine frequencies - 3 responses

    echo ">>> 08_frequencies_merged_bimatrix"

    R33 CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${prefix_data_merging}' freq_outdir='${outdir}/03_frequencies_auto' path_data=c('${RWD_MAIN}/${data_dir1}/01_clustering/${prefix_data1}${prefix_panel}${prefix_pca}${prefix_merging1}${prefix_clsubset}${prefix_cytokines_cutoffs}bimatrix.txt','${RWD_MAIN}/${data_dir2}/01_clustering/${prefix_data2}${prefix_panel}${prefix_pca}${prefix_merging2}${prefix_clsubset}${prefix_cytokines_cutoffs}bimatrix.txt') path_clustering_observables=c('${RWD_MAIN}/${data_dir1}/01_clustering/${prefix_data1}${prefix_panel}${prefix_pca}${prefix_merging1}${prefix_clsubset}${prefix_cytokines_cutoffs}clustering_observables.xls','${RWD_MAIN}/${data_dir2}/01_clustering/${prefix_data2}${prefix_panel}${prefix_pca}${prefix_merging2}${prefix_clsubset}${prefix_cytokines_cutoffs}clustering_observables.xls') path_metadata=c('${METADATA}/${file_metadata1}','${METADATA}/${file_metadata2}') data_name=c('${data_name1}','${data_name2}') path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_2datasets_3responses.R' path_fun_plot_frequencies='$RCODE/00_plot_frequencies.R' path_fun_plot_heatmaps='$RCODE/00_plot_heatmaps_for_sign_freqs.R'  pdf_hight=4 FDR_cutoff=${FDR_cutoff[$j]} suffix='${suffix[$j]}'" $RCODE/08_frequencies_merged_bimatrix.R $ROUT/08_frequencies_merged_bimatrix.Rout
    tail $ROUT/08_frequencies_merged_bimatrix.Rout


    ### Analysis of individual cytokine frequencies - 2 responses

    echo ">>> 08_frequencies_merged_bimatrix"

    R33 CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${prefix_data_merging}' freq_outdir='${outdir}/03_frequencies_auto_2responses' path_data=c('${RWD_MAIN}/${data_dir1}/01_clustering/${prefix_data1}${prefix_panel}${prefix_pca}${prefix_merging1}${prefix_clsubset}${prefix_cytokines_cutoffs}bimatrix.txt','${RWD_MAIN}/${data_dir2}/01_clustering/${prefix_data2}${prefix_panel}${prefix_pca}${prefix_merging2}${prefix_clsubset}${prefix_cytokines_cutoffs}bimatrix.txt') path_clustering_observables=c('${RWD_MAIN}/${data_dir1}/01_clustering/${prefix_data1}${prefix_panel}${prefix_pca}${prefix_merging1}${prefix_clsubset}${prefix_cytokines_cutoffs}clustering_observables.xls','${RWD_MAIN}/${data_dir2}/01_clustering/${prefix_data2}${prefix_panel}${prefix_pca}${prefix_merging2}${prefix_clsubset}${prefix_cytokines_cutoffs}clustering_observables.xls') path_metadata=c('${METADATA}/${file_metadata1}','${METADATA}/${file_metadata2}') data_name=c('${data_name1}','${data_name2}') path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_2datasets_2responses.R' path_fun_plot_frequencies='$RCODE/00_plot_frequencies.R' path_fun_plot_heatmaps='$RCODE/00_plot_heatmaps_for_sign_freqs.R' pdf_hight=4 FDR_cutoff=${FDR_cutoff[$j]} suffix='${suffix[$j]}'" $RCODE/08_frequencies_merged_bimatrix.R $ROUT/08_frequencies_merged_bimatrix.Rout
    tail $ROUT/08_frequencies_merged_bimatrix.Rout

  done


fi














#
