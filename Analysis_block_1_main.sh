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
    --data_normalization)
    data_normalization=${2}
    shift
    ;;
    --pcascores)
    pcascores=${2}
    shift
    ;;
    --select_observables)
    select_observables=${2}
    shift
    ;;
    --flowsom)
    flowsom=${2}
    shift
    ;;
    --flowsom_validation)
    flowsom_validation=${2}
    shift
    ;;
    --heatmaps)
    heatmaps=${2}
    shift
    ;;
    --runtsne)
    runtsne=${2}
    shift
    ;;
    --plottsne)
    plottsne=${2}
    shift
    ;;
    --plottsne_expr)
    plottsne_expr=${2}
    shift
    ;;
    --frequencies)
    frequencies=${2}
    shift
    ;;
    --expression)
    expression=${2}
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
    --prefix_pca)
    prefix_pca=${2}
    shift
    ;;
    --prefix_clust)
    prefix_clust=${2}
    shift
    ;;
    --pca_score_cutoff)
    pca_score_cutoff=${2}
    shift
    ;;
    --rand_seed_consensus)
    rand_seed_consensus=${2}
    shift
    ;;
    --nmetaclusts)
    nmetaclusts=${2}
    shift
    ;;
    --tsne_pmin)
    tsne_pmin=${2}
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

RWD=$RWD_MAIN/${data_dir}
ROUT=$RWD/Rout
mkdir -p $ROUT

if ${data_normalization} || ${pcascores} || ${select_observables} || ${flowsom} || ${flowsom_validation} || ${heatmaps} || ${runtsne} || ${plottsne} || ${plottsne_expr} || ${frequencies} || ${expression}
then
  echo "$RWD"
  echo ">>> reclustering into ${nmetaclusts} groups"
fi

### Data normalization
if ${data_normalization}; then
  echo ">>> 01_data_normalization"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' data_prefix='${prefix_data}${prefix_panel}' data_outdir='010_data' path_metadata='${METADATA}/${file_metadata}' path_panel='${PANELS}/${file_panel}'" $RCODE/01_data_normalization.R $ROUT/01_data_normalization.Rout
  tail $ROUT/01_data_normalization.Rout
fi

### PCA scores based on raw data
if ${pcascores}; then
  echo ">>> 01_pcascores"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' pcas_prefix='${prefix_data}${prefix_panel}' pcas_outdir='020_pcascores' path_data='010_data/${prefix_data}${prefix_panel}expr_raw.rds' path_metadata='${METADATA}/${file_metadata}' path_panel='${PANELS}/${file_panel}'" $RCODE/01_pcascores.R $ROUT/01_pcascores.Rout
  tail $ROUT/01_pcascores.Rout
fi

### Select observables for clustering
if ${select_observables}; then
  echo ">>> 02_select_observables"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' observ_prefix='${prefix_data}${prefix_panel}${prefix_pca}' observ_outdir='030_heatmaps' path_pca_score='020_pcascores/${prefix_data}${prefix_panel}princompscore_by_sample.xls' pca_score_cutoff=${pca_score_cutoff} pca_skip_top=0" $RCODE/02_select_observables.R $ROUT/02_select_observables.Rout
  tail $ROUT/02_select_observables.Rout
fi

### FlowSOM clustering based on raw data - average linkage
if ${flowsom}; then
  echo ">>> 02_flowsom"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' flowsom_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}' flowsom_outdir='030_heatmaps' path_data='010_data/${prefix_data}${prefix_panel}expr_raw.rds' path_clustering_observables='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' nmetaclusts=${nmetaclusts} rand_seed_consensus=${rand_seed_consensus}" $RCODE/02_flowsom.R $ROUT/02_flowsom.Rout
  tail $ROUT/02_flowsom.Rout

  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' codes_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}' codes_outdir='030_heatmaps' path_codes='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}codes.xls' pdf_width=22 pdf_height=7" $RCODE/02_som_codes.R $ROUT/02_som_codes.Rout
  tail $ROUT/02_som_codes.Rout
fi

### FlowSOM validation
if ${flowsom_validation}; then
  echo ">>> 02_flowsom_validation"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' flowsom_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}' flowsom_outdir='030_flowsom_validation' path_data='010_data/${prefix_data}${prefix_panel}expr_raw.rds' path_clustering_observables='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' nmetaclusts=${nmetaclusts} rand_seed_consensus=${rand_seed_consensus}" $RCODE/02_flowsom_validation.R $ROUT/02_flowsom_validation.Rout
  tail $ROUT/02_flowsom_validation.Rout
fi

### Heatmaps
if ${heatmaps}; then
  echo ">>> 02_heatmaps"
  # based on raw data
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' heatmap_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}raw_' heatmap_outdir='030_heatmaps' path_data='010_data/${prefix_data}${prefix_panel}expr_raw.rds'    path_clustering_observables='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' path_clustering='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering.xls'  path_clustering_labels='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering_labels.xls' path_marker_selection='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}marker_selection.txt' aggregate_fun='median' pheatmap_palette='YlGnBu' pheatmap_palette_rev=FALSE pheatmap_scale=TRUE extra_path_data='010_data/${prefix_data}${prefix_panel}expr_norm.rds' extra_pheatmap_palette='RdYlBu' extra_pheatmap_palette_rev=TRUE extra_suffix='_norm'" $RCODE/02_heatmaps.R $ROUT/02_heatmaps.Rout
  tail $ROUT/02_heatmaps.Rout
fi


### Run tSNE
if ${runtsne}; then
  echo ">>> 03_runtsne"
  # on raw data
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsne_prefix='${prefix_data}${prefix_panel}${prefix_pca}raw_' tsne_outdir='040_tsnemaps' path_data='010_data/${prefix_data}${prefix_panel}expr_raw.rds' path_metadata='${METADATA}/${file_metadata}'  path_clustering_observables='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' tsne_pmin=${tsne_pmin}" $RCODE/03_runtsne.R $ROUT/03_runtsne.Rout
  tail $ROUT/03_runtsne.Rout
fi

### Plot tSNE
if ${plottsne}; then
  echo ">>> 03_plottsne"
  ### Based on raw data
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}raw_' tsnep_outdir='040_tsnemaps' path_metadata='${METADATA}/${file_metadata}'  path_rtsne_out='040_tsnemaps/${prefix_data}${prefix_panel}${prefix_pca}raw_rtsne_out.rda' path_rtsne_data='040_tsnemaps/${prefix_data}${prefix_panel}${prefix_pca}raw_rtsne_data.xls' path_clustering='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering.xls' path_clustering_labels='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering_labels.xls'  pdf_width=22 pdf_height=7 tsne_distse=1 tsne_quantse=0.9" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
  tail $ROUT/03_plottsne.Rout
fi


### Run and plot tSNE with different perplexity

# for perplexity in 50 100 200 1000
# do
#
#   if ${runtsne}; then
#     echo ">>> 03_runtsne"
#     R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsne_prefix='${prefix_data}${prefix_panel}${prefix_pca}raw_perp${perplexity}_' tsne_outdir='040_tsnemaps' path_data='010_data/${prefix_data}${prefix_panel}expr_raw.rds' path_metadata='${METADATA}/${file_metadata}'  path_clustering_observables='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' tsne_pmin=${tsne_pmin} perplexity=${perplexity}" $RCODE/03_runtsne.R $ROUT/03_runtsne.Rout
#     tail $ROUT/03_runtsne.Rout
#   fi
#
#   if ${plottsne}; then
#     echo ">>> 03_plottsne"
#     R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}raw_perp${perplexity}_' tsnep_outdir='040_tsnemaps' path_metadata='${METADATA}/${file_metadata}'  path_rtsne_out='040_tsnemaps/${prefix_data}${prefix_panel}${prefix_pca}raw_perp${perplexity}_rtsne_out.rda' path_rtsne_data='040_tsnemaps/${prefix_data}${prefix_panel}${prefix_pca}raw_perp${perplexity}_rtsne_data.xls' path_clustering='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering.xls' path_clustering_labels='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering_labels.xls'  pdf_width=22 pdf_height=7 tsne_distse=1 tsne_quantse=0.9" $RCODE/03_plottsne.R $ROUT/03_plottsne.Rout
#     tail $ROUT/03_plottsne.Rout
#   fi
#
# done


### Plot tSNE with marker expression as a heat
if ${plottsne_expr}; then
  echo ">>> 03_plottsne_expr"

  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' tsnep_prefix='${prefix_data}${prefix_panel}${prefix_pca}raw_mark_raw_' tsnep_outdir='040_tsnemaps_expr' path_data='010_data/${prefix_data}${prefix_panel}expr_raw.rds' path_metadata='${METADATA}/${file_metadata}'  path_rtsne_out='040_tsnemaps/${prefix_data}${prefix_panel}${prefix_pca}raw_rtsne_out.rda' path_rtsne_data='040_tsnemaps/${prefix_data}${prefix_panel}${prefix_pca}raw_rtsne_data.xls' path_clustering_observables='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' pdf_width=22 pdf_height=7" $RCODE/03_plottsne_expr.R $ROUT/03_plottsne_expr.Rout
  tail $ROUT/03_plottsne_expr.Rout

fi

### Get cluster frequencies
if ${frequencies}; then
  echo ">>> 04_frequencies"
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' freq_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}' freq_outdir='050_frequencies_auto' path_metadata='${METADATA}/${file_metadata}'  path_clustering='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering.xls' path_clustering_labels='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering_labels.xls' path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_1dataset_3responses.R'" $RCODE/04_frequencies.R $ROUT/04_frequencies.Rout
  tail $ROUT/04_frequencies.Rout
fi


### Explore the expression of markers
if ${expression}; then
  echo ">>> 04_expression"

  ## Based on raw data
  # For all the cells
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' expr_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}raw_' expr_outdir='080_expression' path_data='010_data/${prefix_data}${prefix_panel}expr_raw.rds' path_metadata='${METADATA}/${file_metadata}'  path_clustering_observables='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' path_clustering='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering.xls'  path_clustering_labels='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering_labels.xls'  path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_1dataset_3responses.R' analysis_type='all'" $RCODE/04_expression.R $ROUT/04_expression.Rout
  tail $ROUT/04_expression.Rout

  # Per cluster
  R CMD BATCH --no-save --no-restore "--args rwd='$RWD' expr_prefix='${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}raw_' expr_outdir='080_expression' path_data='010_data/${prefix_data}${prefix_panel}expr_raw.rds' path_metadata='${METADATA}/${file_metadata}'  path_clustering_observables='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' path_clustering='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering.xls'  path_clustering_labels='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}clustering_labels.xls'  path_fun_models='$RCODE/00_models.R' path_fun_formulas='$RCODE/00_formulas_1dataset_3responses.R' analysis_type='clust'" $RCODE/04_expression.R $ROUT/04_expression.Rout
  tail $ROUT/04_expression.Rout

fi




















#
