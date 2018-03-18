#!/bin/bash

###############################################################################################################
## Define paths to software and reference files

RWD_MAIN=../../carsten_cytof/PD1_project
MAIN_FIGURE_DIR=${RWD_MAIN}/CK_figures_changed

mkdir -p MAIN_FIGURE_DIR

###############################################################################################################


#######################################
### Figure 5
#######################################


DIR_FIGURE=${MAIN_FIGURE_DIR}/Figure5
mkdir -p ${DIR_FIGURE}

# Figure 5a - Frequencies

rsync -v -r -t ${RWD_MAIN}/CK_2016-06-merged_23_29/03all/050_frequencies/3responses_both/03v2_23allmerging5_29all3merging1_frequencies_plot_base3.pdf ${DIR_FIGURE}/

rsync -v -r -t ${RWD_MAIN}/CK_2016-06-merged_23_29/03all/050_frequencies/3responses_both/03v2_23allmerging5_29all3merging1_frequencies_plot_tx3.pdf ${DIR_FIGURE}/


# Figure 5c - Heatmaps

rsync -v -r -t ${RWD_MAIN}/CK_2016-06-merged_23_29/03/080_expression/3responses_both/03_23merging5_29merging2_all_expr_lmer_interglht_pheatmap2_top10_ex.pdf ${DIR_FIGURE}/


#######################################
### Figure 6
#######################################


DIR_FIGURE=${MAIN_FIGURE_DIR}/Figure6
mkdir -p ${DIR_FIGURE}

# Figure 6a - t-SNE

rsync -v -r -t ${RWD_MAIN}/CK_2016-06-23_03/040_tsnemaps/23_03_pca1_merging5_tSNEone.pdf ${DIR_FIGURE}/

rsync -v -r -t ${RWD_MAIN}/CK_2016-06-23_03/030_heatmaps/23_03_pca1_merging5_pheatmap_all_sel_no_clust_scale.pdf ${DIR_FIGURE}/

# P-values

rsync -v -r -t ${RWD_MAIN}/CK_2016-06-merged_23_29/03/050_frequencies/3responses_both/03_23merging5_29merging2_frequencies_glmer_binomial_interglht_pheatmap3pvs_top10.pdf ${DIR_FIGURE}/



#######################################
### Supplementary Figure 7
#######################################


DIR_FIGURE=${MAIN_FIGURE_DIR}/SuppFigure7
mkdir -p ${DIR_FIGURE}


rsync -v -r -t ${RWD_MAIN}/CK_2016-06-merged_23_29/03all/050_frequencies/3responses_both/03v2_23allmerging5_29all3merging1_frequencies_plot_base4.pdf ${DIR_FIGURE}/

rsync -v -r -t ${RWD_MAIN}/CK_2016-06-merged_23_29/03all/050_frequencies/3responses_both/03v2_23allmerging5_29all3merging1_frequencies_plot_tx4.pdf ${DIR_FIGURE}/


# P-values

rsync -v -r -t ${RWD_MAIN}/CK_2016-06-merged_23_29/03all/050_frequencies/3responses_both/03v2_23allmerging5_29all3merging1_frequencies_glmer_binomial_interglht_pheatmap3pvs_top10.pdf ${DIR_FIGURE}/


#######################################
### Supplementary Figure 9
### CD14+ at higher resolution 
#######################################


DIR_FIGURE=${MAIN_FIGURE_DIR}/SuppFigure9
mkdir -p ${DIR_FIGURE}


# Heatmap 

rsync -v -r -t ${RWD_MAIN}/CK_2016-06-23_03/030_codes/23_03_pca1_cl20_merging5_ComplexHeatmap_codes_sel_no_clust_scale_NRvsR_baseANDtx.pdf ${DIR_FIGURE}/


# File with frequencies so Carsten can redo the plot of frequencies 

rsync -v -r -t ${RWD_MAIN}/CK_2016-06-23_03/050_frequencies_codes/23_03_pca1_cl20_frequencies.xls ${DIR_FIGURE}/


#######################################
### Supplementary Figure 11
### CellCnn 
#######################################


DIR_FIGURE=${MAIN_FIGURE_DIR}/SuppFigure11
mkdir -p ${DIR_FIGURE}

rsync -v -r -t ${RWD_MAIN}/cellcnn_tsne/23all_03v2_pca1_merging5_data23_tSNEgroup_filter_0_type2.pdf ${DIR_FIGURE}/



#######################################
### Supplementary Figure 17
### Cox regression
#######################################


DIR_FIGURE=${MAIN_FIGURE_DIR}/SuppFigure17
mkdir -p ${DIR_FIGURE}

# Figure 17 A

rsync -v -r -t ${RWD_MAIN}/PD-1_PatientClinicalData/ck_analysis_cox4/COXmodel_CyTOF_FACS.univariate.pdf ${DIR_FIGURE}/

# Figure 17 B

rsync -v -r -t ${RWD_MAIN}/PD-1_PatientClinicalData/ck_analysis_cox4/COXmodel_CyTOF_FACS.multivariate.pdf ${DIR_FIGURE}/




























