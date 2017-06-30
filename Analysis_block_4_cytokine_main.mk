### --------------------------------------------------------------------------
### Those parameters have to be defined
### --------------------------------------------------------------------------

R := R CMD BATCH --no-restore --no-save

RWD_MAIN := ../carsten_cytof/PD1_project
RCODE := .
METADATA := $(RWD_MAIN)/CK_metadata
PANELS := $(RWD_MAIN)/CK_panels

## arguments for the specific data set
data_dir := CK_2016-06-23_02_CD4_merging3_Tmem_merging1_CD69
file_metadata := $(METADATA)/metadata_23_02.xlsx
file_panel_cytokines := $(PANELS)/panel2CD4_23_cytokines.xlsx

data := 23CD4TmemCD69
cytokines := 02CD4

som_dim := 7
nmetaclusts := 49


### --------------------------------------------------------------------------
### Those parameters do not need to be defined
### --------------------------------------------------------------------------

RWD := $(RWD_MAIN)/$(data_dir)
ROUT := $(RWD)/Rout

rand_seed_consensus := 1234
make_file := Analysis_block_4_cytokine_main.mk

### --------------------------------------------------------------------------
## Define the default rule (makefiles are usually written so that the first target is for compiling the entire program)
.PHONY: all
all: mkdir_rout data_normalization_goal data_plot_goal \
	cytokine_expression_calculate_goal cytokine_expression_goal \
	cytokine_bimatrix_calculate_goal cytokine_bimatrix_frequencies_calculate_goal \
	cytokine_bimatrix_frequencies_overall_goal \
	cytokine_bimatrix_frequencies_clustering_goal \
	cytokine_bimatrix_frequencies_clustering_calculate_goal cytokine_bimatrix_frequencies_clustering_analysis_goal \
	cytokine_bimatrix_frequencies_clustering_plot_significant_goal



### Make sure no intermediate files are deleted
.SECONDARY:

.PHONY: mkdir_rout
mkdir_rout:
	mkdir -p $(ROUT)


### --------------------------------------------------------------------------
### Data normalization
### --------------------------------------------------------------------------

### Normalize (arcsinh) cytokines

.PHONY: data_normalization_goal
data_normalization_goal: $(RWD)/010_data/$(data)_$(cytokines)_expr_raw.rds

$(RWD)/010_data/$(data)_$(cytokines)_expr_raw.rds: $(RCODE)/01_data_normalization.R $(file_metadata) $(file_panel_cytokines) $(RWD)/010_cleanfcs/*.fcs
	echo "\n>> $(make_file)\n>>> 01_data_normalization"
	$(R) "--args dir_fcs='$(RWD)/010_cleanfcs' path_metadata='$(file_metadata)' path_panel='$(file_panel_cytokines)' outdir='$(RWD)/010_data' prefix='$(data)_$(cytokines)_'" $(RCODE)/01_data_normalization.R $(ROUT)/01_data_normalization.Rout


.PHONY: data_plot_goal
data_plot_goal: $(RWD)/010_data/$(data)_$(cytokines)_distrosgrp_raw.pdf $(RWD)/010_data/$(data)_$(cytokines)_distrosgrp_norm.pdf

### For raw data
$(RWD)/010_data/$(data)_$(cytokines)_distrosgrp_raw.pdf: $(RCODE)/01_data_plot.R $(file_metadata) $(file_panel_cytokines) $(RWD)/010_data/$(data)_$(cytokines)_expr_raw.rds
	echo "\n>> $(make_file)\n>>> 01_data_plot"
	$(R) "--args path_data='$(RWD)/010_data/$(data)_$(cytokines)_expr_raw.rds' path_metadata='$(file_metadata)' path_panel='$(file_panel_cytokines)' outdir='$(RWD)/010_data' prefix='$(data)_$(cytokines)_' suffix='_raw'" $(RCODE)/01_data_plot.R $(ROUT)/01_data_plot.Rout



### For norm data
$(RWD)/010_data/$(data)_$(cytokines)_distrosgrp_norm.pdf: $(RCODE)/01_data_plot.R $(file_metadata) $(file_panel_cytokines) $(RWD)/010_data/$(data)_$(cytokines)_expr_norm.rds
	echo "\n>> $(make_file)\n>>> 01_data_plot"
	$(R) "--args path_data='$(RWD)/010_data/$(data)_$(cytokines)_expr_norm.rds' path_metadata='$(file_metadata)' path_panel='$(file_panel_cytokines)' outdir='$(RWD)/010_data' prefix='$(data)_$(cytokines)_' suffix='_norm'" $(RCODE)/01_data_plot.R $(ROUT)/01_data_plot.Rout


### --------------------------------------------------------------------------
### Cytokine expression analysis
### --------------------------------------------------------------------------

.PHONY: cytokine_expression_calculate_goal
cytokine_expression_calculate_goal: $(RWD)/090_cytokine_expression/$(data)_$(cytokines)_all_expr.xls

$(RWD)/090_cytokine_expression/$(data)_$(cytokines)_all_expr.xls: $(RCODE)/05_cytokine_expression_calculate.R $(RWD)/010_data/$(data)_$(cytokines)_expr_raw.rds $(file_panel_cytokines)
	echo "\n>> $(make_file)\n>>> 05_cytokine_expression_calculate"
	$(R) "--args prefix='$(data)_$(cytokines)_' outdir='$(RWD)/090_cytokine_expression' path_data='$(RWD)/010_data/$(data)_$(cytokines)_expr_raw.rds' path_panel_cytokines='$(file_panel_cytokines)'" $(RCODE)/05_cytokine_expression_calculate.R $(ROUT)/05_cytokine_expression_calculate.Rout


analysis_type := all

.PHONY: cytokine_expression_goal
cytokine_expression_goal: $(foreach i,$(analysis_type),$(RWD)/090_cytokine_expression/$(data)_$(cytokines)_$(i)_expr_plot_both2.pdf) \
	$(foreach i,$(analysis_type),$(RWD)/090_cytokine_expression/$(data)_$(cytokines)_$(i)_expr_lmer_interglht_pheatmap3pvs_top10.pdf)

define 04_expression_plot_rule
$(RWD)/090_cytokine_expression/$(data)_$(cytokines)_$(1)_expr_plot_both2.pdf: $(RCODE)/04_expression_plot.R $(RCODE)/00_plot_expression.R $(RWD)/090_cytokine_expression/$(data)_$(cytokines)_$(1)_expr.xls
	echo "\n>> $(make_file)\n>>> 04_expression_plot"
	$(R) "--args prefix='$(data)_$(cytokines)_$(1)_' outdir='$(RWD)/090_cytokine_expression' path_expression='$(RWD)/090_cytokine_expression/$(data)_$(cytokines)_$(1)_expr.xls' path_metadata='$(file_metadata)' path_fun_plot_expression='$(RCODE)/00_plot_expression.R'" $(RCODE)/04_expression_plot.R $(ROUT)/04_expression_plot.Rout
endef
$(foreach i,$(analysis_type),$(eval $(call 04_expression_plot_rule,$(i))))


define 04_expression_analysis_rule
$(RWD)/090_cytokine_expression/$(data)_$(cytokines)_$(1)_expr_lmer_interglht_pheatmap3pvs_top10.pdf: $(RCODE)/04_expression_analysis.R $(RCODE)/00_models.R $(RCODE)/00_formulas_1dataset_3responses_both.R $(RCODE)/00_plot_heatmaps_for_sign_expr.R $(RWD)/090_cytokine_expression/$(data)_$(cytokines)_$(1)_expr.xls $(wildcard $(RWD)/010_helpfiles/$(data)_$(cytokines)_$(1)_marker_exclusion.txt)
	echo "\n>> $(make_file)\n>>> 04_expression_analysis"
	$(R) "--args prefix='$(data)_$(cytokines)_$(1)_' outdir='$(RWD)/090_cytokine_expression' path_expression='$(RWD)/090_cytokine_expression/$(data)_$(cytokines)_$(1)_expr.xls' path_metadata='$(file_metadata)' path_fun_models='$(RCODE)/00_models.R' path_fun_formulas='$(RCODE)/00_formulas_1dataset_3responses_both.R' \
	path_fun_plot_heatmaps='$(RCODE)/00_plot_heatmaps_for_sign_expr.R' path_marker_exclusion='$(RWD)/010_helpfiles/$(data)_$(cytokines)_$(1)_marker_exclusion.txt' FDR_cutoff='10'" $(RCODE)/04_expression_analysis.R $(ROUT)/04_expression_analysis.Rout
endef
$(foreach i,$(analysis_type),$(eval $(call 04_expression_analysis_rule,$(i))))


### --------------------------------------------------------------------------
### Cytokine bimatrix analysis
### --------------------------------------------------------------------------

### ---------------------------------------
### Generate the cytokine bimatrix
### ---------------------------------------

.PHONY: cytokine_bimatrix_calculate_goal
cytokine_bimatrix_calculate_goal: $(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_bimatrix.rds

$(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_bimatrix.rds: $(RCODE)/05_cytokine_bimatrix_calculate.R $(RWD)/010_data/$(data)_$(cytokines)_expr_raw.rds $(file_panel_cytokines)
	echo "\n>> $(make_file)\n>>> 05_cytokine_bimatrix_calculate"
	$(R) "--args prefix='$(data)_$(cytokines)_' outdir='$(RWD)/090_cytokine_bimatrix' path_data='$(RWD)/010_data/$(data)_$(cytokines)_expr_raw.rds' path_panel_cytokines='$(file_panel_cytokines)'" $(RCODE)/05_cytokine_bimatrix_calculate.R $(ROUT)/05_cytokine_bimatrix_calculate.Rout


### ---------------------------------------
### Calculate the frequencies of bimatrix
### ---------------------------------------

.PHONY: cytokine_bimatrix_frequencies_calculate_goal
cytokine_bimatrix_frequencies_calculate_goal: $(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_frequencies.xls

$(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_frequencies.xls: $(RCODE)/05_cytokine_bimatrix_frequencies_calculate.R $(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_bimatrix.rds
	echo "\n>> $(make_file)\n>>> 05_cytokine_bimatrix_frequencies_calculate"
	$(R) "--args prefix='$(data)_$(cytokines)_' outdir='$(RWD)/090_cytokine_bimatrix' path_bimatrix='$(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_bimatrix.rds'" $(RCODE)/05_cytokine_bimatrix_frequencies_calculate.R $(ROUT)/05_cytokine_bimatrix_frequencies_calculate.Rout


### ---------------------------------------
### Cytokine bimatrix frequencies overall
### ---------------------------------------

.PHONY: cytokine_bimatrix_frequencies_overall_goal
cytokine_bimatrix_frequencies_overall_goal: $(RWD)/090_cytokine_bimatrix_frequencies_overall/$(data)_$(cytokines)_frequencies_plot_both2.pdf  $(RWD)/090_cytokine_bimatrix_frequencies_overall/$(data)_$(cytokines)_frequencies_glmer_binomial_interglht_pheatmap3pvs_top10.pdf


$(RWD)/090_cytokine_bimatrix_frequencies_overall/$(data)_$(cytokines)_frequencies_plot_both2.pdf: $(RCODE)/04_frequencies_plot.R $(RCODE)/00_plot_frequencies.R $(file_metadata) $(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_frequencies.xls
	echo "\n>> $(make_file)\n>>> 04_frequencies_plot"
	$(R) "--args prefix='$(data)_$(cytokines)_' outdir='$(RWD)/090_cytokine_bimatrix_frequencies_overall' path_metadata='$(file_metadata)' path_frequencies='$(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_frequencies.xls' path_fun_plot_frequencies='$(RCODE)/00_plot_frequencies.R'" $(RCODE)/04_frequencies_plot.R $(ROUT)/04_frequencies_plot.Rout


$(RWD)/090_cytokine_bimatrix_frequencies_overall/$(data)_$(cytokines)_frequencies_glmer_binomial_interglht_pheatmap3pvs_top10.pdf: $(RCODE)/05_cytokine_bimatrix_overall_frequencies_analysis.R $(RCODE)/00_models.R $(RCODE)/00_formulas_1dataset_3responses_both.R $(RCODE)/00_plot_heatmaps_for_sign_freqs.R $(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_counts.xls $(file_metadata)
	echo "\n>> $(make_file)\n>>> 05_cytokine_bimatrix_overall_frequencies_analysis"
	$(R) "--args prefix='$(data)_$(cytokines)_' outdir='$(RWD)/090_cytokine_bimatrix_frequencies_overall' path_metadata='$(file_metadata)' path_counts='$(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_counts.xls' path_fun_models='$(RCODE)/00_models.R' path_fun_formulas='$(RCODE)/00_formulas_1dataset_3responses_both.R' path_fun_plot_heatmaps='$(RCODE)/00_plot_heatmaps_for_sign_freqs.R' FDR_cutoff='10'" $(RCODE)/05_cytokine_bimatrix_overall_frequencies_analysis.R $(ROUT)/05_cytokine_bimatrix_overall_frequencies_analysis.Rout


### ---------------------------------------
### Cytokine bimatrix frequencies clustering
### ---------------------------------------


.PHONY: cytokine_bimatrix_frequencies_clustering_goal
cytokine_bimatrix_frequencies_clustering_goal: $(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_clustering.xls $(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_pheatmap_all_all_row_clust_raw.pdf


# FlowSOM clustering of the bimatrix

$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_clustering.xls: $(RCODE)/02_flowsom_bimatrix.R $(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_bimatrix.rds $(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_clustering_observables.xls
	echo "\n>> $(make_file)\n>>> 02_flowsom_bimatrix"
	$(R) "--args prefix='$(data)_$(cytokines)_cl$(nmetaclusts)_' outdir='$(RWD)/090_cytokine_bimatrix_frequencies_clustering' path_data='$(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_bimatrix.rds' path_clustering_observables='$(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_clustering_observables.xls' nmetaclusts=$(nmetaclusts) som_dim=$(som_dim)" $(RCODE)/02_flowsom_bimatrix.R $(ROUT)/02_flowsom_bimatrix.Rout

# Plot a heatmap of the bimatrix clustering

$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_pheatmap_all_all_row_clust_raw.pdf: $(RCODE)/02_heatmaps.R $(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_bimatrix.rds $(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_clustering.xls
	echo "\n>> $(make_file)\n>>> 02_heatmaps"
	$(R) "--args prefix='$(data)_$(cytokines)_cl$(nmetaclusts)_' outdir='$(RWD)/090_cytokine_bimatrix_frequencies_clustering' path_data='$(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_bimatrix.rds' path_data_norm=NULL \
	path_clustering_observables='$(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_clustering_observables.xls' path_clustering='$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_clustering.xls'  path_clustering_labels='$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_clustering_labels.xls' \
	path_marker_selection=NULL path_cluster_merging=NULL aggregate_fun='mean' scale=FALSE" $(RCODE)/02_heatmaps.R $(ROUT)/02_heatmaps.Rout



# Calculate frequancies of the bimatrix clusters

.PHONY: cytokine_bimatrix_frequencies_clustering_calculate_goal
cytokine_bimatrix_frequencies_clustering_calculate_goal: $(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_frequencies.xls


$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_frequencies.xls: $(RCODE)/04_frequencies_calculate.R $(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_clustering.xls
	echo "\n>> $(make_file)\n>>> 04_frequencies_calculate"
	$(R) "--args prefix='$(data)_$(cytokines)_cl$(nmetaclusts)_' outdir='$(RWD)/090_cytokine_bimatrix_frequencies_clustering' path_clustering='$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_clustering.xls' path_clustering_labels='$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_clustering_labels.xls'" $(RCODE)/04_frequencies_calculate.R $(ROUT)/04_frequencies_calculate.Rout



# Differential analysis of the frequancies of the bimatrix clusters

model := glmer_binomial_interglht

.PHONY: cytokine_bimatrix_frequencies_clustering_analysis_goal
cytokine_bimatrix_frequencies_clustering_analysis_goal: $(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_frequencies_$(model)_pheatmap3pvs_top10.pdf


$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_frequencies_$(model)_pheatmap3pvs_top10.pdf: $(RCODE)/04_frequencies_analysis.R $(RCODE)/00_models.R $(RCODE)/00_formulas_1dataset_3responses_both.R $(RCODE)/00_plot_heatmaps_for_sign_freqs.R $(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_counts.xls $(file_metadata)
	echo "\n>> $(make_file)\n>>> 04_frequencies_analysis"
	$(R) "--args prefix='$(data)_$(cytokines)_cl$(nmetaclusts)_' outdir='$(RWD)/090_cytokine_bimatrix_frequencies_clustering' path_metadata='$(file_metadata)' path_counts='$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_counts.xls' path_fun_models='$(RCODE)/00_models.R' path_fun_formulas='$(RCODE)/00_formulas_1dataset_3responses_both.R' path_fun_plot_heatmaps='$(RCODE)/00_plot_heatmaps_for_sign_freqs.R' FDR_cutoff='10'" $(RCODE)/04_frequencies_analysis.R $(ROUT)/04_frequencies_analysis.Rout



# Plot frequencies of the significant bimatrix clusters

.PHONY: cytokine_bimatrix_frequencies_clustering_plot_significant_goal
cytokine_bimatrix_frequencies_clustering_plot_significant_goal:  $(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_top10_$(model)_NRvsR_base_frequencies_plot_both2.pdf $(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_top10_$(model)_NRvsR_base_pheatmap.pdf


$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_top10_$(model)_NRvsR_base_frequencies_plot_both2.pdf: $(RCODE)/04_frequencies_plot_significant.R $(RCODE)/00_plot_frequencies.R $(file_metadata) $(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_frequencies.xls $(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_frequencies_pvs_$(model)_top10.xls
	echo "\n>> $(make_file)\n>>> 04_frequencies_plot_significant"
	$(R) "--args prefix='$(data)_$(cytokines)_cl$(nmetaclusts)_top10_$(model)_' outdir='$(RWD)/090_cytokine_bimatrix_frequencies_clustering' path_metadata='$(file_metadata)' path_frequencies='$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_frequencies.xls' path_pvs='$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_frequencies_pvs_$(model)_top10.xls' path_fun_plot_frequencies='$(RCODE)/00_plot_frequencies.R' FDR_cutoff='10'" $(RCODE)/04_frequencies_plot_significant.R $(ROUT)/04_frequencies_plot_significant.Rout


# Plot a heatmap of the significant bimatrix clusters

$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_top10_$(model)_NRvsR_base_pheatmap.pdf: $(RCODE)/02_heatmaps_significant.R $(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_bimatrix.rds $(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_clustering.xls
	echo "\n>> $(make_file)\n>>> 02_heatmaps_significant"
	$(R) "--args prefix='$(data)_$(cytokines)_cl$(nmetaclusts)_top10_$(model)_' outdir='$(RWD)/090_cytokine_bimatrix_frequencies_clustering' path_data='$(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_bimatrix.rds' \
	path_clustering_observables='$(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_clustering_observables.xls' path_clustering='$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_clustering.xls'  path_clustering_labels='$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_clustering_labels.xls' \
	aggregate_fun='mean' path_pvs='$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_frequencies_pvs_$(model)_top10.xls' FDR_cutoff='10'" $(RCODE)/02_heatmaps_significant.R $(ROUT)/02_heatmaps_significant.Rout






### Cytokine bimatrix frequencies top selection









































#
