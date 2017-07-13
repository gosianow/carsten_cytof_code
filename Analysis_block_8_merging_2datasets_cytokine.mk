### --------------------------------------------------------------------------
### Those parameters have to be defined
### --------------------------------------------------------------------------

R := R CMD BATCH --no-restore --no-save

RWD_MAIN := ../carsten_cytof/PD1_project
RCODE := .
METADATA := $(RWD_MAIN)/CK_metadata
PANELS := $(RWD_MAIN)/CK_panels

## arguments for the specific data set
data_dir_merged := CK_2016-06-merged_23_29/02_CD4

data_dir1 := CK_2016-06-23_02_CD4_merging2_Tmem_merging2_CD69
data_dir2 := CK_2016-06-29_02_CD4_merging_Tmem_merging3_CD69
file_metadata1 := $(METADATA)/metadata_23_02.xlsx
file_metadata2 := $(METADATA)/metadata_29_02.xlsx

data1 := 23CD4TmemCD69
data2 := 29CD4TmemCD69

cytokines := 02CD4

som_dim := 7
nmetaclusts := 49

data_dir_back1 := CK_2016-06-23_02_CD4_merging2
data_dir_back2 := CK_2016-06-29_02_CD4_merging
data_back1 := 23CD4
data_back2 := 29CD4
panel_back1 := 02CD4
panel_back2 := 02CD4
pca_back1 := pca1
pca_back2 := pca1
merging_back1 := merging2
merging_back2 := merging3

### --------------------------------------------------------------------------
### Those parameters do not need to be defined
### --------------------------------------------------------------------------

RWD1 := $(RWD_MAIN)/$(data_dir1)
RWD2 := $(RWD_MAIN)/$(data_dir2)
RWD_MERGED := $(RWD_MAIN)/$(data_dir_merged)
ROUT := $(RWD_MERGED)/Rout

RWD_BACK1 := $(RWD_MAIN)/$(data_dir_back1)
RWD_BACK2 := $(RWD_MAIN)/$(data_dir_back2)

data_back := $(data_back1) $(data_back2)

rand_seed_consensus := 1234
make_file := Analysis_block_8_merging_2datasets_cytokine.mk

### --------------------------------------------------------------------------
## Define the default rule (makefiles are usually written so that the first target is for compiling the entire program)
.PHONY: all
all: mkdir_rout cytokine_expression_goal \
	cytokine_bimatrix_frequencies_calculate_goal \
	cytokine_bimatrix_frequencies_overall_goal cytokine_bimatrix_frequencies_overall_cytokine_profiles_get_clustering cytokine_bimatrix_frequencies_overall_cytokine_profiles_heatmaps1 cytokine_bimatrix_frequencies_overall_cytokine_profiles_heatmaps2 \
	cytokine_bimatrix_frequencies_clustering_goal cytokine_bimatrix_frequencies_clustering_calculate_goal cytokine_bimatrix_frequencies_clustering_analysis_goal cytokine_bimatrix_frequencies_clustering_plot_significant_goal \
	cytokine_bimatrix_frequencies_clustering_cytokine_profiles_get_clustering cytokine_bimatrix_frequencies_clustering_cytokine_profiles_heatmaps1 cytokine_bimatrix_frequencies_clustering_cytokine_profiles_heatmaps2


### Make sure no intermediate files are deleted
.SECONDARY:

.PHONY: mkdir_rout
mkdir_rout:
	mkdir -p $(ROUT)


### --------------------------------------------------------------------------
### Cytokine expression analysis
### --------------------------------------------------------------------------


analysis_type := all


.PHONY: cytokine_expression_goal
cytokine_expression_goal: $(foreach i,$(analysis_type),$(RWD_MERGED)/090_cytokine_expression/3responses_both/$(data1)_$(data2)_$(cytokines)_$(i)_expr_plot_both2.pdf) \
	$(foreach i,$(analysis_type),$(RWD_MERGED)/090_cytokine_expression/3responses_both/$(data1)_$(data2)_$(cytokines)_$(i)_expr_lmer_interglht_pheatmap3pvs_top10.pdf)


define 04_expression_plot_rule
$(RWD_MERGED)/090_cytokine_expression/3responses_both/$(data1)_$(data2)_$(cytokines)_$(1)_expr_plot_both2.pdf: $(RCODE)/04_expression_plot.R $(RCODE)/00_plot_expression.R $(RWD1)/090_cytokine_expression/$(data1)_$(cytokines)_$(1)_expr.xls $(RWD2)/090_cytokine_expression/$(data2)_$(cytokines)_$(1)_expr.xls
	echo "\n>> $(make_file)\n>>> 04_expression_plot"
	$(R) "--args prefix='$(data1)_$(data2)_$(cytokines)_$(1)_' outdir='$(RWD_MERGED)/090_cytokine_expression/3responses_both' path_expression=c('$(RWD1)/090_cytokine_expression/$(data1)_$(cytokines)_$(1)_expr.xls','$(RWD2)/090_cytokine_expression/$(data2)_$(cytokines)_$(1)_expr.xls') path_metadata=c('$(file_metadata1)','$(file_metadata2)') path_fun_plot_expression='$(RCODE)/00_plot_expression.R'" $(RCODE)/04_expression_plot.R $(ROUT)/04_expression_plot.Rout
endef
$(foreach i,$(analysis_type),$(eval $(call 04_expression_plot_rule,$(i))))



define 04_expression_analysis_rule
$(RWD_MERGED)/090_cytokine_expression/3responses_both/$(data1)_$(data2)_$(cytokines)_$(1)_expr_lmer_interglht_pheatmap3pvs_top10.pdf: $(RCODE)/04_expression_analysis.R $(RCODE)/00_models.R $(RCODE)/00_formulas_2datasets_3responses_both.R $(RCODE)/00_plot_heatmaps_for_sign_expr.R $(RWD1)/090_cytokine_expression/$(data1)_$(cytokines)_$(1)_expr.xls $(RWD2)/090_cytokine_expression/$(data2)_$(cytokines)_$(1)_expr.xls
	echo "\n>> $(make_file)\n>>> 04_expression_analysis"
	$(R) "--args prefix='$(data1)_$(data2)_$(cytokines)_$(1)_' outdir='$(RWD_MERGED)/090_cytokine_expression/3responses_both' path_expression=c('$(RWD1)/090_cytokine_expression/$(data1)_$(cytokines)_$(1)_expr.xls','$(RWD2)/090_cytokine_expression/$(data2)_$(cytokines)_$(1)_expr.xls') path_metadata=c('$(file_metadata1)','$(file_metadata2)') path_fun_models='$(RCODE)/00_models.R' path_fun_formulas='$(RCODE)/00_formulas_2datasets_3responses_both.R' \
	path_fun_plot_heatmaps='$(RCODE)/00_plot_heatmaps_for_sign_expr.R' path_marker_exclusion=NULL FDR_cutoff='10'" $(RCODE)/04_expression_analysis.R $(ROUT)/04_expression_analysis.Rout
endef
$(foreach i,$(analysis_type),$(eval $(call 04_expression_analysis_rule,$(i))))


### --------------------------------------------------------------------------
### Cytokine bimatrix analysis
### --------------------------------------------------------------------------

### ---------------------------------------
### Calculate the frequencies of bimatrix
### ---------------------------------------

.PHONY: cytokine_bimatrix_frequencies_calculate_goal
cytokine_bimatrix_frequencies_calculate_goal: $(RWD_MERGED)/090_cytokine_bimatrix/$(data1)_$(data2)_$(cytokines)_frequencies.xls

$(RWD_MERGED)/090_cytokine_bimatrix/$(data1)_$(data2)_$(cytokines)_frequencies.xls: $(RCODE)/05_cytokine_bimatrix_frequencies_calculate.R $(RWD1)/090_cytokine_bimatrix/$(data1)_$(cytokines)_bimatrix.rds $(RWD2)/090_cytokine_bimatrix/$(data2)_$(cytokines)_bimatrix.rds
	echo "\n>> $(make_file)\n>>> 05_cytokine_bimatrix_frequencies_calculate"
	$(R) "--args prefix='$(data1)_$(data2)_$(cytokines)_' outdir='$(RWD_MERGED)/090_cytokine_bimatrix' path_bimatrix=c('$(RWD1)/090_cytokine_bimatrix/$(data1)_$(cytokines)_bimatrix.rds','$(RWD2)/090_cytokine_bimatrix/$(data2)_$(cytokines)_bimatrix.rds')" $(RCODE)/05_cytokine_bimatrix_frequencies_calculate.R $(ROUT)/05_cytokine_bimatrix_frequencies_calculate.Rout


### ---------------------------------------
### Cytokine bimatrix frequencies overall
### ---------------------------------------

.PHONY: cytokine_bimatrix_frequencies_overall_goal
cytokine_bimatrix_frequencies_overall_goal: $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/3responses_both/$(data1)_$(data2)_$(cytokines)_frequencies_plot_both2.pdf  $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/3responses_both/$(data1)_$(data2)_$(cytokines)_frequencies_glmer_binomial_interglht_pheatmap3pvs_top10.pdf


$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/3responses_both/$(data1)_$(data2)_$(cytokines)_frequencies_plot_both2.pdf: $(RCODE)/04_frequencies_plot.R $(RCODE)/00_plot_frequencies.R $(file_metadata1) $(file_metadata2) $(RWD_MERGED)/090_cytokine_bimatrix/$(data1)_$(data2)_$(cytokines)_frequencies.xls
	echo "\n>> $(make_file)\n>>> 04_frequencies_plot"
	$(R) "--args prefix='$(data1)_$(data2)_$(cytokines)_' outdir='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/3responses_both' path_metadata=c('$(file_metadata1)','$(file_metadata2)') path_frequencies='$(RWD_MERGED)/090_cytokine_bimatrix/$(data1)_$(data2)_$(cytokines)_frequencies.xls' path_fun_plot_frequencies='$(RCODE)/00_plot_frequencies.R'" $(RCODE)/04_frequencies_plot.R $(ROUT)/04_frequencies_plot.Rout


$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/3responses_both/$(data1)_$(data2)_$(cytokines)_frequencies_glmer_binomial_interglht_pheatmap3pvs_top10.pdf: $(RCODE)/05_cytokine_bimatrix_overall_frequencies_analysis.R $(RCODE)/00_models.R $(RCODE)/00_formulas_2datasets_3responses_both.R $(RCODE)/00_plot_heatmaps_for_sign_freqs.R $(RWD_MERGED)/090_cytokine_bimatrix/$(data1)_$(data2)_$(cytokines)_counts.xls $(file_metadata)
	echo "\n>> $(make_file)\n>>> 05_cytokine_bimatrix_overall_frequencies_analysis"
	$(R) "--args prefix='$(data1)_$(data2)_$(cytokines)_' outdir='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/3responses_both' path_metadata=c('$(file_metadata1)','$(file_metadata2)') path_counts='$(RWD_MERGED)/090_cytokine_bimatrix/$(data1)_$(data2)_$(cytokines)_counts.xls' path_fun_models='$(RCODE)/00_models.R' path_fun_formulas='$(RCODE)/00_formulas_2datasets_3responses_both.R' path_fun_plot_heatmaps='$(RCODE)/00_plot_heatmaps_for_sign_freqs.R' FDR_cutoff='10'" $(RCODE)/05_cytokine_bimatrix_overall_frequencies_analysis.R $(ROUT)/05_cytokine_bimatrix_overall_frequencies_analysis.Rout


### ----------------------------
### Cytokine profiles for bimatrix frequencies overall
### ----------------------------

# Get the Tmem clustering for each of the positive cytokines
.PHONY: cytokine_bimatrix_frequencies_overall_cytokine_profiles_get_clustering
cytokine_bimatrix_frequencies_overall_cytokine_profiles_get_clustering: $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_clustering_data$(data_back1)_clPD1.txt


$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_clustering_data$(data_back1)_clPD1.txt: $(RCODE)/06_cytokine_profiles_get_clustering_overall.R $(RWD_MERGED)/090_cytokine_bimatrix/$(data1)_$(data2)_$(cytokines)_bimatrix.rds
	echo "\n>> $(make_file)\n>>> 06_cytokine_profiles_get_clustering_overall"
	$(R) "--args prefix='$(data1)_$(data2)_$(cytokines)_' outdir='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/cytokine_profiles' \
	path_bimatrix='$(RWD_MERGED)/090_cytokine_bimatrix/$(data1)_$(data2)_$(cytokines)_bimatrix.rds' \
	path_clustering=c('$(RWD_BACK1)/030_heatmaps/$(data_back1)_$(panel_back1)_$(pca_back1)_$(merging_back1)_clustering.xls','$(RWD_BACK2)/030_heatmaps/$(data_back2)_$(panel_back2)_$(pca_back2)_$(merging_back2)_clustering.xls') \
	path_cells2keep=c('$(RWD1)/010_cleanfcs/cells2keep.txt','$(RWD2)/010_cleanfcs/cells2keep.txt') data=c('$(data_back1)','$(data_back2)')" $(RCODE)/06_cytokine_profiles_get_clustering_overall.R $(ROUT)/06_cytokine_profiles_get_clustering_overall.Rout


### Plot heatmaps with cytokine profiles

# For data_back1 RWD_BACK1
generated_clustering := $(wildcard $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_clustering_data$(data_back1)_cl*.txt)

generated_clustering := $(foreach i,$(generated_clustering),$(subst $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_clustering_data$(data_back1)_cl,,$(i)))

generated_clustering := $(foreach i,$(generated_clustering),$(subst .txt,,$(i)))


.PHONY: cytokine_bimatrix_frequencies_overall_cytokine_profiles_heatmaps1
cytokine_bimatrix_frequencies_overall_cytokine_profiles_heatmaps1: $(foreach i,$(generated_clustering),$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_clustering_data$(data_back1)_cl$(i)_pheatmap_all_all_no_clust_raw.pdf)


define cytokine_bimatrix_frequencies_overall_cytokine_profiles_heatmaps_rule
$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_clustering_data$(data_back1)_cl$(1)_pheatmap_all_all_no_clust_raw.pdf: $(RCODE)/02_heatmaps.R $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_clustering_data$(data_back1)_cl$(1).txt $(RWD_BACK1)/010_data/$(data_back1)_$(panel_back1)_expr_raw.rds
	echo "\n>> $(make_file)\n>>> 02_heatmaps"
	$(R) "--args prefix='$(data1)_$(data2)_$(cytokines)_clustering_data$(data_back1)_cl$(1)_' outdir='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/cytokine_profiles' path_data='$(RWD_BACK1)/010_data/$(data_back1)_$(panel_back1)_expr_raw.rds' path_data_norm='$(RWD_BACK1)/010_data/$(data_back1)_$(panel_back1)_expr_norm.rds' \
	path_clustering_observables='$(RWD_BACK1)/030_heatmaps/$(data_back1)_$(panel_back1)_$(pca_back1)_clustering_observables.xls' path_clustering='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_clustering_data$(data_back1)_cl$(1).txt'  path_clustering_labels='$(RWD_BACK1)/030_heatmaps/$(data_back1)_$(panel_back1)_$(pca_back1)_$(merging_back1)_clustering_labels.xls' \
	path_marker_selection='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/marker_selection.txt' path_cluster_merging=NULL aggregate_fun='median' scale=FALSE" $(RCODE)/02_heatmaps.R $(ROUT)/02_heatmaps.Rout
endef
$(foreach i,$(generated_clustering),$(eval $(call cytokine_bimatrix_frequencies_overall_cytokine_profiles_heatmaps_rule,$(i))))



# For data_back2 RWD_BACK2
generated_clustering := $(wildcard $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_clustering_data$(data_back2)_cl*.txt)

generated_clustering := $(foreach i,$(generated_clustering),$(subst $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_clustering_data$(data_back2)_cl,,$(i)))

generated_clustering := $(foreach i,$(generated_clustering),$(subst .txt,,$(i)))


.PHONY: cytokine_bimatrix_frequencies_overall_cytokine_profiles_heatmaps2
cytokine_bimatrix_frequencies_overall_cytokine_profiles_heatmaps2: $(foreach i,$(generated_clustering),$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_clustering_data$(data_back2)_cl$(i)_pheatmap_all_all_no_clust_raw.pdf)


define cytokine_bimatrix_frequencies_overall_cytokine_profiles_heatmaps_rule
$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_clustering_data$(data_back2)_cl$(1)_pheatmap_all_all_no_clust_raw.pdf: $(RCODE)/02_heatmaps.R $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_clustering_data$(data_back2)_cl$(1).txt $(RWD_BACK2)/010_data/$(data_back2)_$(panel_back2)_expr_raw.rds
	echo "\n>> $(make_file)\n>>> 02_heatmaps"
	$(R) "--args prefix='$(data1)_$(data2)_$(cytokines)_clustering_data$(data_back2)_cl$(1)_' outdir='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/cytokine_profiles' path_data='$(RWD_BACK2)/010_data/$(data_back2)_$(panel_back2)_expr_raw.rds' path_data_norm='$(RWD_BACK2)/010_data/$(data_back2)_$(panel_back2)_expr_norm.rds' \
	path_clustering_observables='$(RWD_BACK2)/030_heatmaps/$(data_back2)_$(panel_back2)_$(pca_back2)_clustering_observables.xls' path_clustering='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_clustering_data$(data_back2)_cl$(1).txt'  path_clustering_labels='$(RWD_BACK2)/030_heatmaps/$(data_back2)_$(panel_back2)_$(pca_back2)_$(merging_back2)_clustering_labels.xls' \
	path_marker_selection='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_overall/marker_selection.txt' path_cluster_merging=NULL aggregate_fun='median' scale=FALSE" $(RCODE)/02_heatmaps.R $(ROUT)/02_heatmaps.Rout
endef
$(foreach i,$(generated_clustering),$(eval $(call cytokine_bimatrix_frequencies_overall_cytokine_profiles_heatmaps_rule,$(i))))







### ------------------------------------------------------------------------------
### Cytokine bimatrix frequencies clustering
### ------------------------------------------------------------------------------


.PHONY: cytokine_bimatrix_frequencies_clustering_goal
cytokine_bimatrix_frequencies_clustering_goal: $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering.xls $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_pheatmap_all_all_row_clust_raw.pdf


# FlowSOM clustering of the bimatrix

$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering.xls: $(RCODE)/02_flowsom_bimatrix.R $(RWD_MERGED)/090_cytokine_bimatrix/$(data1)_$(data2)_$(cytokines)_bimatrix.rds $(RWD_MERGED)/090_cytokine_bimatrix/$(data1)_$(data2)_$(cytokines)_clustering_observables.xls
	echo "\n>> $(make_file)\n>>> 02_flowsom_bimatrix"
	$(R) "--args prefix='$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_' outdir='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering' path_data='$(RWD_MERGED)/090_cytokine_bimatrix/$(data1)_$(data2)_$(cytokines)_bimatrix.rds' path_clustering_observables='$(RWD_MERGED)/090_cytokine_bimatrix/$(data1)_$(data2)_$(cytokines)_clustering_observables.xls' nmetaclusts=$(nmetaclusts) som_dim=$(som_dim)" $(RCODE)/02_flowsom_bimatrix.R $(ROUT)/02_flowsom_bimatrix.Rout

# Plot a heatmap of the bimatrix clustering

$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_pheatmap_all_all_row_clust_raw.pdf: $(RCODE)/02_heatmaps.R $(RWD_MERGED)/090_cytokine_bimatrix/$(data1)_$(data2)_$(cytokines)_bimatrix.rds $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering.xls
	echo "\n>> $(make_file)\n>>> 02_heatmaps"
	$(R) "--args prefix='$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_' outdir='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering' path_data='$(RWD_MERGED)/090_cytokine_bimatrix/$(data1)_$(data2)_$(cytokines)_bimatrix.rds' path_data_norm=NULL \
	path_clustering_observables='$(RWD_MERGED)/090_cytokine_bimatrix/$(data1)_$(data2)_$(cytokines)_clustering_observables.xls' path_clustering='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering.xls'  path_clustering_labels='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering_labels.xls' \
	path_marker_selection=NULL path_cluster_merging=NULL aggregate_fun='mean' scale=FALSE" $(RCODE)/02_heatmaps.R $(ROUT)/02_heatmaps.Rout


# Calculate frequancies of the bimatrix clusters


.PHONY: cytokine_bimatrix_frequencies_clustering_calculate_goal
cytokine_bimatrix_frequencies_clustering_calculate_goal: $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/3responses_both/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_frequencies.xls


$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/3responses_both/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_frequencies.xls: $(RCODE)/04_frequencies_calculate.R $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering.xls
	echo "\n>> $(make_file)\n>>> 04_frequencies_calculate"
	$(R) "--args prefix='$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_' outdir='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/3responses_both' path_clustering='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering.xls' path_clustering_labels='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering_labels.xls'" $(RCODE)/04_frequencies_calculate.R $(ROUT)/04_frequencies_calculate.Rout



# Differential analysis of the frequancies of the bimatrix clusters

model := glmer_binomial_interglht

.PHONY: cytokine_bimatrix_frequencies_clustering_analysis_goal
cytokine_bimatrix_frequencies_clustering_analysis_goal: $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/3responses_both/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_frequencies_$(model)_pheatmap3pvs_top10.pdf


$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/3responses_both/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_frequencies_$(model)_pheatmap3pvs_top10.pdf: $(RCODE)/04_frequencies_analysis.R $(RCODE)/00_models.R $(RCODE)/00_formulas_2datasets_3responses_both.R $(RCODE)/00_plot_heatmaps_for_sign_freqs.R $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/3responses_both/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_counts.xls $(file_metadata)
	echo "\n>> $(make_file)\n>>> 04_frequencies_analysis"
	$(R) "--args prefix='$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_' outdir='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/3responses_both' path_metadata=c('$(file_metadata1)','$(file_metadata2)') path_counts='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/3responses_both/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_counts.xls' path_fun_models='$(RCODE)/00_models.R' path_fun_formulas='$(RCODE)/00_formulas_2datasets_3responses_both.R' path_fun_plot_heatmaps='$(RCODE)/00_plot_heatmaps_for_sign_freqs.R' FDR_cutoff='10'" $(RCODE)/04_frequencies_analysis.R $(ROUT)/04_frequencies_analysis.Rout




# Plot the significant bimatrix clusters

.PHONY: cytokine_bimatrix_frequencies_clustering_plot_significant_goal
cytokine_bimatrix_frequencies_clustering_plot_significant_goal: $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/3responses_both/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_top10_$(model)_NRvsR_base_frequencies_plot_both2.pdf $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/3responses_both/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_top10_$(model)_NRvsR_base_pheatmap.pdf


$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/3responses_both/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_top10_$(model)_NRvsR_base_frequencies_plot_both2.pdf: $(RCODE)/04_frequencies_plot_significant.R $(RCODE)/00_plot_frequencies.R $(file_metadata) $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/3responses_both/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_frequencies.xls $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/3responses_both/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_frequencies_pvs_$(model)_top10.xls
	echo "\n>> $(make_file)\n>>> 04_frequencies_plot_significant"
	$(R) "--args prefix='$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_top10_$(model)_' outdir='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/3responses_both' path_metadata=c('$(file_metadata1)','$(file_metadata2)') path_frequencies='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/3responses_both/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_frequencies.xls' path_pvs='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/3responses_both/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_frequencies_pvs_$(model)_top10.xls' path_fun_plot_frequencies='$(RCODE)/00_plot_frequencies.R' FDR_cutoff='10'" $(RCODE)/04_frequencies_plot_significant.R $(ROUT)/04_frequencies_plot_significant.Rout


# Plot a heatmap of the significant bimatrix clusters

$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/3responses_both/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_top10_$(model)_NRvsR_base_pheatmap.pdf: $(RCODE)/02_heatmaps_significant.R $(RWD_MERGED)/090_cytokine_bimatrix/$(data1)_$(data2)_$(cytokines)_bimatrix.rds $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering.xls
	echo "\n>> $(make_file)\n>>> 02_heatmaps_significant"
	$(R) "--args prefix='$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_top10_$(model)_' outdir='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/3responses_both' path_data='$(RWD_MERGED)/090_cytokine_bimatrix/$(data1)_$(data2)_$(cytokines)_bimatrix.rds' \
	path_clustering_observables='$(RWD_MERGED)/090_cytokine_bimatrix/$(data1)_$(data2)_$(cytokines)_clustering_observables.xls' path_clustering='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering.xls'  path_clustering_labels='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering_labels.xls' \
	aggregate_fun='mean' path_pvs='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/3responses_both/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_frequencies_pvs_$(model)_top10.xls' FDR_cutoff='10'" $(RCODE)/02_heatmaps_significant.R $(ROUT)/02_heatmaps_significant.Rout


### ----------------------------
### Cytokine profiles for bimatrix frequencies clustering
### ----------------------------

# Get the Tmem clustering for each of the bimatrix clusters
.PHONY: cytokine_bimatrix_frequencies_clustering_cytokine_profiles_get_clustering
cytokine_bimatrix_frequencies_clustering_cytokine_profiles_get_clustering: $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering_data$(data_back1)_cl1.txt


$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering_data$(data_back1)_cl1.txt: $(RCODE)/06_cytokine_profiles_get_clustering.R $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering.xls
	echo "\n>> $(make_file)\n>>> 06_cytokine_profiles_get_clustering"
	$(R) "--args prefix='$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_' outdir='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/cytokine_profiles' \
	path_bimatrix_clustering='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering.xls' \
	path_clustering=c('$(RWD_BACK1)/030_heatmaps/$(data_back1)_$(panel_back1)_$(pca_back1)_$(merging_back1)_clustering.xls','$(RWD_BACK2)/030_heatmaps/$(data_back2)_$(panel_back2)_$(pca_back2)_$(merging_back2)_clustering.xls') \
	path_cells2keep=c('$(RWD1)/010_cleanfcs/cells2keep.txt','$(RWD2)/010_cleanfcs/cells2keep.txt') data=c('$(data_back1)','$(data_back2)')" $(RCODE)/06_cytokine_profiles_get_clustering.R $(ROUT)/06_cytokine_profiles_get_clustering.Rout


### Plot heatmaps with cytokine profiles

# For data_back1 RWD_BACK1
generated_clustering := $(wildcard $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering_data$(data_back1)_cl*.txt)

generated_clustering := $(foreach i,$(generated_clustering),$(subst $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering_data$(data_back1)_cl,,$(i)))

generated_clustering := $(foreach i,$(generated_clustering),$(subst .txt,,$(i)))


.PHONY: cytokine_bimatrix_frequencies_clustering_cytokine_profiles_heatmaps1
cytokine_bimatrix_frequencies_clustering_cytokine_profiles_heatmaps1: $(foreach i,$(generated_clustering),$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering_data$(data_back1)_cl$(i)_pheatmap_all_all_no_clust_raw.pdf)


define cytokine_bimatrix_frequencies_clustering_cytokine_profiles_heatmaps_rule
$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering_data$(data_back1)_cl$(1)_pheatmap_all_all_no_clust_raw.pdf: $(RCODE)/02_heatmaps.R $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering_data$(data_back1)_cl$(1).txt $(RWD_BACK1)/010_data/$(data_back1)_$(panel_back1)_expr_raw.rds
	echo "\n>> $(make_file)\n>>> 02_heatmaps"
	$(R) "--args prefix='$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering_data$(data_back1)_cl$(1)_' outdir='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/cytokine_profiles' path_data='$(RWD_BACK1)/010_data/$(data_back1)_$(panel_back1)_expr_raw.rds' path_data_norm='$(RWD_BACK1)/010_data/$(data_back1)_$(panel_back1)_expr_norm.rds' \
	path_clustering_observables='$(RWD_BACK1)/030_heatmaps/$(data_back1)_$(panel_back1)_$(pca_back1)_clustering_observables.xls' path_clustering='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering_data$(data_back1)_cl$(1).txt'  path_clustering_labels='$(RWD_BACK1)/030_heatmaps/$(data_back1)_$(panel_back1)_$(pca_back1)_$(merging_back1)_clustering_labels.xls' \
	path_marker_selection='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/marker_selection.txt' path_cluster_merging=NULL aggregate_fun='median' scale=FALSE" $(RCODE)/02_heatmaps.R $(ROUT)/02_heatmaps.Rout
endef
$(foreach i,$(generated_clustering),$(eval $(call cytokine_bimatrix_frequencies_clustering_cytokine_profiles_heatmaps_rule,$(i))))


# For data_back2 RWD_BACK2
generated_clustering := $(wildcard $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering_data$(data_back2)_cl*.txt)

generated_clustering := $(foreach i,$(generated_clustering),$(subst $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering_data$(data_back2)_cl,,$(i)))

generated_clustering := $(foreach i,$(generated_clustering),$(subst .txt,,$(i)))


.PHONY: cytokine_bimatrix_frequencies_clustering_cytokine_profiles_heatmaps2
cytokine_bimatrix_frequencies_clustering_cytokine_profiles_heatmaps2: $(foreach i,$(generated_clustering),$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering_data$(data_back2)_cl$(i)_pheatmap_all_all_no_clust_raw.pdf)


define cytokine_bimatrix_frequencies_clustering_cytokine_profiles_heatmaps_rule
$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering_data$(data_back2)_cl$(1)_pheatmap_all_all_no_clust_raw.pdf: $(RCODE)/02_heatmaps.R $(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering_data$(data_back2)_cl$(1).txt $(RWD_BACK2)/010_data/$(data_back2)_$(panel_back2)_expr_raw.rds
	echo "\n>> $(make_file)\n>>> 02_heatmaps"
	$(R) "--args prefix='$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering_data$(data_back2)_cl$(1)_' outdir='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/cytokine_profiles' path_data='$(RWD_BACK2)/010_data/$(data_back2)_$(panel_back2)_expr_raw.rds' path_data_norm='$(RWD_BACK2)/010_data/$(data_back2)_$(panel_back2)_expr_norm.rds' \
	path_clustering_observables='$(RWD_BACK2)/030_heatmaps/$(data_back2)_$(panel_back2)_$(pca_back2)_clustering_observables.xls' path_clustering='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/cytokine_profiles/$(data1)_$(data2)_$(cytokines)_cl$(nmetaclusts)_clustering_data$(data_back2)_cl$(1).txt'  path_clustering_labels='$(RWD_BACK2)/030_heatmaps/$(data_back2)_$(panel_back2)_$(pca_back2)_$(merging_back2)_clustering_labels.xls' \
	path_marker_selection='$(RWD_MERGED)/090_cytokine_bimatrix_frequencies_clustering/marker_selection.txt' path_cluster_merging=NULL aggregate_fun='median' scale=FALSE" $(RCODE)/02_heatmaps.R $(ROUT)/02_heatmaps.Rout
endef
$(foreach i,$(generated_clustering),$(eval $(call cytokine_bimatrix_frequencies_clustering_cytokine_profiles_heatmaps_rule,$(i))))



### ------------------------------------------------------------------------------
### Cytokine bimatrix frequencies top selection
### ------------------------------------------------------------------------------








































#
