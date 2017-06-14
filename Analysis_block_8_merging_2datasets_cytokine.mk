### --------------------------------------------------------------------------
### Those parameters have to be defined
### --------------------------------------------------------------------------

R := R CMD BATCH --no-restore --no-save

RWD_MAIN := ../carsten_cytof/PD1_project
RCODE := .
METADATA := $(RWD_MAIN)/CK_metadata
PANELS := $(RWD_MAIN)/CK_panels

## arguments for the specific data set
data_dir_merged := CK_2016-06-merged_23_29/02

data_dir1 := CK_2016-06-23_02_CD4_merging2_Tmem_merging2_CD69
data_dir2 := CK_2016-06-23_02_CD4_merging2_Tmem_merging2_CD69
file_metadata1 := $(METADATA)/metadata_23_02.xlsx
file_metadata2 := $(METADATA)/metadata_23_02.xlsx

data1 := 23CD4TmemCD69
data2 := 23CD4TmemCD69

cytokines := 02CD4

som_dim := 7
nmetaclusts := 49


### --------------------------------------------------------------------------
### Those parameters do not need to be defined
### --------------------------------------------------------------------------

RWD1 := $(RWD_MAIN)/$(data_dir1)
RWD2 := $(RWD_MAIN)/$(data_dir2)

RWD_MERGED := $(RWD_MAIN)/$(data_dir_merged)

ROUT := $(RWD)/Rout

rand_seed_consensus := 1234

### --------------------------------------------------------------------------
## Define the default rule (makefiles are usually written so that the first target is for compiling the entire program)
.PHONY: all
all: mkdir_rout cytokine_expression_goal cytokine_bimatrix_frequencies_overall_goal cytokine_bimatrix_frequencies_clustering_goal


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
cytokine_expression_goal: $(foreach i,$(analysis_type),$(RWD)/090_cytokine_expression/$(data)_$(cytokines)_$(i)_expr_plot_both2.pdf) \
	$(foreach i,$(analysis_type),$(RWD)/090_cytokine_expression/$(data)_$(cytokines)_$(i)_expr_lmer_interglht_pheatmap3pvs_top05.pdf)

define 04_expression_plot_rule
$(RWD)/090_cytokine_expression/$(data)_$(cytokines)_$(1)_expr_plot_both2.pdf: $(RCODE)/04_expression_plot.R $(RCODE)/00_plot_expression.R $(RWD)/090_cytokine_expression/$(data)_$(cytokines)_$(1)_expr.xls
	echo "\n>>> 04_expression_plot"
	$(R) "--args prefix='$(data)_$(cytokines)_$(1)_' outdir='$(RWD)/090_cytokine_expression' path_expression='$(RWD)/090_cytokine_expression/$(data)_$(cytokines)_$(1)_expr.xls' path_metadata='$(file_metadata)' path_fun_plot_expression='$(RCODE)/00_plot_expression.R'" $(RCODE)/04_expression_plot.R $(ROUT)/04_expression_plot.Rout
endef
$(foreach i,$(analysis_type),$(eval $(call 04_expression_plot_rule,$(i))))


define 04_expression_analysis_rule
$(RWD)/090_cytokine_expression/$(data)_$(cytokines)_$(1)_expr_lmer_interglht_pheatmap3pvs_top05.pdf: $(RCODE)/04_expression_analysis.R $(RCODE)/00_models.R $(RCODE)/00_formulas_2datasets_3responses_both.R $(RCODE)/00_plot_heatmaps_for_sign_expr.R $(RWD)/090_cytokine_expression/$(data)_$(cytokines)_$(1)_expr.xls $(wildcard $(RWD)/010_helpfiles/$(data)_$(cytokines)_$(1)_marker_exclusion.txt)
	echo "\n>>> 04_expression_analysis"
	$(R) "--args prefix='$(data)_$(cytokines)_$(1)_' outdir='$(RWD)/090_cytokine_expression' path_expression='$(RWD)/090_cytokine_expression/$(data)_$(cytokines)_$(1)_expr.xls' path_metadata='$(file_metadata)' path_fun_models='$(RCODE)/00_models.R' path_fun_formulas='$(RCODE)/00_formulas_2datasets_3responses_both.R' \
	path_fun_plot_heatmaps='$(RCODE)/00_plot_heatmaps_for_sign_expr.R' path_marker_exclusion='$(RWD)/010_helpfiles/$(data)_$(cytokines)_$(1)_marker_exclusion.txt' FDR_cutoff='05'" $(RCODE)/04_expression_analysis.R $(ROUT)/04_expression_analysis.Rout
endef
$(foreach i,$(analysis_type),$(eval $(call 04_expression_analysis_rule,$(i))))


### --------------------------------------------------------------------------
### Cytokine bimatrix analysis
### --------------------------------------------------------------------------


### ---------------------------------------
### Cytokine bimatrix frequencies overall
### ---------------------------------------

.PHONY: cytokine_bimatrix_frequencies_overall_goal
cytokine_bimatrix_frequencies_overall_goal: $(RWD)/090_cytokine_bimatrix_frequencies_overall/$(data)_$(cytokines)_frequencies_plot_both2.pdf  $(RWD)/090_cytokine_bimatrix_frequencies_overall/$(data)_$(cytokines)_frequencies_glmer_binomial_interglht_pheatmap3pvs_top05.pdf


$(RWD)/090_cytokine_bimatrix_frequencies_overall/$(data)_$(cytokines)_frequencies_plot_both2.pdf: $(RCODE)/04_frequencies_plot.R $(RCODE)/00_plot_frequencies.R $(file_metadata) $(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_frequencies.xls
	echo "\n>>> 04_frequencies_plot"
	$(R) "--args prefix='$(data)_$(cytokines)_' outdir='$(RWD)/090_cytokine_bimatrix_frequencies_overall' path_metadata='$(file_metadata)' path_frequencies='$(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_frequencies.xls' path_fun_plot_frequencies='$(RCODE)/00_plot_frequencies.R'" $(RCODE)/04_frequencies_plot.R $(ROUT)/04_frequencies_plot.Rout


$(RWD)/090_cytokine_bimatrix_frequencies_overall/$(data)_$(cytokines)_frequencies_glmer_binomial_interglht_pheatmap3pvs_top05.pdf: $(RCODE)/05_cytokine_bimatrix_overall_frequencies_analysis.R $(RCODE)/00_models.R $(RCODE)/00_formulas_2datasets_3responses_both.R $(RCODE)/00_plot_heatmaps_for_sign_freqs.R $(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_counts.xls $(file_metadata)
	echo "\n>>> 05_cytokine_bimatrix_overall_frequencies_analysis"
	$(R) "--args prefix='$(data)_$(cytokines)_' outdir='$(RWD)/090_cytokine_bimatrix_frequencies_overall' path_metadata='$(file_metadata)' path_counts='$(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_counts.xls' path_fun_models='$(RCODE)/00_models.R' path_fun_formulas='$(RCODE)/00_formulas_2datasets_3responses_both.R' path_fun_plot_heatmaps='$(RCODE)/00_plot_heatmaps_for_sign_freqs.R' FDR_cutoff='05'" $(RCODE)/05_cytokine_bimatrix_overall_frequencies_analysis.R $(ROUT)/05_cytokine_bimatrix_overall_frequencies_analysis.Rout


### ---------------------------------------
### Cytokine bimatrix frequencies clustering
### ---------------------------------------


.PHONY: cytokine_bimatrix_frequencies_clustering_goal
cytokine_bimatrix_frequencies_clustering_goal: $(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_clustering.xls $(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_pheatmap_all_all_row_clust_raw.pdf frequencies_calculate_goal1 frequencies_goal1


# FlowSOM clustering of the bimatrix

$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_clustering.xls: $(RCODE)/02_flowsom.R $(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_bimatrix.rds $(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_clustering_observables.xls
	echo "\n>>> 02_flowsom"
	$(R) "--args prefix='$(data)_$(cytokines)_cl$(nmetaclusts)_' outdir='$(RWD)/090_cytokine_bimatrix_frequencies_clustering' path_data='$(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_bimatrix.rds' path_clustering_observables='$(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_clustering_observables.xls' nmetaclusts=$(nmetaclusts) rand_seed_consensus=$(rand_seed_consensus) som_dim=$(som_dim)" $(RCODE)/02_flowsom.R $(ROUT)/02_flowsom.Rout

# Plot a heatmap of the bimatrix clustering

$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_pheatmap_all_all_row_clust_raw.pdf: $(RCODE)/02_heatmaps.R $(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_bimatrix.rds $(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_clustering.xls
	echo "\n>>> 02_heatmaps"
	$(R) "--args prefix='$(data)_$(cytokines)_cl$(nmetaclusts)_' outdir='$(RWD)/090_cytokine_bimatrix_frequencies_clustering' path_data='$(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_bimatrix.rds' path_data_norm=NULL \
	path_clustering_observables='$(RWD)/090_cytokine_bimatrix/$(data)_$(cytokines)_clustering_observables.xls' path_clustering='$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_clustering.xls'  path_clustering_labels='$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_clustering_labels.xls' \
	path_marker_selection=NULL path_cluster_merging=NULL aggregate_fun='mean' scale=FALSE" $(RCODE)/02_heatmaps.R $(ROUT)/02_heatmaps.Rout



.PHONY: frequencies_calculate_goal1
frequencies_calculate_goal1: $(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_frequencies.xls

# Calculate frequancies of the bimatrix clusters


$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_frequencies.xls: $(RCODE)/04_frequencies_calculate.R $(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_clustering.xls
	echo "\n>>> 04_frequencies_calculate"
	$(R) "--args prefix='$(data)_$(cytokines)_cl$(nmetaclusts)_' outdir='$(RWD)/090_cytokine_bimatrix_frequencies_clustering' path_clustering='$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_clustering.xls' path_clustering_labels='$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_clustering_labels.xls'" $(RCODE)/04_frequencies_calculate.R $(ROUT)/04_frequencies_calculate.Rout



model := glmer_binomial_interglht


.PHONY: frequencies_goal1
frequencies_goal1: $(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_frequencies_$(model)_pheatmap3pvs_top05.pdf $(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_top_$(model)_NRvsR_base_frequencies_plot_both2.pdf


# Differential analysis of the frequancies of the bimatrix clusters

$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_frequencies_$(model)_pheatmap3pvs_top05.pdf: $(RCODE)/04_frequencies_analysis.R $(RCODE)/00_models.R $(RCODE)/00_formulas_2datasets_3responses_both.R $(RCODE)/00_plot_heatmaps_for_sign_freqs.R $(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_counts.xls $(file_metadata)
	echo "\n>>> 04_frequencies_analysis"
	$(R) "--args prefix='$(data)_$(cytokines)_cl$(nmetaclusts)_' outdir='$(RWD)/090_cytokine_bimatrix_frequencies_clustering' path_metadata='$(file_metadata)' path_counts='$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_counts.xls' path_fun_models='$(RCODE)/00_models.R' path_fun_formulas='$(RCODE)/00_formulas_2datasets_3responses_both.R' path_fun_plot_heatmaps='$(RCODE)/00_plot_heatmaps_for_sign_freqs.R' FDR_cutoff='05'" $(RCODE)/04_frequencies_analysis.R $(ROUT)/04_frequencies_analysis.Rout


# Plot the significant bimatrix clusters

$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_top_$(model)_NRvsR_base_frequencies_plot_both2.pdf: $(RCODE)/04_frequencies_plot_significant.R $(RCODE)/00_plot_frequencies.R $(file_metadata) $(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_frequencies.xls $(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_frequencies_pvs_$(model)_top05.xls
	echo "\n>>> 04_frequencies_plot_significant"
	$(R) "--args prefix='$(data)_$(cytokines)_cl$(nmetaclusts)_top_$(model)_' outdir='$(RWD)/090_cytokine_bimatrix_frequencies_clustering' path_metadata='$(file_metadata)' path_frequencies='$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_frequencies.xls' path_pvs='$(RWD)/090_cytokine_bimatrix_frequencies_clustering/$(data)_$(cytokines)_cl$(nmetaclusts)_frequencies_pvs_$(model)_top05.xls' path_fun_plot_frequencies='$(RCODE)/00_plot_frequencies.R' FDR_cutoff='05'" $(RCODE)/04_frequencies_plot_significant.R $(ROUT)/04_frequencies_plot_significant.Rout









### Cytokine bimatrix frequencies top selection









































#
