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

data_dir1 := CK_2016-06-23_02_CD4_merging3_Tmem_merging1_CD69
data_dir2 := CK_2016-06-29_02_CD4_merging3_Tmem_merging1_CD69
file_metadata1 := $(METADATA)/metadata_23_02.xlsx
file_metadata2 := $(METADATA)/metadata_29_02.xlsx

data1 := 23CD4TmemCD69
data2 := 29CD4TmemCD69

cytokines := 02CD4

som_dim := 7
nmetaclusts := 49


### --------------------------------------------------------------------------
### Those parameters do not need to be defined
### --------------------------------------------------------------------------

RWD1 := $(RWD_MAIN)/$(data_dir1)
RWD2 := $(RWD_MAIN)/$(data_dir2)
RWD_MERGED := $(RWD_MAIN)/$(data_dir_merged)
ROUT := $(RWD_MERGED)/Rout

rand_seed_consensus := 1234
make_file := Analysis_block_8_merging_2datasets_cytokine.mk

### --------------------------------------------------------------------------
## Define the default rule (makefiles are usually written so that the first target is for compiling the entire program)
.PHONY: all
all: mkdir_rout cytokine_expression_goal cytokine_bimatrix_frequencies_calculate_goal cytokine_bimatrix_frequencies_overall_goal cytokine_bimatrix_frequencies_clustering_goal cytokine_bimatrix_frequencies_clustering_calculate_goal cytokine_bimatrix_frequencies_clustering_analysis_goal cytokine_bimatrix_frequencies_clustering_plot_significant_goal


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


### ---------------------------------------
### Cytokine bimatrix frequencies clustering
### ---------------------------------------


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





### Cytokine bimatrix frequencies top selection









































#
