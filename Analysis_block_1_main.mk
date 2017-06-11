### --------------------------------------------------------------------------
### Those parameters have to be defined
### --------------------------------------------------------------------------

R := R CMD BATCH --no-restore --no-save

RWD_MAIN := ../carsten_cytof/PD1_project
RCODE := .
METADATA := $(RWD_MAIN)/CK_metadata
PANELS := $(RWD_MAIN)/CK_panels


## arguments for the specific data set
data_dir := CK_2016-06-23_01
file_panel := $(PANELS)/panel1.xlsx
file_metadata := $(METADATA)/metadata_23_01.xlsx

data := 23
panel := 01
pca := pca1
nmetaclusts := 20

rand_seed_consensus := 123

### --------------------------------------------------------------------------
### Those parameters do not need to be defined
### --------------------------------------------------------------------------


## Global parameters
tsne_pmin := 1500

RWD := $(RWD_MAIN)/$(data_dir)
ROUT := $(RWD)/Rout


### --------------------------------------------------------------------------
## Define the default rule (makefiles are usually written so that the first target is for compiling the entire program)
.PHONY: all
all: mkdir_rout data_normalization_goal data_plot_goal pcascores_goal select_observables_goal flowsom_goal heatmaps_goal runtsne_goal plottsne_goal frequencies_calculate_goal frequencies_goal expression_calculate_goal expression_goal


### Make sure no intermediate files are deleted
.SECONDARY:

.PHONY: mkdir_rout
mkdir_rout:
	mkdir -p $(ROUT)


### --------------------------------------------------------------------------
### Data normalization
### --------------------------------------------------------------------------

.PHONY: data_normalization_goal
data_normalization_goal: $(RWD)/010_data/$(data)_$(panel)_expr_raw.rds

$(RWD)/010_data/$(data)_$(panel)_expr_raw.rds: $(RCODE)/01_data_normalization.R $(file_metadata) $(file_panel) $(RWD)/010_cleanfcs/*.fcs
	echo "\n>>> 01_data_normalization"
	$(R) "--args dir_fcs='$(RWD)/010_cleanfcs' path_metadata='$(file_metadata)' path_panel='$(file_panel)' outdir='$(RWD)/010_data' prefix='$(data)_$(panel)_'" $(RCODE)/01_data_normalization.R $(ROUT)/01_data_normalization.Rout


.PHONY: data_plot_goal
data_plot_goal: $(RWD)/010_data/$(data)_$(panel)_distrosgrp_raw.pdf $(RWD)/010_data/$(data)_$(panel)_distrosgrp_norm.pdf

### For raw data
$(RWD)/010_data/$(data)_$(panel)_distrosgrp_raw.pdf: $(RCODE)/01_data_plot.R $(file_metadata) $(file_panel) $(RWD)/010_data/$(data)_$(panel)_expr_raw.rds
	echo "\n>>> 01_data_plot"
	$(R) "--args path_data='$(RWD)/010_data/$(data)_$(panel)_expr_raw.rds' path_metadata='$(file_metadata)' path_panel='$(file_panel)' outdir='$(RWD)/010_data' prefix='$(data)_$(panel)_' suffix='_raw'" $(RCODE)/01_data_plot.R $(ROUT)/01_data_plot.Rout

### For norm data
$(RWD)/010_data/$(data)_$(panel)_distrosgrp_norm.pdf: $(RCODE)/01_data_plot.R $(file_metadata) $(file_panel) $(RWD)/010_data/$(data)_$(panel)_expr_norm.rds
	echo "\n>>> 01_data_plot"
	$(R) "--args path_data='$(RWD)/010_data/$(data)_$(panel)_expr_norm.rds' path_metadata='$(file_metadata)' path_panel='$(file_panel)' outdir='$(RWD)/010_data' prefix='$(data)_$(panel)_' suffix='_norm'" $(RCODE)/01_data_plot.R $(ROUT)/01_data_plot.Rout


### --------------------------------------------------------------------------
### PCA scores based on raw data
### --------------------------------------------------------------------------

.PHONY: pcascores_goal
pcascores_goal: $(RWD)/020_pcascores/$(data)_$(panel)_princompscore_by_sample.xls

$(RWD)/020_pcascores/$(data)_$(panel)_princompscore_by_sample.xls: $(RCODE)/01_pcascores.R $(RWD)/010_data/$(data)_$(panel)_expr_raw.rds
	echo "\n>>> 01_pcascores"
	$(R) "--args prefix='$(data)_$(panel)_' outdir='$(RWD)/020_pcascores' path_data='$(RWD)/010_data/$(data)_$(panel)_expr_raw.rds' path_metadata='$(file_metadata)' path_panel='$(file_panel)'" $(RCODE)/01_pcascores.R $(ROUT)/01_pcascores.Rout


### --------------------------------------------------------------------------
### Select observables for clustering
### --------------------------------------------------------------------------

.PHONY: select_observables_goal
select_observables_goal: $(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_clustering_observables.xls


$(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_clustering_observables.xls: $(RCODE)/01_select_observables.R $(RWD)/020_pcascores/$(data)_$(panel)_princompscore_by_sample.xls $(RWD)/010_helpfiles/$(data)_$(panel)_$(pca)_princompscore_cutoff.txt
	echo "\n>>> 02_select_observables"
	$(R) "--args prefix='$(data)_$(panel)_$(pca)_' outdir='$(RWD)/030_heatmaps' path_pca_score='$(RWD)/020_pcascores/$(data)_$(panel)_princompscore_by_sample.xls' path_pca_score_cutoff='$(RWD)/010_helpfiles/$(data)_$(panel)_$(pca)_princompscore_cutoff.txt'" $(RCODE)/01_select_observables.R $(ROUT)/01_select_observables.Rout


### --------------------------------------------------------------------------
### FlowSOM clustering based on raw data - average linkage
### --------------------------------------------------------------------------

.PHONY: flowsom_goal
flowsom_goal: $(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_clustering.xls

$(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_clustering.xls: $(RCODE)/02_flowsom.R $(RWD)/010_data/$(data)_$(panel)_expr_raw.rds $(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_clustering_observables.xls
	echo "\n>>> 02_flowsom"
	$(R) "--args prefix='$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_' outdir='$(RWD)/030_heatmaps' path_data='$(RWD)/010_data/$(data)_$(panel)_expr_raw.rds' path_clustering_observables='$(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_clustering_observables.xls' nmetaclusts=$(nmetaclusts) rand_seed_consensus=$(rand_seed_consensus)" $(RCODE)/02_flowsom.R $(ROUT)/02_flowsom.Rout



### --------------------------------------------------------------------------
### Heatmaps
### --------------------------------------------------------------------------


.PHONY: heatmaps_goal
heatmaps_goal: $(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_pheatmap_all_all_row_clust_raw.pdf

$(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_pheatmap_all_all_row_clust_raw.pdf: $(RCODE)/02_heatmaps.R $(RWD)/010_data/$(data)_$(panel)_expr_raw.rds $(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_clustering.xls $(wildcard $(RWD)/010_helpfiles/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_marker_selection.txt)
	echo "\n>>> 02_heatmaps"
	$(R) "--args prefix='$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_' outdir='$(RWD)/030_heatmaps' path_data='$(RWD)/010_data/$(data)_$(panel)_expr_raw.rds' path_data_norm='$(RWD)/010_data/$(data)_$(panel)_expr_norm.rds' \
	path_clustering_observables='$(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_clustering_observables.xls' path_clustering='$(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_clustering.xls'  path_clustering_labels='$(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_clustering_labels.xls' \
	path_marker_selection='$(RWD)/010_helpfiles/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_marker_selection.txt' path_cluster_merging=NULL" $(RCODE)/02_heatmaps.R $(ROUT)/02_heatmaps.Rout


### --------------------------------------------------------------------------
### t-SNE
### --------------------------------------------------------------------------

.PHONY: runtsne_goal
runtsne_goal: $(RWD)/040_tsnemaps/$(data)_$(panel)_$(pca)_rtsne_out.rds

$(RWD)/040_tsnemaps/$(data)_$(panel)_$(pca)_rtsne_out.rds: $(RCODE)/03_runtsne.R $(RWD)/010_data/$(data)_$(panel)_expr_raw.rds $(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_clustering_observables.xls
	echo "\n>>> 03_runtsne"
	$(R) "--args prefix='$(data)_$(panel)_$(pca)_' outdir='$(RWD)/040_tsnemaps' path_data='$(RWD)/010_data/$(data)_$(panel)_expr_raw.rds' path_clustering_observables='$(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_clustering_observables.xls' tsne_pmin=$(tsne_pmin)" $(RCODE)/03_runtsne.R $(ROUT)/03_runtsne.Rout

.PHONY: plottsne_goal
plottsne_goal: $(RWD)/040_tsnemaps/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_tSNEone.pdf

$(RWD)/040_tsnemaps/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_tSNEone.pdf: $(RCODE)/03_plottsne.R $(RWD)/040_tsnemaps/$(data)_$(panel)_$(pca)_rtsne_out.rds $(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_clustering.xls $(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_clustering_labels.xls
	echo "\n>>> 03_plottsne"
	$(R) "--args prefix='$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_' outdir='$(RWD)/040_tsnemaps' path_metadata='$(file_metadata)'  path_rtsne_out='$(RWD)/040_tsnemaps/$(data)_$(panel)_$(pca)_rtsne_out.rds' \
	path_rtsne_data='$(RWD)/040_tsnemaps/$(data)_$(panel)_$(pca)_rtsne_data.xls' path_clustering='$(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_clustering.xls' path_clustering_labels='$(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_clustering_labels.xls'" $(RCODE)/03_plottsne.R $(ROUT)/03_plottsne.Rout


### --------------------------------------------------------------------------
### Frequency analysis
### --------------------------------------------------------------------------

.PHONY: frequencies_calculate_goal
frequencies_calculate_goal: $(RWD)/050_frequencies_auto/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_frequencies.xls

$(RWD)/050_frequencies_auto/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_frequencies.xls: $(RCODE)/04_frequencies_calculate.R $(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_clustering.xls
	echo "\n>>> 04_frequencies_calculate"
	$(R) "--args prefix='$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_' outdir='$(RWD)/050_frequencies_auto' path_clustering='$(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_clustering.xls' path_clustering_labels='$(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_clustering_labels.xls'" $(RCODE)/04_frequencies_calculate.R $(ROUT)/04_frequencies_calculate.Rout


.PHONY: frequencies_goal
frequencies_goal: $(RWD)/050_frequencies_auto/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_frequencies_plot_both2.pdf $(RWD)/050_frequencies_auto/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_frequencies_glmer_binomial_interglht_pheatmap3pvs_top05.pdf


$(RWD)/050_frequencies_auto/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_frequencies_plot_both2.pdf: $(RCODE)/04_frequencies_plot.R $(RCODE)/00_plot_frequencies.R $(file_metadata) $(RWD)/050_frequencies_auto/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_frequencies.xls
	echo "\n>>> 04_frequencies_plot"
	$(R) "--args prefix='$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_' outdir='$(RWD)/050_frequencies_auto' path_metadata='$(file_metadata)' path_frequencies='$(RWD)/050_frequencies_auto/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_frequencies.xls' path_fun_plot_frequencies='$(RCODE)/00_plot_frequencies.R'" $(RCODE)/04_frequencies_plot.R $(ROUT)/04_frequencies_plot.Rout


$(RWD)/050_frequencies_auto/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_frequencies_glmer_binomial_interglht_pheatmap3pvs_top05.pdf: $(RCODE)/04_frequencies_analysis.R $(RCODE)/00_models.R $(RCODE)/00_formulas_1dataset_3responses_both.R $(RCODE)/00_plot_heatmaps_for_sign_freqs.R $(RWD)/050_frequencies_auto/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_counts.xls $(file_metadata)
	echo "\n>>> 04_frequencies_analysis"
	$(R) "--args prefix='$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_' outdir='$(RWD)/050_frequencies_auto' path_metadata='$(file_metadata)' path_counts='$(RWD)/050_frequencies_auto/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_counts.xls' path_fun_models='$(RCODE)/00_models.R' path_fun_formulas='$(RCODE)/00_formulas_1dataset_3responses_both.R' path_fun_plot_heatmaps='$(RCODE)/00_plot_heatmaps_for_sign_freqs.R' FDR_cutoff='05'" $(RCODE)/04_frequencies_analysis.R $(ROUT)/04_frequencies_analysis.Rout


### --------------------------------------------------------------------------
### Marker expression analysis
### --------------------------------------------------------------------------

.PHONY: expression_calculate_goal
expression_calculate_goal: $(RWD)/080_expression_auto/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_all_expr.xls

$(RWD)/080_expression_auto/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_all_expr.xls: $(RCODE)/04_expression_calculate.R $(RWD)/010_data/$(data)_$(panel)_expr_raw.rds $(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_clustering_observables.xls $(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_clustering.xls
	echo "\n>>> 04_expression_calculate"
	$(R) "--args prefix='$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_' outdir='$(RWD)/080_expression_auto' path_data='$(RWD)/010_data/$(data)_$(panel)_expr_raw.rds' \
	path_clustering_observables='$(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_clustering_observables.xls' path_clustering='$(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_clustering.xls'  path_clustering_labels='$(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_clustering_labels.xls'" $(RCODE)/04_expression_calculate.R $(ROUT)/04_expression_calculate.Rout


analysis_type := all clust

.PHONY: expression_goal
expression_goal: $(foreach i,$(analysis_type),$(RWD)/080_expression_auto/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_$(i)_expr_plot_both2.pdf) \
	$(foreach i,$(analysis_type),$(RWD)/080_expression_auto/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_$(i)_expr_lmer_interglht_pheatmap3pvs_top05.pdf)


define 04_expression_plot_rule
$(RWD)/080_expression_auto/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_$(1)_expr_plot_both2.pdf: $(RCODE)/04_expression_plot.R $(RCODE)/00_plot_expression.R $(RWD)/080_expression_auto/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_$(1)_expr.xls
	echo "\n>>> 04_expression_plot"
	$(R) "--args prefix='$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_$(1)_' outdir='$(RWD)/080_expression_auto' path_expression='$(RWD)/080_expression_auto/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_$(1)_expr.xls' path_metadata='$(file_metadata)' path_fun_plot_expression='$(RCODE)/00_plot_expression.R'" $(RCODE)/04_expression_plot.R $(ROUT)/04_expression_plot.Rout
endef
$(foreach i,$(analysis_type),$(eval $(call 04_expression_plot_rule,$(i))))


define 04_expression_analysis_rule
$(RWD)/080_expression_auto/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_$(1)_expr_lmer_interglht_pheatmap3pvs_top05.pdf: $(RCODE)/04_expression_analysis.R $(RCODE)/00_models.R $(RCODE)/00_formulas_1dataset_3responses_both.R $(RCODE)/00_plot_heatmaps_for_sign_expr.R $(RWD)/080_expression_auto/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_$(1)_expr.xls $(wildcard $(RWD)/010_helpfiles/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_marker_exclusion.txt)
	echo "\n>>> 04_expression_analysis"
	$(R) "--args prefix='$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_$(1)_' outdir='$(RWD)/080_expression_auto' path_expression='$(RWD)/080_expression_auto/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_$(1)_expr.xls' path_metadata='$(file_metadata)' path_fun_models='$(RCODE)/00_models.R' path_fun_formulas='$(RCODE)/00_formulas_1dataset_3responses_both.R' \
	path_fun_plot_heatmaps='$(RCODE)/00_plot_heatmaps_for_sign_expr.R' path_marker_exclusion='$(RWD)/010_helpfiles/$(data)_$(panel)_$(pca)_cl$(nmetaclusts)_marker_exclusion.txt' FDR_cutoff='05'" $(RCODE)/04_expression_analysis.R $(ROUT)/04_expression_analysis.Rout
endef
$(foreach i,$(analysis_type),$(eval $(call 04_expression_analysis_rule,$(i))))




















#
