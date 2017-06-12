### --------------------------------------------------------------------------
### Those parameters have to be defined
### --------------------------------------------------------------------------

R := R CMD BATCH --no-restore --no-save

RWD_MAIN := ../carsten_cytof/PD1_project
RCODE := .
METADATA := $(RWD_MAIN)/CK_metadata
PANELS := $(RWD_MAIN)/CK_panels

## arguments for the specific data set
data_dir_merged := CK_2016-06-merged_23_29/01

data_dir1 := CK_2016-06-23_01
data_dir2 := CK_2016-06-29_01
file_metadata1 := $(METADATA)/metadata_23_01.xlsx
file_metadata2 := $(METADATA)/metadata_29_01.xlsx

data1 := 23
data2 := 29
pca1 := pca1
pca2 := pca1
merging1 := merging6
merging2 := merging4

panel := 01

### --------------------------------------------------------------------------
### Those parameters do not need to be defined
### --------------------------------------------------------------------------

RWD1 := $(RWD_MAIN)/$(data_dir1)
RWD2 := $(RWD_MAIN)/$(data_dir2)
RWD_MERGED := $(RWD_MAIN)/$(data_dir_merged)
ROUT := $(RWD_MERGED)/Rout

responses := 3responses 2responses
days := both base

### --------------------------------------------------------------------------
## Define the default rule (makefiles are usually written so that the first target is for compiling the entire program)
.PHONY: all
all: mkdir_rout frequencies_goal expression_goal


### Make sure no intermediate files are deleted
.SECONDARY:

.PHONY: mkdir_rout
mkdir_rout:
	mkdir -p $(ROUT)

### --------------------------------------------------------------------------
### Frequency analysis
### --------------------------------------------------------------------------


.PHONY: frequencies_goal
frequencies_goal: $(RWD_MERGED)/050_frequencies_2responses_both/$(panel)_$(data1)$(merging1)_$(data2)$(merging2)_frequencies_plot_both2.pdf $(RWD_MERGED)/050_frequencies_2responses_both/$(panel)_$(data1)$(merging1)_$(data2)$(merging2)_frequencies_glmer_binomial_interglht_pheatmap3pvs_top05.pdf


$(RWD_MERGED)/050_frequencies_2responses_both/$(panel)_$(data1)$(merging1)_$(data2)$(merging2)_frequencies_plot_both2.pdf: $(RCODE)/04_frequencies_plot.R $(RCODE)/00_plot_frequencies.R $(file_metadata1) $(file_metadata2) $(RWD1)/050_frequencies/$(data1)_$(panel)_$(pca1)_$(merging1)_frequencies.xls $(RWD2)/050_frequencies/$(data2)_$(panel)_$(pca2)_$(merging2)_frequencies.xls
	echo "\n>>> 04_frequencies_plot"
	$(R) "--args prefix='$(panel)_$(data1)$(merging1)_$(data2)$(merging2)_' outdir='$(RWD_MERGED)/050_frequencies_2responses_both' path_metadata=c('$(file_metadata1)','$(file_metadata2)') path_frequencies=c('$(RWD1)/050_frequencies/$(data1)_$(panel)_$(pca1)_$(merging1)_frequencies.xls','$(RWD2)/050_frequencies/$(data2)_$(panel)_$(pca2)_$(merging2)_frequencies.xls') path_fun_plot_frequencies='$(RCODE)/00_plot_frequencies.R'" $(RCODE)/04_frequencies_plot.R $(ROUT)/04_frequencies_plot.Rout


$(RWD_MERGED)/050_frequencies_2responses_both/$(panel)_$(data1)$(merging1)_$(data2)$(merging2)_frequencies_glmer_binomial_interglht_pheatmap3pvs_top05.pdf: $(RCODE)/04_frequencies_analysis.R $(RCODE)/00_models.R $(RCODE)/00_formulas_2datasets_2responses_both.R $(RCODE)/00_plot_heatmaps_for_sign_freqs.R $(file_metadata1) $(file_metadata2) $(RWD1)/050_frequencies/$(data1)_$(panel)_$(pca1)_$(merging1)_counts.xls $(RWD2)/050_frequencies/$(data2)_$(panel)_$(pca2)_$(merging2)_counts.xls
	echo "\n>>> 04_frequencies_analysis"
	$(R) "--args prefix='$(panel)_$(data1)$(merging1)_$(data2)$(merging2)_' outdir='$(RWD_MERGED)/050_frequencies_2responses_both' path_metadata=c('$(file_metadata1)','$(file_metadata2)') path_counts=c('$(RWD1)/050_frequencies/$(data1)_$(panel)_$(pca1)_$(merging1)_counts.xls','$(RWD2)/050_frequencies/$(data2)_$(panel)_$(pca2)_$(merging2)_counts.xls') path_fun_models='$(RCODE)/00_models.R' path_fun_formulas='$(RCODE)/00_formulas_2datasets_2responses_both.R' path_fun_plot_heatmaps='$(RCODE)/00_plot_heatmaps_for_sign_freqs.R' FDR_cutoff='05'" $(RCODE)/04_frequencies_analysis.R $(ROUT)/04_frequencies_analysis.Rout


### --------------------------------------------------------------------------
### Marker expression analysis
### --------------------------------------------------------------------------

analysis_type := all clust

.PHONY: expression_goal
expression_goal: $(foreach i,$(analysis_type),$(RWD_MERGED)/080_expression_2responses_both/$(panel)_$(data1)$(merging1)_$(data2)$(merging2)_$(i)_expr_plot_both2.pdf) \
	$(foreach i,$(analysis_type),$(RWD_MERGED)/080_expression_2responses_both/$(panel)_$(data1)$(merging1)_$(data2)$(merging2)_$(i)_expr_lmer_interglht_pheatmap3pvs_top05.pdf)


define 04_expression_plot_rule
$(RWD_MERGED)/080_expression_2responses_both/$(panel)_$(data1)$(merging1)_$(data2)$(merging2)_$(1)_expr_plot_both2.pdf: $(RCODE)/04_expression_plot.R $(RCODE)/00_plot_expression.R $(file_metadata1) $(file_metadata2) $(RWD1)/080_expression/$(data1)_$(panel)_$(pca1)_$(merging1)_$(1)_expr.xls $(RWD2)/080_expression/$(data2)_$(panel)_$(pca2)_$(merging2)_$(1)_expr.xls
	echo "\n>>> 04_expression_plot"
	$(R) "--args prefix='$(panel)_$(data1)$(merging1)_$(data2)$(merging2)_$(1)_' outdir='$(RWD_MERGED)/080_expression_2responses_both' path_expression=c('$(RWD1)/080_expression/$(data1)_$(panel)_$(pca1)_$(merging1)_$(1)_expr.xls','$(RWD2)/080_expression/$(data2)_$(panel)_$(pca2)_$(merging2)_$(1)_expr.xls') path_metadata=c('$(file_metadata1)','$(file_metadata2)') path_fun_plot_expression='$(RCODE)/00_plot_expression.R'" $(RCODE)/04_expression_plot.R $(ROUT)/04_expression_plot.Rout
endef
$(foreach i,$(analysis_type),$(eval $(call 04_expression_plot_rule,$(i))))


define 04_expression_analysis_rule
$(RWD_MERGED)/080_expression_2responses_both/$(panel)_$(data1)$(merging1)_$(data2)$(merging2)_$(1)_expr_lmer_interglht_pheatmap3pvs_top05.pdf: $(RCODE)/04_expression_analysis.R $(RCODE)/00_models.R $(RCODE)/00_formulas_2datasets_2responses_both.R $(RCODE)/00_plot_heatmaps_for_sign_expr.R \
$(file_metadata1) $(file_metadata2) \
$(RWD1)/080_expression/$(data1)_$(panel)_$(pca1)_$(merging1)_$(1)_expr.xls $(RWD2)/080_expression/$(data2)_$(panel)_$(pca2)_$(merging2)_$(1)_expr.xls \
$(wildcard $(RWD_MERGED)/010_helpfiles/$(panel)_$(data1)$(merging1)_$(data2)$(merging2)_marker_exclusion.txt)
	echo "\n>>> 04_expression_analysis"
	$(R) "--args prefix='$(panel)_$(data1)$(merging1)_$(data2)$(merging2)_$(1)_' outdir='$(RWD_MERGED)/080_expression_2responses_both' path_expression=c('$(RWD1)/080_expression/$(data1)_$(panel)_$(pca1)_$(merging1)_$(1)_expr.xls','$(RWD2)/080_expression/$(data2)_$(panel)_$(pca2)_$(merging2)_$(1)_expr.xls') path_metadata=c('$(file_metadata1)','$(file_metadata2)') path_fun_models='$(RCODE)/00_models.R' path_fun_formulas='$(RCODE)/00_formulas_2datasets_2responses_both.R' \
	path_fun_plot_heatmaps='$(RCODE)/00_plot_heatmaps_for_sign_expr.R' path_marker_exclusion='$(RWD_MERGED)/010_helpfiles/$(panel)_$(data1)$(merging1)_$(data2)$(merging2)_marker_exclusion.txt' FDR_cutoff='05'" $(RCODE)/04_expression_analysis.R $(ROUT)/04_expression_analysis.Rout
endef
$(foreach i,$(analysis_type),$(eval $(call 04_expression_analysis_rule,$(i))))















#
