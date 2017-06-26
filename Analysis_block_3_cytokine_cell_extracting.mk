### --------------------------------------------------------------------------
### Those parameters have to be defined
### --------------------------------------------------------------------------

R := R CMD BATCH --no-restore --no-save

RWD_MAIN := ../carsten_cytof/PD1_project
RCODE := .
METADATA := $(RWD_MAIN)/CK_metadata
PANELS := $(RWD_MAIN)/CK_panels

## arguments for the specific data set
data_dir := CK_2016-06-23_02_CD4_merging3
file_metadata := $(METADATA)/metadata_23_02.xlsx
file_extract_marker := $(PANELS)/panel2CD4_23_cytokines_CD69.xlsx

data := 23CD4
panel := 02CD4v2
pca := pca0
merging := merging1
extract_cluster := Tmem
extract_marker := CD69


### --------------------------------------------------------------------------
### Those parameters do not need to be defined
### --------------------------------------------------------------------------

RWD := $(RWD_MAIN)/$(data_dir)
ROUT := $(RWD)/Rout

make_file := Analysis_block_3_cytokine_cell_extracting.mk

### --------------------------------------------------------------------------
## Define the default rule (makefiles are usually written so that the first target is for compiling the entire program)
.PHONY: all
all: cytokine_cell_extracting_goal


### Make sure no intermediate files are deleted
.SECONDARY:


### --------------------------------------------------------------------------
### Cytokine cell extracting
### --------------------------------------------------------------------------


.PHONY: cytokine_cell_extracting_goal
cytokine_cell_extracting_goal: $(RWD)_$(extract_cluster)_$(merging)_$(extract_marker)/010_cleanfcs/*.fcs

$(RWD)_$(extract_cluster)_$(merging)_$(extract_marker)/010_cleanfcs/*.fcs: $(RCODE)/02_cytokine_cell_extracting.R $(file_metadata) $(RWD)/010_cleanfcs/*.fcs $(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_$(merging)_clustering.xls $(file_extract_marker)
	echo "\n>> $(make_file)\n>>> 02_cytokine_cell_extracting"
	$(R) "--args outdir='$(RWD)_$(extract_cluster)_$(merging)_$(extract_marker)/010_cleanfcs' dir_fcs='$(RWD)/010_cleanfcs' path_metadata='$(file_metadata)' path_clustering='$(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_$(merging)_clustering.xls'  path_clustering_labels='$(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_$(merging)_clustering_labels.xls' path_extract_cluster='$(RWD)/010_helpfiles/$(data)_$(panel)_$(pca)_$(merging)_extract_cluster_$(extract_cluster).txt' path_extract_marker='$(file_extract_marker)'" $(RCODE)/02_cytokine_cell_extracting.R $(ROUT)/02_cytokine_cell_extracting.Rout





























#
