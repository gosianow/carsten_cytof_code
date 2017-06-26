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
file_metadata := $(METADATA)/metadata_23_01.xlsx

data := 23
panel := 01
pca := pca1
merging := mergingNEW2
extract_cluster := CD4


### --------------------------------------------------------------------------
### Those parameters do not need to be defined
### --------------------------------------------------------------------------

RWD := $(RWD_MAIN)/$(data_dir)
ROUT := $(RWD)/Rout

make_file := Analysis_block_3_cluster_extracting.mk

### --------------------------------------------------------------------------
## Define the default rule (makefiles are usually written so that the first target is for compiling the entire program)
.PHONY: all
all: cluster_extracting_goal


### Make sure no intermediate files are deleted
.SECONDARY:


### --------------------------------------------------------------------------
### Cluster extracting
### --------------------------------------------------------------------------

.PHONY: cluster_extracting_goal
cluster_extracting_goal: $(RWD)_$(extract_cluster)_$(merging)/010_cleanfcs/*.fcs

$(RWD)_$(extract_cluster)_$(merging)/010_cleanfcs/*.fcs: $(RCODE)/02_cluster_extracting.R $(file_metadata) $(RWD)/010_cleanfcs/*.fcs $(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_$(merging)_clustering.xls
	echo "\n>> $(make_file)\n>>> 02_cluster_extracting"
	$(R) "--args outdir='$(RWD)_$(extract_cluster)_$(merging)/010_cleanfcs' dir_fcs='$(RWD)/010_cleanfcs' path_metadata='$(file_metadata)'  path_clustering='$(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_$(merging)_clustering.xls'  path_clustering_labels='$(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_$(merging)_clustering_labels.xls' extract_cluster='$(extract_cluster)'" $(RCODE)/02_cluster_extracting.R $(ROUT)/02_cluster_extracting.Rout

















#
