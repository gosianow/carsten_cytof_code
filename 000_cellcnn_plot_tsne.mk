### --------------------------------------------------------------------------
### Those parameters have to be defined
### --------------------------------------------------------------------------

R := R CMD BATCH --no-restore --no-save

RWD_MAIN := ../carsten_cytof/PD1_project
RWD_CELLCNN := ../PD1_CellCnn_Lukas_June2017
RCODE := .
METADATA := $(RWD_MAIN)/CK_metadata
PANELS := $(RWD_MAIN)/CK_panels


## arguments for the specific data set
data_dir := CK_2016-06-23_01
file_metadata := $(METADATA)/metadata_23_01.xlsx

data := 23
panel := 01
pca := pca1
merging := merging6

cellcnn_files_data := panel1_base
cellcnn_files_type := combined

cellcnn_files_dir := $(cellcnn_files_data)_$(cellcnn_files_type)

day := base

### --------------------------------------------------------------------------
### Those parameters do not need to be defined
### --------------------------------------------------------------------------

RWD := $(RWD_MAIN)/$(data_dir)
ROUT := $(RWD)/Rout

### --------------------------------------------------------------------------
## Define the default rule (makefiles are usually written so that the first target is for compiling the entire program)
.PHONY: all
all: mkdir_rout plottsne_goal


### Make sure no intermediate files are deleted
.SECONDARY:

.PHONY: mkdir_rout
mkdir_rout:
	mkdir -p $(ROUT)


### --------------------------------------------------------------------------
### t-SNE
### --------------------------------------------------------------------------


.PHONY: plottsne_goal
plottsne_goal: $(RWD_MAIN)/cellcnn_tsne/$(data)_$(panel)_$(pca)_$(merging)_$(cellcnn_files_type)_tSNEone_filter_0.pdf

$(RWD_MAIN)/cellcnn_tsne/$(data)_$(panel)_$(pca)_$(merging)_$(cellcnn_files_type)_tSNEone_filter_0.pdf: $(RCODE)/000_cellcnn_plot_tsne.R $(RWD)/040_tsnemaps/$(data)_$(panel)_$(pca)_rtsne_out.rds $(file_metadata) $(RWD_CELLCNN)/$(cellcnn_files_dir)/selected_cells/*.csv
	echo "\n>> $(make_file)\n>>> 000_cellcnn_plot_tsne"
	$(R) "--args prefix='$(data)_$(panel)_$(pca)_$(merging)_$(cellcnn_files_type)_' outdir='$(RWD_MAIN)/cellcnn_tsne' path_metadata='$(file_metadata)' \
	path_rtsne_out='$(RWD)/040_tsnemaps/$(data)_$(panel)_$(pca)_rtsne_out.rds' path_rtsne_data='$(RWD)/040_tsnemaps/$(data)_$(panel)_$(pca)_rtsne_data.xls' \
	path_clustering='$(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_$(merging)_clustering.xls' path_clustering_labels='$(RWD)/030_heatmaps/$(data)_$(panel)_$(pca)_$(merging)_clustering_labels.xls' dir_cellcnn_files='$(RWD_CELLCNN)/$(cellcnn_files_dir)/selected_cells' day='$(day)'" $(RCODE)/000_cellcnn_plot_tsne.R $(ROUT)/000_cellcnn_plot_tsne.Rout













#
