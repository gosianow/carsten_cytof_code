
MAKEARGS :=

R := R CMD BATCH --no-restore --no-save

RWD_MAIN := ../carsten_cytof/PD1_project
RCODE := .
RWD_CELLCNN := ../PD1_CellCnn_Lukas_June2017

METADATA := $(RWD_MAIN)/CK_metadata
PANELS := $(RWD_MAIN)/CK_panels


## Define the default rule (makefiles are usually written so that the first target is for compiling the entire program)
.PHONY: all

all: cnn_23_01_panel1_base_goal cnn_23CD4_01CD4_panel1_CD4_Tcells_base_goal cnn_23CD8_01CD8_panel1_CD8_Tcells_base_goal \
	cnn_23_02_panel2_base_goal cnn_23CD4_02CD4_panel2_CD4_Tcells_base_goal cnn_23CD8_02CD8_panel2_CD8_Tcells_base_goal \
	cnn_23all_03v2_panel3v3_base_goal


### Make sure no intermediate files are deleted
.SECONDARY:


###############################################################################################################
# Analysis of CK_2016-06-23_01 data
# Use Analysis block 1
###############################################################################################################


data_dir := CK_2016-06-23_01
file_metadata := $(METADATA)/metadata_23_01.xlsx

data := 23
panel := 01
pca := pca1
merging := merging6

cellcnn_files_data := panel1_base
cellcnn_files_type := combined data23
day := base


.PHONY: cnn_$(data)_$(panel)_$(cellcnn_files_data)_goal $(foreach i,$(cellcnn_files_type),cnn_$(data)_$(panel)_$(cellcnn_files_data)_$(i)_goal)
cnn_$(data)_$(panel)_$(cellcnn_files_data)_goal: $(foreach i,$(cellcnn_files_type),cnn_$(data)_$(panel)_$(cellcnn_files_data)_$(i)_goal)

define cnn_rule
cnn_$(data)_$(panel)_$(cellcnn_files_data)_$(1)_goal:
	echo "\n>> make"
	make $(MAKEARGS) -f 000_cellcnn_plot_tsne.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" RWD_CELLCNN="$(RWD_CELLCNN)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" merging="$(merging)" cellcnn_files_data="$(cellcnn_files_data)" cellcnn_files_type="$(1)" day="$(day)"
endef
$(foreach i,$(cellcnn_files_type),$(eval $(call cnn_rule,$(i))))


# -----------------------------------------------------------------------------
# CD4

data_dir := CK_2016-06-23_01_CD4_mergingNEW2
file_metadata := $(METADATA)/metadata_23_01.xlsx

data := 23CD4
panel := 01CD4
pca := pca1
merging := merging5

cellcnn_files_data := panel1_CD4_Tcells_base
cellcnn_files_type := combined data23
day := base


.PHONY: cnn_$(data)_$(panel)_$(cellcnn_files_data)_goal $(foreach i,$(cellcnn_files_type),cnn_$(data)_$(panel)_$(cellcnn_files_data)_$(i)_goal)
cnn_$(data)_$(panel)_$(cellcnn_files_data)_goal: $(foreach i,$(cellcnn_files_type),cnn_$(data)_$(panel)_$(cellcnn_files_data)_$(i)_goal)

define cnn_rule
cnn_$(data)_$(panel)_$(cellcnn_files_data)_$(1)_goal:
	echo "\n>> make"
	make $(MAKEARGS) -f 000_cellcnn_plot_tsne.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" RWD_CELLCNN="$(RWD_CELLCNN)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" merging="$(merging)" cellcnn_files_data="$(cellcnn_files_data)" cellcnn_files_type="$(1)" day="$(day)"
endef
$(foreach i,$(cellcnn_files_type),$(eval $(call cnn_rule,$(i))))


# -----------------------------------------------------------------------------
# CD8

data_dir := CK_2016-06-23_01_CD8_mergingNEW2
file_metadata := $(METADATA)/metadata_23_01.xlsx

data := 23CD8
panel := 01CD8
pca := pca1
merging := merging6

cellcnn_files_data := panel1_CD8_Tcells_base
cellcnn_files_type := combined data23
day := base


.PHONY: cnn_$(data)_$(panel)_$(cellcnn_files_data)_goal $(foreach i,$(cellcnn_files_type),cnn_$(data)_$(panel)_$(cellcnn_files_data)_$(i)_goal)
cnn_$(data)_$(panel)_$(cellcnn_files_data)_goal: $(foreach i,$(cellcnn_files_type),cnn_$(data)_$(panel)_$(cellcnn_files_data)_$(i)_goal)

define cnn_rule
cnn_$(data)_$(panel)_$(cellcnn_files_data)_$(1)_goal:
	echo "\n>> make"
	make $(MAKEARGS) -f 000_cellcnn_plot_tsne.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" RWD_CELLCNN="$(RWD_CELLCNN)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" merging="$(merging)" cellcnn_files_data="$(cellcnn_files_data)" cellcnn_files_type="$(1)" day="$(day)"
endef
$(foreach i,$(cellcnn_files_type),$(eval $(call cnn_rule,$(i))))



###############################################################################################################
# Analysis of CK_2016-06-23_02 data
# Use Analysis block 1
###############################################################################################################


data_dir := CK_2016-06-23_02
file_metadata := $(METADATA)/metadata_23_02.xlsx

data := 23
panel := 02
pca := pca1
merging := merging2

cellcnn_files_data := panel2_base
cellcnn_files_type := combined data23
day := base


.PHONY: cnn_$(data)_$(panel)_$(cellcnn_files_data)_goal $(foreach i,$(cellcnn_files_type),cnn_$(data)_$(panel)_$(cellcnn_files_data)_$(i)_goal)
cnn_$(data)_$(panel)_$(cellcnn_files_data)_goal: $(foreach i,$(cellcnn_files_type),cnn_$(data)_$(panel)_$(cellcnn_files_data)_$(i)_goal)

define cnn_rule
cnn_$(data)_$(panel)_$(cellcnn_files_data)_$(1)_goal:
	echo "\n>> make"
	make $(MAKEARGS) -f 000_cellcnn_plot_tsne.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" RWD_CELLCNN="$(RWD_CELLCNN)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" merging="$(merging)" cellcnn_files_data="$(cellcnn_files_data)" cellcnn_files_type="$(1)" day="$(day)"
endef
$(foreach i,$(cellcnn_files_type),$(eval $(call cnn_rule,$(i))))


# -----------------------------------------------------------------------------
# CD4

data_dir := CK_2016-06-23_02_CD4_merging2
file_metadata := $(METADATA)/metadata_23_02.xlsx

data := 23CD4
panel := 02CD4
pca := pca1
merging := merging2

cellcnn_files_data := panel2_CD4_Tcells_base
cellcnn_files_type := combined data23
day := base


.PHONY: cnn_$(data)_$(panel)_$(cellcnn_files_data)_goal $(foreach i,$(cellcnn_files_type),cnn_$(data)_$(panel)_$(cellcnn_files_data)_$(i)_goal)
cnn_$(data)_$(panel)_$(cellcnn_files_data)_goal: $(foreach i,$(cellcnn_files_type),cnn_$(data)_$(panel)_$(cellcnn_files_data)_$(i)_goal)

define cnn_rule
cnn_$(data)_$(panel)_$(cellcnn_files_data)_$(1)_goal:
	echo "\n>> make"
	make $(MAKEARGS) -f 000_cellcnn_plot_tsne.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" RWD_CELLCNN="$(RWD_CELLCNN)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" merging="$(merging)" cellcnn_files_data="$(cellcnn_files_data)" cellcnn_files_type="$(1)" day="$(day)"
endef
$(foreach i,$(cellcnn_files_type),$(eval $(call cnn_rule,$(i))))


# -----------------------------------------------------------------------------
# CD8

data_dir := CK_2016-06-23_02_CD8_merging2
file_metadata := $(METADATA)/metadata_23_02.xlsx

data := 23CD8
panel := 02CD8
pca := pca1
merging := merging2

cellcnn_files_data := panel2_CD8_Tcells_base
cellcnn_files_type := combined data23
day := base


.PHONY: cnn_$(data)_$(panel)_$(cellcnn_files_data)_goal $(foreach i,$(cellcnn_files_type),cnn_$(data)_$(panel)_$(cellcnn_files_data)_$(i)_goal)
cnn_$(data)_$(panel)_$(cellcnn_files_data)_goal: $(foreach i,$(cellcnn_files_type),cnn_$(data)_$(panel)_$(cellcnn_files_data)_$(i)_goal)

define cnn_rule
cnn_$(data)_$(panel)_$(cellcnn_files_data)_$(1)_goal:
	echo "\n>> make"
	make $(MAKEARGS) -f 000_cellcnn_plot_tsne.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" RWD_CELLCNN="$(RWD_CELLCNN)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" merging="$(merging)" cellcnn_files_data="$(cellcnn_files_data)" cellcnn_files_type="$(1)" day="$(day)"
endef
$(foreach i,$(cellcnn_files_type),$(eval $(call cnn_rule,$(i))))


###############################################################################################################
# Analysis of CK_2016-06-23_03all data
# Use Analysis block 1
###############################################################################################################

data_dir := CK_2016-06-23_03all
file_metadata := $(METADATA)/metadata_23_03all.xlsx

data := 23all
panel := 03v2
pca := pca1
merging := merging5

cellcnn_files_data := panel3v3_base
cellcnn_files_type := combined data23
day := base


.PHONY: cnn_$(data)_$(panel)_$(cellcnn_files_data)_goal $(foreach i,$(cellcnn_files_type),cnn_$(data)_$(panel)_$(cellcnn_files_data)_$(i)_goal)
cnn_$(data)_$(panel)_$(cellcnn_files_data)_goal: $(foreach i,$(cellcnn_files_type),cnn_$(data)_$(panel)_$(cellcnn_files_data)_$(i)_goal)

define cnn_rule
cnn_$(data)_$(panel)_$(cellcnn_files_data)_$(1)_goal:
	echo "\n>> make"
	make $(MAKEARGS) -f 000_cellcnn_plot_tsne.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" RWD_CELLCNN="$(RWD_CELLCNN)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" merging="$(merging)" cellcnn_files_data="$(cellcnn_files_data)" cellcnn_files_type="$(1)" day="$(day)"
endef
$(foreach i,$(cellcnn_files_type),$(eval $(call cnn_rule,$(i))))




































#
