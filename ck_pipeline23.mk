
MAKEARGS :=

R := R CMD BATCH --no-restore --no-save

RWD_MAIN := ../carsten_cytof/PD1_project
RCODE := .

METADATA := $(RWD_MAIN)/CK_metadata
PANELS := $(RWD_MAIN)/CK_panels


## Define the default rule (makefiles are usually written so that the first target is for compiling the entire program)
.PHONY: all analysis_panel1 analysis_panel2 analysis_panel3 \
	Analysis_block_1_main_23_01_goal Analysis_block_2_cluster_merging_23_01_cl20_goal Analysis_block_3_cluster_extracting_23_01_goal \
		Analysis_block_1_main_23CD4_01CD4_goal Analysis_block_1_main_23CD8_01CD8_goal \
		Analysis_block_2_cluster_merging_23CD4_01CD4_cl20_goal Analysis_block_2_cluster_merging_23CD8_01CD8_cl20_goal \
		Analysis_block_2_cluster_merging_23CD4_01CD4_cl5_goal Analysis_block_2_cluster_merging_23CD8_01CD8_cl5_goal \
		Analysis_block_1_main_23_02_goal Analysis_block_2_cluster_merging_23_02_cl20_goal Analysis_block_3_cluster_extracting_23_02_goal \
		Analysis_block_1_main_23CD4_02CD4_goal Analysis_block_1_main_23CD8_02CD8_goal \
		Analysis_block_2_cluster_merging_23CD4_02CD4_cl6_goal Analysis_block_2_cluster_merging_23CD8_02CD8_cl7_goal \
		Analysis_block_3_cytokine_cell_extracting_23CD4_02CD4_goal \
		Analysis_block_4_cytokine_main_23CD4TmemCD69_02CD4_goal \
		Analysis_block_3_cytokine_cell_extracting_23CD8_02CD8_goal \
		Analysis_block_4_cytokine_main_23CD8TmemCD69_02CD8_goal \
		Analysis_block_1_main_23_02v2_goal Analysis_block_2_cluster_merging_23_02v2_cl20_goal Analysis_block_3_cluster_extracting_23_02v2_goal \
		Analysis_block_1_main_23CD4_02CD4v2_goal Analysis_block_1_main_23CD8_02CD8v2_goal \
		Analysis_block_2_cluster_merging_23CD4_02CD4v2_cl20_goal Analysis_block_2_cluster_merging_23CD8_02CD8v2_cl20_goal \
		Analysis_block_3_cytokine_cell_extracting_23CD4_02CD4v2_goal \
		Analysis_block_4_cytokine_main_02CD4v2_goal \
		Analysis_block_3_cytokine_cell_extracting_23CD8_02CD8v2_goal \
		Analysis_block_4_cytokine_main_02CD8v2_goal \
		Analysis_block_1_main_23_03_goal Analysis_block_2_cluster_merging_23_03_cl20_goal \
		Analysis_block_1_main_23all_03v2_goal Analysis_block_2_cluster_merging_23all_03v2_cl20_goal

all: analysis_panel1 analysis_panel2 analysis_panel2v2 analysis_panel3

analysis_panel1: Analysis_block_1_main_23_01_goal Analysis_block_2_cluster_merging_23_01_cl20_goal Analysis_block_3_cluster_extracting_23_01_goal \
	Analysis_block_1_main_23CD4_01CD4_goal Analysis_block_1_main_23CD8_01CD8_goal \
	Analysis_block_2_cluster_merging_23CD4_01CD4_cl20_goal Analysis_block_2_cluster_merging_23CD8_01CD8_cl20_goal \
	Analysis_block_2_cluster_merging_23CD4_01CD4_cl5_goal Analysis_block_2_cluster_merging_23CD8_01CD8_cl5_goal

analysis_panel2: Analysis_block_1_main_23_02_goal Analysis_block_2_cluster_merging_23_02_cl20_goal Analysis_block_3_cluster_extracting_23_02_goal \
	Analysis_block_1_main_23CD4_02CD4_goal Analysis_block_1_main_23CD8_02CD8_goal \
	Analysis_block_2_cluster_merging_23CD4_02CD4_cl6_goal Analysis_block_2_cluster_merging_23CD8_02CD8_cl7_goal \
	Analysis_block_3_cytokine_cell_extracting_23CD4_02CD4_goal \
	Analysis_block_4_cytokine_main_23CD4TmemCD69_02CD4_goal \
	Analysis_block_3_cytokine_cell_extracting_23CD8_02CD8_goal \
	Analysis_block_4_cytokine_main_23CD8TmemCD69_02CD8_goal

analysis_panel2v2: Analysis_block_1_main_23_02v2_goal Analysis_block_2_cluster_merging_23_02v2_cl20_goal Analysis_block_3_cluster_extracting_23_02v2_goal \
	Analysis_block_1_main_23CD4_02CD4v2_goal Analysis_block_1_main_23CD8_02CD8v2_goal \
	Analysis_block_2_cluster_merging_23CD4_02CD4v2_cl20_goal Analysis_block_2_cluster_merging_23CD8_02CD8v2_cl20_goal \
	Analysis_block_3_cytokine_cell_extracting_23CD4_02CD4v2_goal \
	Analysis_block_4_cytokine_main_02CD4v2_goal \
	Analysis_block_3_cytokine_cell_extracting_23CD8_02CD8v2_goal \
	Analysis_block_4_cytokine_main_02CD8v2_goal \
	Analysis_block_3_cytokine_cell_extracting_23_02v2_goal \
	Analysis_block_4_cytokine_main_02Tregsv2_goal

analysis_panel2: Analysis_block_1_main_23_03_goal Analysis_block_2_cluster_merging_23_03_cl20_goal \
	Analysis_block_1_main_23all_03v2_goal Analysis_block_2_cluster_merging_23all_03v2_cl20_goal



### Make sure no intermediate files are deleted
.SECONDARY:


###############################################################################################################
# Analysis of CK_2016-06-23_01 data
# Use Analysis block 1
###############################################################################################################


data_dir := CK_2016-06-23_01
file_panel := $(PANELS)/panel1.xlsx
file_metadata := $(METADATA)/metadata_23_01.xlsx

data := 23
panel := 01
pca := pca1
nmetaclusts := 20

rand_seed_consensus := 123


.PHONY: Analysis_block_1_main_$(data)_$(panel)_goal
Analysis_block_1_main_$(data)_$(panel)_goal: $(foreach i,$(nmetaclusts),Analysis_block_1_main_$(data)_$(panel)_goal$(i))

define Analysis_block_1_main_rule
Analysis_block_1_main_$(data)_$(panel)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_1_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_panel="$(file_panel)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" rand_seed_consensus="$(rand_seed_consensus)" nmetaclusts="$(1)"
endef
$(foreach i,$(nmetaclusts),$(eval $(call Analysis_block_1_main_rule,$(i))))


# --------------------------------------------------
# Analysis of CK_2016-06-23_01 cluster_merging
# Use Analysis block 2
# --------------------------------------------------


nmetaclusts := 20
merging := mergingNEW2 merging6 merging7


.PHONY: Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal: $(foreach i,$(merging),Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(i))


define Analysis_block_2_cluster_merging_rule
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_2_cluster_merging.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" nmetaclusts="$(nmetaclusts)" merging="$(1)"
endef
$(foreach i,$(merging),$(eval $(call Analysis_block_2_cluster_merging_rule,$(i))))



# --------------------------------------------------
# CK_2016-06-23_01 - CD4 and CD8 cluster extracting from mergingNEW2
# Use Analysis block 3
# --------------------------------------------------

merging := mergingNEW2
extract_cluster := CD4 CD8


.PHONY: Analysis_block_3_cluster_extracting_$(data)_$(panel)_goal
Analysis_block_3_cluster_extracting_$(data)_$(panel)_goal: $(foreach i,$(extract_cluster),Analysis_block_3_cluster_extracting_$(data)_$(panel)_goal$(i))


define Analysis_block_3_cluster_extracting_rule
Analysis_block_3_cluster_extracting_$(data)_$(panel)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_3_cluster_extracting.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" merging="$(merging)" extract_cluster="$(1)"
endef
$(foreach i,$(extract_cluster),$(eval $(call Analysis_block_3_cluster_extracting_rule,$(i))))



# ----------------------------------------------------------------------------------------------------
# Analysis of CK_2016-06-23_01_CD4_mergingNEW2 using panel1CD4.xlsx
# Use Analysis block 1
# ----------------------------------------------------------------------------------------------------


data_dir := CK_2016-06-23_01_CD4_mergingNEW2
file_panel := $(PANELS)/panel1CD4.xlsx
file_metadata := $(METADATA)/metadata_23_01.xlsx

data := 23CD4
panel := 01CD4
pca := pca1
nmetaclusts := 20 5

rand_seed_consensus := 1234


.PHONY: Analysis_block_1_main_$(data)_$(panel)_goal
Analysis_block_1_main_$(data)_$(panel)_goal: $(foreach i,$(nmetaclusts),Analysis_block_1_main_$(data)_$(panel)_goal$(i))

define Analysis_block_1_main_rule
Analysis_block_1_main_$(data)_$(panel)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_1_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_panel="$(file_panel)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" rand_seed_consensus="$(rand_seed_consensus)" nmetaclusts="$(1)"
endef
$(foreach i,$(nmetaclusts),$(eval $(call Analysis_block_1_main_rule,$(i))))


# --------------------------------------------------
# Analysis of CK_2016-06-23_01_CD4_mergingNEW2 cluster merging
# Use Analysis block 2
# --------------------------------------------------

nmetaclusts := 20
merging := merging6

.PHONY: Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal: $(foreach i,$(merging),Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(i))

define Analysis_block_2_cluster_merging_rule
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_2_cluster_merging.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" nmetaclusts="$(nmetaclusts)" merging="$(1)"
endef
$(foreach i,$(merging),$(eval $(call Analysis_block_2_cluster_merging_rule,$(i))))


nmetaclusts := 5
merging := merging5

.PHONY: Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal: $(foreach i,$(merging),Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(i))

define Analysis_block_2_cluster_merging_rule
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_2_cluster_merging.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" nmetaclusts="$(nmetaclusts)" merging="$(1)"
endef
$(foreach i,$(merging),$(eval $(call Analysis_block_2_cluster_merging_rule,$(i))))


# ----------------------------------------------------------------------------------------------------
# Analysis of CK_2016-06-23_01_CD8_mergingNEW2 using panel1CD8.xlsx
# Use Analysis block 1
# ----------------------------------------------------------------------------------------------------


data_dir := CK_2016-06-23_01_CD8_mergingNEW2
file_panel := $(PANELS)/panel1CD8.xlsx
file_metadata := $(METADATA)/metadata_23_01.xlsx

data := 23CD8
panel := 01CD8
pca := pca1
nmetaclusts := 20 5

rand_seed_consensus := 1234


.PHONY: Analysis_block_1_main_$(data)_$(panel)_goal
Analysis_block_1_main_$(data)_$(panel)_goal: $(foreach i,$(nmetaclusts),Analysis_block_1_main_$(data)_$(panel)_goal$(i))

define Analysis_block_1_main_rule
Analysis_block_1_main_$(data)_$(panel)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_1_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_panel="$(file_panel)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" rand_seed_consensus="$(rand_seed_consensus)" nmetaclusts="$(1)"
endef
$(foreach i,$(nmetaclusts),$(eval $(call Analysis_block_1_main_rule,$(i))))


# --------------------------------------------------
# Analysis of CK_2016-06-23_01_CD8_mergingNEW2 cluster merging
# Use Analysis block 2
# --------------------------------------------------

nmetaclusts := 20
merging := merging7

.PHONY: Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal: $(foreach i,$(merging),Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(i))

define Analysis_block_2_cluster_merging_rule
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_2_cluster_merging.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" nmetaclusts="$(nmetaclusts)" merging="$(1)"
endef
$(foreach i,$(merging),$(eval $(call Analysis_block_2_cluster_merging_rule,$(i))))


nmetaclusts := 5
merging := merging6

.PHONY: Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal: $(foreach i,$(merging),Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(i))

define Analysis_block_2_cluster_merging_rule
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_2_cluster_merging.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" nmetaclusts="$(nmetaclusts)" merging="$(1)"
endef
$(foreach i,$(merging),$(eval $(call Analysis_block_2_cluster_merging_rule,$(i))))



###############################################################################################################
# Analysis of CK_2016-06-23_02 data using panel 02v2
# Use Analysis block 1
###############################################################################################################


data_dir := CK_2016-06-23_02
file_panel := $(PANELS)/panel2_v2.xlsx
file_metadata := $(METADATA)/metadata_23_02.xlsx

data := 23
panel := 02v2
pca := pca0
nmetaclusts := 20

rand_seed_consensus := 1234


.PHONY: Analysis_block_1_main_$(data)_$(panel)_goal
Analysis_block_1_main_$(data)_$(panel)_goal: $(foreach i,$(nmetaclusts),Analysis_block_1_main_$(data)_$(panel)_goal$(i))

define Analysis_block_1_main_rule
Analysis_block_1_main_$(data)_$(panel)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_1_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_panel="$(file_panel)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" rand_seed_consensus="$(rand_seed_consensus)" nmetaclusts="$(1)"
endef
$(foreach i,$(nmetaclusts),$(eval $(call Analysis_block_1_main_rule,$(i))))


# --------------------------------------------------
# Analysis of CK_2016-06-23_02 cluster merging
# Use Analysis block 2
# --------------------------------------------------

nmetaclusts := 20
merging := merging3

.PHONY: Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal: $(foreach i,$(merging),Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(i))

define Analysis_block_2_cluster_merging_rule
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_2_cluster_merging.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" nmetaclusts="$(nmetaclusts)" merging="$(1)"
endef
$(foreach i,$(merging),$(eval $(call Analysis_block_2_cluster_merging_rule,$(i))))


# --------------------------------------------------
# CK_2016-06-23_02 - CD4 and CD8 cluster extracting from merging
# Use Analysis block 3
# --------------------------------------------------

merging := merging3
extract_cluster := CD4 CD8

.PHONY: Analysis_block_3_cluster_extracting_$(data)_$(panel)_goal
Analysis_block_3_cluster_extracting_$(data)_$(panel)_goal: $(foreach i,$(extract_cluster),Analysis_block_3_cluster_extracting_$(data)_$(panel)_goal$(i))

define Analysis_block_3_cluster_extracting_rule
Analysis_block_3_cluster_extracting_$(data)_$(panel)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_3_cluster_extracting.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" merging="$(merging)" extract_cluster="$(1)"
endef
$(foreach i,$(extract_cluster),$(eval $(call Analysis_block_3_cluster_extracting_rule,$(i))))



# ----------------------------------------------------------------------------------------------------
# Analysis of CK_2016-06-23_02_CD4_merging3 using panel2CD4_v2.xlsx
# Use Analysis block 1
# ----------------------------------------------------------------------------------------------------


data_dir := CK_2016-06-23_02_CD4_merging3
file_panel := $(PANELS)/panel2CD4_v2.xlsx
file_metadata := $(METADATA)/metadata_23_02.xlsx

data := 23CD4
panel := 02CD4v2
pca := pca0
nmetaclusts := 20

rand_seed_consensus := 1234


.PHONY: Analysis_block_1_main_$(data)_$(panel)_goal
Analysis_block_1_main_$(data)_$(panel)_goal: $(foreach i,$(nmetaclusts),Analysis_block_1_main_$(data)_$(panel)_goal$(i))

define Analysis_block_1_main_rule
Analysis_block_1_main_$(data)_$(panel)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_1_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_panel="$(file_panel)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" rand_seed_consensus="$(rand_seed_consensus)" nmetaclusts="$(1)"
endef
$(foreach i,$(nmetaclusts),$(eval $(call Analysis_block_1_main_rule,$(i))))


# --------------------------------------------------
# Analysis of CK_2016-06-23_02_CD4_merging2 cluster merging
# Use Analysis block 2
# --------------------------------------------------

nmetaclusts := 20
merging := merging1

.PHONY: Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal: $(foreach i,$(merging),Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(i))

define Analysis_block_2_cluster_merging_rule
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_2_cluster_merging.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" nmetaclusts="$(nmetaclusts)" merging="$(1)"
endef
$(foreach i,$(merging),$(eval $(call Analysis_block_2_cluster_merging_rule,$(i))))


####################################################
# CYTOKINE ANALYSIS


# --------------------------------------------------
# Analysis of CK_2016-06-23_02_CD4_merging3 cytokine extracting
# Use Analysis block 3
# --------------------------------------------------

file_extract_marker := $(PANELS)/panel2CD4_23_cytokines_all.xlsx
extract_marker := all

merging := merging1
extract_cluster := all naive CM


.PHONY: Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_goal
Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_goal: $(foreach i,$(extract_cluster),Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_goal$(i))

define Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_rule
Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_3_cytokine_cell_extracting.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" merging="$(merging)" extract_cluster="$(1)" file_extract_marker="$(file_extract_marker)" extract_marker="$(extract_marker)"
endef
$(foreach i,$(extract_cluster),$(eval $(call Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_rule,$(i))))


# --------------------------------------------------
# Analysis of CK_2016-06-23_02_CD4_merging3_all_merging1_all cytokine analysis
# Use Analysis block 4
# --------------------------------------------------

data_dir := CK_2016-06-23_02_CD4_merging3
file_metadata := $(METADATA)/metadata_23_02.xlsx
file_panel_cytokines := $(PANELS)/panel2CD4_23_cytokines.xlsx

data := 23CD4
cytokines := 02CD4v2

som_dim := 7
nmetaclusts := 49


.PHONY: Analysis_block_4_cytokine_main_$(cytokines)_goal
Analysis_block_4_cytokine_main_$(cytokines)_goal: $(foreach i,$(extract_cluster),Analysis_block_4_cytokine_main_$(cytokines)_goal$(i))

define Analysis_block_4_cytokine_main_$(cytokines)_rule
Analysis_block_4_cytokine_main_$(cytokines)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_4_cytokine_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)_$(1)_$(merging)_$(extract_marker)" file_metadata="$(file_metadata)" file_panel_cytokines="$(file_panel_cytokines)" data="$(data)$(1)$(extract_marker)" cytokines="$(cytokines)" som_dim="$(som_dim)" nmetaclusts="$(nmetaclusts)"
endef
$(foreach i,$(extract_cluster),$(eval $(call Analysis_block_4_cytokine_main_$(cytokines)_rule,$(i))))



####################################################


# ----------------------------------------------------------------------------------------------------
# Analysis of CK_2016-06-23_02_CD8_merging3 using panel2CD8_v2.xlsx
# Use Analysis block 1
# ----------------------------------------------------------------------------------------------------


data_dir := CK_2016-06-23_02_CD8_merging3
file_panel := $(PANELS)/panel2CD8_v2.xlsx
file_metadata := $(METADATA)/metadata_23_02.xlsx

data := 23CD8
panel := 02CD8v2
pca := pca0
nmetaclusts := 20

rand_seed_consensus := 1234


.PHONY: Analysis_block_1_main_$(data)_$(panel)_goal
Analysis_block_1_main_$(data)_$(panel)_goal: $(foreach i,$(nmetaclusts),Analysis_block_1_main_$(data)_$(panel)_goal$(i))

define Analysis_block_1_main_rule
Analysis_block_1_main_$(data)_$(panel)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_1_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_panel="$(file_panel)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" rand_seed_consensus="$(rand_seed_consensus)" nmetaclusts="$(1)"
endef
$(foreach i,$(nmetaclusts),$(eval $(call Analysis_block_1_main_rule,$(i))))


# --------------------------------------------------
# Analysis of CK_2016-06-23_02_CD8_merging3 cluster merging
# Use Analysis block 2
# --------------------------------------------------

nmetaclusts := 20
merging := merging1

.PHONY: Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal: $(foreach i,$(merging),Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(i))

define Analysis_block_2_cluster_merging_rule
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_2_cluster_merging.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" nmetaclusts="$(nmetaclusts)" merging="$(1)"
endef
$(foreach i,$(merging),$(eval $(call Analysis_block_2_cluster_merging_rule,$(i))))


####################################################
# CYTOKINE ANALYSIS


# --------------------------------------------------
# Analysis of CK_2016-06-23_02_CD8_merging3 cytokine extracting
# Use Analysis block 3
# --------------------------------------------------

file_extract_marker := $(PANELS)/panel2CD8_23_cytokines_all.xlsx
extract_marker := all

merging := merging1
extract_cluster := all naive CM

.PHONY: Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_goal
Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_goal: $(foreach i,$(extract_cluster),Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_goal$(i))

define Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_rule
Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_3_cytokine_cell_extracting.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" merging="$(merging)" extract_cluster="$(1)" file_extract_marker="$(file_extract_marker)" extract_marker="$(extract_marker)"
endef
$(foreach i,$(extract_cluster),$(eval $(call Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_rule,$(i))))


# --------------------------------------------------
# Analysis of CK_2016-06-23_02_CD8_merging3_all_merging1_all cytokine analysis
# Use Analysis block 4
# --------------------------------------------------

data_dir := CK_2016-06-23_02_CD8_merging3
file_metadata := $(METADATA)/metadata_23_02.xlsx
file_panel_cytokines := $(PANELS)/panel2CD8_23_cytokines.xlsx

data := 23CD8
cytokines := 02CD8v2

som_dim := 7
nmetaclusts := 49


.PHONY: Analysis_block_4_cytokine_main_$(cytokines)_goal
Analysis_block_4_cytokine_main_$(cytokines)_goal: $(foreach i,$(extract_cluster),Analysis_block_4_cytokine_main_$(cytokines)_goal$(i))

define Analysis_block_4_cytokine_main_$(cytokines)_rule
Analysis_block_4_cytokine_main_$(cytokines)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_4_cytokine_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)_$(1)_$(merging)_$(extract_marker)" file_metadata="$(file_metadata)" file_panel_cytokines="$(file_panel_cytokines)" data="$(data)$(1)$(extract_marker)" cytokines="$(cytokines)" som_dim="$(som_dim)" nmetaclusts="$(nmetaclusts)"
endef
$(foreach i,$(extract_cluster),$(eval $(call Analysis_block_4_cytokine_main_$(cytokines)_rule,$(i))))



####################################################


# ----------------------------------------------------------------------------------------------------
# Analysis of CK_2016-06-23_02 Tregs using panel2CD4_v2.xlsx
# Use Analysis block 1
# ----------------------------------------------------------------------------------------------------

data_dir := CK_2016-06-23_02
file_panel := $(PANELS)/panel2CD4.xlsx
file_metadata := $(METADATA)/metadata_23_02.xlsx

data := 23
panel := 02v2
pca := pca0

####################################################
# CYTOKINE ANALYSIS


# --------------------------------------------------
# Analysis of CK_2016-06-23_02 Tregs cytokine extracting
# Use Analysis block 3
# --------------------------------------------------

file_extract_marker := $(PANELS)/panel2CD4_23_cytokines_all.xlsx
extract_marker := all

merging := merging3
extract_cluster := Tregs


.PHONY: Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_goal
Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_goal: $(foreach i,$(extract_cluster),Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_goal$(i))

define Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_rule
Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_3_cytokine_cell_extracting.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" merging="$(merging)" extract_cluster="$(1)" file_extract_marker="$(file_extract_marker)" extract_marker="$(extract_marker)"
endef
$(foreach i,$(extract_cluster),$(eval $(call Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_rule,$(i))))


# --------------------------------------------------
# Analysis of CK_2016-06-23_02_Tregs_merging3_all cytokine analysis
# Use Analysis block 4
# --------------------------------------------------

data_dir := CK_2016-06-23_02
file_metadata := $(METADATA)/metadata_23_02.xlsx
file_panel_cytokines := $(PANELS)/panel2CD4_23_cytokines.xlsx

data := 23
cytokines := 02Tregsv2

som_dim := 2
nmetaclusts := 4


.PHONY: Analysis_block_4_cytokine_main_$(cytokines)_goal
Analysis_block_4_cytokine_main_$(cytokines)_goal: $(foreach i,$(extract_cluster),Analysis_block_4_cytokine_main_$(cytokines)_goal$(i))

define Analysis_block_4_cytokine_main_$(cytokines)_rule
Analysis_block_4_cytokine_main_$(cytokines)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_4_cytokine_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)_$(1)_$(merging)_$(extract_marker)" file_metadata="$(file_metadata)" file_panel_cytokines="$(file_panel_cytokines)" data="$(data)$(1)$(extract_marker)" cytokines="$(cytokines)" som_dim="$(som_dim)" nmetaclusts="$(nmetaclusts)"
endef
$(foreach i,$(extract_cluster),$(eval $(call Analysis_block_4_cytokine_main_$(cytokines)_rule,$(i))))



####################################################



###############################################################################################################
# Analysis of CK_2016-06-23_02 data
# Use Analysis block 1
###############################################################################################################


data_dir := CK_2016-06-23_02
file_panel := $(PANELS)/panel2.xlsx
file_metadata := $(METADATA)/metadata_23_02.xlsx

data := 23
panel := 02
pca := pca1
nmetaclusts := 20

rand_seed_consensus := 1234


.PHONY: Analysis_block_1_main_$(data)_$(panel)_goal
Analysis_block_1_main_$(data)_$(panel)_goal: $(foreach i,$(nmetaclusts),Analysis_block_1_main_$(data)_$(panel)_goal$(i))

define Analysis_block_1_main_rule
Analysis_block_1_main_$(data)_$(panel)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_1_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_panel="$(file_panel)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" rand_seed_consensus="$(rand_seed_consensus)" nmetaclusts="$(1)"
endef
$(foreach i,$(nmetaclusts),$(eval $(call Analysis_block_1_main_rule,$(i))))


# --------------------------------------------------
# Analysis of CK_2016-06-23_02 cluster merging
# Use Analysis block 2
# --------------------------------------------------

nmetaclusts := 20
merging := merging2

.PHONY: Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal: $(foreach i,$(merging),Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(i))

define Analysis_block_2_cluster_merging_rule
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_2_cluster_merging.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" nmetaclusts="$(nmetaclusts)" merging="$(1)"
endef
$(foreach i,$(merging),$(eval $(call Analysis_block_2_cluster_merging_rule,$(i))))


# --------------------------------------------------
# CK_2016-06-23_02 - CD4 and CD8 cluster extracting from merging
# Use Analysis block 3
# --------------------------------------------------

merging := merging2
extract_cluster := CD4 CD8

.PHONY: Analysis_block_3_cluster_extracting_$(data)_$(panel)_goal
Analysis_block_3_cluster_extracting_$(data)_$(panel)_goal: $(foreach i,$(extract_cluster),Analysis_block_3_cluster_extracting_$(data)_$(panel)_goal$(i))

define Analysis_block_3_cluster_extracting_rule
Analysis_block_3_cluster_extracting_$(data)_$(panel)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_3_cluster_extracting.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" merging="$(merging)" extract_cluster="$(1)"
endef
$(foreach i,$(extract_cluster),$(eval $(call Analysis_block_3_cluster_extracting_rule,$(i))))



# ----------------------------------------------------------------------------------------------------
# Analysis of CK_2016-06-23_02_CD4_merging2 using panel2CD4.xlsx
# Use Analysis block 1
# ----------------------------------------------------------------------------------------------------


data_dir := CK_2016-06-23_02_CD4_merging2
file_panel := $(PANELS)/panel2CD4.xlsx
file_metadata := $(METADATA)/metadata_23_02.xlsx

data := 23CD4
panel := 02CD4
pca := pca1
nmetaclusts := 20 6

rand_seed_consensus := 1234


.PHONY: Analysis_block_1_main_$(data)_$(panel)_goal
Analysis_block_1_main_$(data)_$(panel)_goal: $(foreach i,$(nmetaclusts),Analysis_block_1_main_$(data)_$(panel)_goal$(i))

define Analysis_block_1_main_rule
Analysis_block_1_main_$(data)_$(panel)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_1_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_panel="$(file_panel)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" rand_seed_consensus="$(rand_seed_consensus)" nmetaclusts="$(1)"
endef
$(foreach i,$(nmetaclusts),$(eval $(call Analysis_block_1_main_rule,$(i))))


# --------------------------------------------------
# Analysis of CK_2016-06-23_02_CD4_merging2 cluster merging
# Use Analysis block 2
# --------------------------------------------------

nmetaclusts := 6
merging := merging2

.PHONY: Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal: $(foreach i,$(merging),Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(i))

define Analysis_block_2_cluster_merging_rule
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_2_cluster_merging.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" nmetaclusts="$(nmetaclusts)" merging="$(1)"
endef
$(foreach i,$(merging),$(eval $(call Analysis_block_2_cluster_merging_rule,$(i))))


####################################################
# CYTOKINE ANALYSIS


# --------------------------------------------------
# Analysis of CK_2016-06-23_02_CD4_merging2 cytokine extracting
# Use Analysis block 3
# --------------------------------------------------

file_extract_marker := $(PANELS)/panel2CD4_23_cytokines_CD69.xlsx
extract_marker := CD69

merging := merging2
extract_cluster := Tmem

.PHONY: Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_goal

define Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_rule
Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_goal:
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_3_cytokine_cell_extracting.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" merging="$(merging)" extract_cluster="$(extract_cluster)" file_extract_marker="$(file_extract_marker)" extract_marker="$(extract_marker)"
endef

$(eval $(call Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_rule,))


# --------------------------------------------------
# Analysis of CK_2016-06-23_02_CD4_merging2_Tmem_merging2_CD69 cytokine analysis
# Use Analysis block 4
# --------------------------------------------------

data_dir := CK_2016-06-23_02_CD4_merging2_Tmem_merging2_CD69
file_metadata := $(METADATA)/metadata_23_02.xlsx
file_panel_cytokines := $(PANELS)/panel2CD4_23_cytokines.xlsx

data := 23CD4TmemCD69
cytokines := 02CD4

som_dim := 7
nmetaclusts := 49


.PHONY: Analysis_block_4_cytokine_main_$(data)_$(cytokines)_goal

define Analysis_block_4_cytokine_main_$(data)_$(cytokines)_rule
Analysis_block_4_cytokine_main_$(data)_$(cytokines)_goal:
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_4_cytokine_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" file_panel_cytokines="$(file_panel_cytokines)" data="$(data)" cytokines="$(cytokines)" som_dim="$(som_dim)" nmetaclusts="$(nmetaclusts)"
endef

$(eval $(call Analysis_block_4_cytokine_main_$(data)_$(cytokines)_rule,))



####################################################


# ----------------------------------------------------------------------------------------------------
# Analysis of CK_2016-06-23_02_CD8_merging2 using panel2CD8.xlsx
# Use Analysis block 1
# ----------------------------------------------------------------------------------------------------


data_dir := CK_2016-06-23_02_CD8_merging2
file_panel := $(PANELS)/panel2CD8.xlsx
file_metadata := $(METADATA)/metadata_23_02.xlsx

data := 23CD8
panel := 02CD8
pca := pca1
nmetaclusts := 20 7

rand_seed_consensus := 1234


.PHONY: Analysis_block_1_main_$(data)_$(panel)_goal
Analysis_block_1_main_$(data)_$(panel)_goal: $(foreach i,$(nmetaclusts),Analysis_block_1_main_$(data)_$(panel)_goal$(i))

define Analysis_block_1_main_rule
Analysis_block_1_main_$(data)_$(panel)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_1_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_panel="$(file_panel)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" rand_seed_consensus="$(rand_seed_consensus)" nmetaclusts="$(1)"
endef
$(foreach i,$(nmetaclusts),$(eval $(call Analysis_block_1_main_rule,$(i))))


# --------------------------------------------------
# Analysis of CK_2016-06-23_02_CD8_merging2 cluster merging
# Use Analysis block 2
# --------------------------------------------------

nmetaclusts := 7
merging := merging2

.PHONY: Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal: $(foreach i,$(merging),Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(i))

define Analysis_block_2_cluster_merging_rule
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_2_cluster_merging.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" nmetaclusts="$(nmetaclusts)" merging="$(1)"
endef
$(foreach i,$(merging),$(eval $(call Analysis_block_2_cluster_merging_rule,$(i))))


####################################################
# CYTOKINE ANALYSIS


# --------------------------------------------------
# Analysis of CK_2016-06-23_02_CD8_merging2 cytokine extracting
# Use Analysis block 3
# --------------------------------------------------

file_extract_marker := $(PANELS)/panel2CD8_23_cytokines_CD69.xlsx
extract_marker := CD69

merging := merging2
extract_cluster := Tmem

.PHONY: Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_goal

define Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_rule
Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_goal:
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_3_cytokine_cell_extracting.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" merging="$(merging)" extract_cluster="$(extract_cluster)" file_extract_marker="$(file_extract_marker)" extract_marker="$(extract_marker)"
endef

$(eval $(call Analysis_block_3_cytokine_cell_extracting_$(data)_$(panel)_rule,))


# --------------------------------------------------
# Analysis of CK_2016-06-23_02_CD8_merging2_Tmem_merging2_CD69 cytokine analysis
# Use Analysis block 4
# --------------------------------------------------

data_dir := CK_2016-06-23_02_CD8_merging2_Tmem_merging2_CD69
file_metadata := $(METADATA)/metadata_23_02.xlsx
file_panel_cytokines := $(PANELS)/panel2CD8_23_cytokines.xlsx

data := 23CD8TmemCD69
cytokines := 02CD8

som_dim := 7
nmetaclusts := 49


.PHONY: Analysis_block_4_cytokine_main_$(data)_$(cytokines)_goal

define Analysis_block_4_cytokine_main_$(data)_$(cytokines)_rule
Analysis_block_4_cytokine_main_$(data)_$(cytokines)_goal:
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_4_cytokine_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" file_panel_cytokines="$(file_panel_cytokines)" data="$(data)" cytokines="$(cytokines)" som_dim="$(som_dim)" nmetaclusts="$(nmetaclusts)"
endef

$(eval $(call Analysis_block_4_cytokine_main_$(data)_$(cytokines)_rule,))



####################################################

###############################################################################################################
# Analysis of CK_2016-06-23_03 data
# Use Analysis block 1
###############################################################################################################

data_dir := CK_2016-06-23_03
file_panel := $(PANELS)/panel3.xlsx
file_metadata := $(METADATA)/metadata_23_03.xlsx

data := 23
panel := 03
pca := pca1
nmetaclusts := 20

rand_seed_consensus := 1234


.PHONY: Analysis_block_1_main_$(data)_$(panel)_goal
Analysis_block_1_main_$(data)_$(panel)_goal: $(foreach i,$(nmetaclusts),Analysis_block_1_main_$(data)_$(panel)_goal$(i))

define Analysis_block_1_main_rule
Analysis_block_1_main_$(data)_$(panel)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_1_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_panel="$(file_panel)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" rand_seed_consensus="$(rand_seed_consensus)" nmetaclusts="$(1)"
endef
$(foreach i,$(nmetaclusts),$(eval $(call Analysis_block_1_main_rule,$(i))))


# --------------------------------------------------
# Analysis of CK_2016-06-23_03 cluster merging
# Use Analysis block 2
# --------------------------------------------------

nmetaclusts := 20
merging := merging5


.PHONY: Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal: $(foreach i,$(merging),Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(i))


define Analysis_block_2_cluster_merging_rule
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_2_cluster_merging.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" nmetaclusts="$(nmetaclusts)" merging="$(1)"
endef
$(foreach i,$(merging),$(eval $(call Analysis_block_2_cluster_merging_rule,$(i))))



###############################################################################################################
# Analysis of CK_2016-06-23_03all data
# Use Analysis block 1
###############################################################################################################

data_dir := CK_2016-06-23_03all
file_panel := $(PANELS)/panel3_v2.xlsx
file_metadata := $(METADATA)/metadata_23_03all.xlsx

data := 23all
panel := 03v2
pca := pca1
nmetaclusts := 20

rand_seed_consensus := 1234


.PHONY: Analysis_block_1_main_$(data)_$(panel)_goal
Analysis_block_1_main_$(data)_$(panel)_goal: $(foreach i,$(nmetaclusts),Analysis_block_1_main_$(data)_$(panel)_goal$(i))

define Analysis_block_1_main_rule
Analysis_block_1_main_$(data)_$(panel)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_1_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_panel="$(file_panel)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" rand_seed_consensus="$(rand_seed_consensus)" nmetaclusts="$(1)"
endef
$(foreach i,$(nmetaclusts),$(eval $(call Analysis_block_1_main_rule,$(i))))


# --------------------------------------------------
# Analysis of CK_2016-06-23_03all cluster merging
# Use Analysis block 2
# --------------------------------------------------

nmetaclusts := 20
merging := merging5

.PHONY: Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal: $(foreach i,$(merging),Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(i))

define Analysis_block_2_cluster_merging_rule
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_2_cluster_merging.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" nmetaclusts="$(nmetaclusts)" merging="$(1)"
endef
$(foreach i,$(merging),$(eval $(call Analysis_block_2_cluster_merging_rule,$(i))))





































#
