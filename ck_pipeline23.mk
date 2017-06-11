
make_args := -f

R := R CMD BATCH --no-restore --no-save

RWD_MAIN := ../carsten_cytof/PD1_project
RCODE := .
METADATA := $(RWD_MAIN)/CK_metadata
PANELS := $(RWD_MAIN)/CK_panels


## Define the default rule (makefiles are usually written so that the first target is for compiling the entire program)
.PHONY: all
all: Analysis_block_1_main_23_01_goal Analysis_block_2_cluster_merging_23_01_cl20_goal Analysis_block_3_cluster_extracting_23_01_goal \
	Analysis_block_1_main_23CD4_01CD4_goal Analysis_block_1_main_23CD8_01CD8_goal \
	Analysis_block_1_main_23_02v2_goal Analysis_block_2_cluster_merging_23_02v2_cl20_goal \
	Analysis_block_1_main_23_03_goal Analysis_block_2_cluster_merging_23_03_cl20_goal \
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
	make $(make_args) Analysis_block_1_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_panel="$(file_panel)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" rand_seed_consensus="$(rand_seed_consensus)" nmetaclusts="$(1)"
endef
$(foreach i,$(nmetaclusts),$(eval $(call Analysis_block_1_main_rule,$(i))))


# --------------------------------------------------
# Analysis of CK_2016-06-23_01 cluster_merging
# Use Analysis block 2
# --------------------------------------------------


nmetaclusts := 20
merging := mergingNEW2 merging6


.PHONY: Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal: $(foreach i,$(merging),Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(i))


define Analysis_block_2_cluster_merging_rule
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(1):
	echo "\n>> make"
	make $(make_args) Analysis_block_2_cluster_merging.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" nmetaclusts="$(nmetaclusts)" merging="$(1)"
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
	make $(make_args) Analysis_block_3_cluster_extracting.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" merging="$(merging)" extract_cluster="$(1)"
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
nmetaclusts := 20

rand_seed_consensus := 1234


.PHONY: Analysis_block_1_main_$(data)_$(panel)_goal
Analysis_block_1_main_$(data)_$(panel)_goal: $(foreach i,$(nmetaclusts),Analysis_block_1_main_$(data)_$(panel)_goal$(i))

define Analysis_block_1_main_rule
Analysis_block_1_main_$(data)_$(panel)_goal$(1):
	echo "\n>> make"
	make $(make_args) Analysis_block_1_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_panel="$(file_panel)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" rand_seed_consensus="$(rand_seed_consensus)" nmetaclusts="$(1)"
endef
$(foreach i,$(nmetaclusts),$(eval $(call Analysis_block_1_main_rule,$(i))))


# --------------------------------------------------
# Analysis of CK_2016-06-23_01 cluster_merging
# Use Analysis block 2
# --------------------------------------------------


# nmetaclusts := 20
# merging := mergingNEW2 merging6
#
#
# .PHONY: Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal
# Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal: $(foreach i,$(merging),Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(i))
#
#
# define Analysis_block_2_cluster_merging_rule
# Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(1):
# 	echo "\n>> make"
# 	make $(make_args) Analysis_block_2_cluster_merging.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" nmetaclusts="$(nmetaclusts)" merging="$(1)"
# endef
# $(foreach i,$(merging),$(eval $(call Analysis_block_2_cluster_merging_rule,$(i))))


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
nmetaclusts := 20

rand_seed_consensus := 1234


.PHONY: Analysis_block_1_main_$(data)_$(panel)_goal
Analysis_block_1_main_$(data)_$(panel)_goal: $(foreach i,$(nmetaclusts),Analysis_block_1_main_$(data)_$(panel)_goal$(i))

define Analysis_block_1_main_rule
Analysis_block_1_main_$(data)_$(panel)_goal$(1):
	echo "\n>> make"
	make $(make_args) Analysis_block_1_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_panel="$(file_panel)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" rand_seed_consensus="$(rand_seed_consensus)" nmetaclusts="$(1)"
endef
$(foreach i,$(nmetaclusts),$(eval $(call Analysis_block_1_main_rule,$(i))))


# --------------------------------------------------
# Analysis of CK_2016-06-23_01 cluster_merging
# Use Analysis block 2
# --------------------------------------------------


# nmetaclusts := 20
# merging := mergingNEW2 merging6
#
#
# .PHONY: Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal
# Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal: $(foreach i,$(merging),Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(i))
#
#
# define Analysis_block_2_cluster_merging_rule
# Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(1):
# 	echo "\n>> make"
# 	make $(make_args) Analysis_block_2_cluster_merging.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" nmetaclusts="$(nmetaclusts)" merging="$(1)"
# endef
# $(foreach i,$(merging),$(eval $(call Analysis_block_2_cluster_merging_rule,$(i))))


###############################################################################################################
# Analysis of CK_2016-06-23_02 data
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
	make $(make_args) Analysis_block_1_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_panel="$(file_panel)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" rand_seed_consensus="$(rand_seed_consensus)" nmetaclusts="$(1)"
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
	make $(make_args) Analysis_block_2_cluster_merging.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" nmetaclusts="$(nmetaclusts)" merging="$(1)"
endef
$(foreach i,$(merging),$(eval $(call Analysis_block_2_cluster_merging_rule,$(i))))


# --------------------------------------------------
# CK_2016-06-23_02 - CD4 and CD8 cluster extracting from merging
# Use Analysis block 3
# --------------------------------------------------

# merging := merging
# extract_cluster := CD4 CD8
#
# .PHONY: Analysis_block_3_cluster_extracting_$(data)_$(panel)_goal
# Analysis_block_3_cluster_extracting_$(data)_$(panel)_goal: $(foreach i,$(extract_cluster),Analysis_block_3_cluster_extracting_$(data)_$(panel)_goal$(i))
#
# define Analysis_block_3_cluster_extracting_rule
# Analysis_block_3_cluster_extracting_$(data)_$(panel)_goal$(1):
# 	echo "\n>> make"
# 	make $(make_args) Analysis_block_3_cluster_extracting.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" merging="$(merging)" extract_cluster="$(1)"
# endef
# $(foreach i,$(extract_cluster),$(eval $(call Analysis_block_3_cluster_extracting_rule,$(i))))





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
	make $(make_args) Analysis_block_1_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_panel="$(file_panel)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" rand_seed_consensus="$(rand_seed_consensus)" nmetaclusts="$(1)"
endef
$(foreach i,$(nmetaclusts),$(eval $(call Analysis_block_1_main_rule,$(i))))


# --------------------------------------------------
# Analysis of CK_2016-06-23_03 cluster merging
# Use Analysis block 2
# --------------------------------------------------

nmetaclusts := 20
merging := merging3 merging4


.PHONY: Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal: $(foreach i,$(merging),Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(i))


define Analysis_block_2_cluster_merging_rule
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(1):
	echo "\n>> make"
	make $(make_args) Analysis_block_2_cluster_merging.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" nmetaclusts="$(nmetaclusts)" merging="$(1)"
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
	make $(make_args) Analysis_block_1_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_panel="$(file_panel)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" rand_seed_consensus="$(rand_seed_consensus)" nmetaclusts="$(1)"
endef
$(foreach i,$(nmetaclusts),$(eval $(call Analysis_block_1_main_rule,$(i))))


# --------------------------------------------------
# Analysis of CK_2016-06-23_03all cluster merging
# Use Analysis block 2
# --------------------------------------------------

nmetaclusts := 20
merging := merging4

.PHONY: Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal: $(foreach i,$(merging),Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(i))

define Analysis_block_2_cluster_merging_rule
Analysis_block_2_cluster_merging_$(data)_$(panel)_cl$(nmetaclusts)_goal$(1):
	echo "\n>> make"
	make $(make_args) Analysis_block_2_cluster_merging.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir="$(data_dir)" file_metadata="$(file_metadata)" data="$(data)" panel="$(panel)" pca="$(pca)" nmetaclusts="$(nmetaclusts)" merging="$(1)"
endef
$(foreach i,$(merging),$(eval $(call Analysis_block_2_cluster_merging_rule,$(i))))





































#
