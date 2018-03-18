
MAKEARGS :=

R := R CMD BATCH --no-restore --no-save

RWD_MAIN := ../carsten_cytof/PD1_project
RCODE := .

METADATA := $(RWD_MAIN)/CK_metadata
PANELS := $(RWD_MAIN)/CK_panels


## Define the default rule (makefiles are usually written so that the first target is for compiling the entire program)
.PHONY: all
all: Analysis_block_8_merging_2datasets_main_01_merging6_merging4_goal Analysis_block_8_merging_2datasets_main_01_merging7_merging4_goal \
	Analysis_block_8_merging_2datasets_main_01CD4_merging5_merging5_goal \
	Analysis_block_8_merging_2datasets_main_01CD8_merging6_merging6_goal \
	Analysis_block_8_merging_2datasets_main_01CD4_merging5_merging5_goal_noHD10 \
	Analysis_block_8_merging_2datasets_main_01CD8_merging6_merging6_goal_noHD10 \
	Analysis_block_8_merging_2datasets_main_03_merging5_merging2_goal \
	Analysis_block_8_merging_2datasets_main_03v2_merging5_merging1_goal \
	Analysis_block_8_merging_2datasets_cytokine_23CD4TmemCD69_29CD4TmemCD69_02CD4_goal \
	Analysis_block_8_merging_2datasets_cytokine_23CD8TmemCD69_29CD8TmemCD69_02CD8_goal \
	Analysis_block_8_merging_2datasets_cytokine_02CD4v2_goal \
	Analysis_block_8_merging_2datasets_cytokine_02CD8v2_goal \
	Analysis_block_8_merging_2datasets_cytokine_02Tregsv2_goal




### Make sure no intermediate files are deleted
.SECONDARY:


###############################################################################################################
# Merging data 23 and 29 for panel 1
###############################################################################################################

data_dir_merged := CK_2016-06-merged_23_29/01

data_dir1 := CK_2016-06-23_01
data_dir2 := CK_2016-06-29_01
file_metadata1 := $(METADATA)/metadata_23_01.xlsx
file_metadata2 := $(METADATA)/metadata_29_01.xlsx

data1 := 23
data2 := 29
pca1 := pca1
pca2 := pca1

panel := 01

merging1 := merging6
merging2 := merging4


.PHONY: Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_goal

define Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_rule
Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_goal:
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_8_merging_2datasets_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir_merged="$(data_dir_merged)" data_dir1="$(data_dir1)" data_dir2="$(data_dir2)" file_metadata1="$(file_metadata1)" file_metadata2="$(file_metadata2)" data1="$(data1)" data2="$(data2)" panel="$(panel)" pca1="$(pca1)" pca2="$(pca2)" merging1="$(merging1)" merging2="$(merging2)"
endef

$(eval $(call Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_rule,))


merging1 := merging7
merging2 := merging4

.PHONY: Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_goal

define Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_rule
Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_goal:
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_8_merging_2datasets_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir_merged="$(data_dir_merged)" data_dir1="$(data_dir1)" data_dir2="$(data_dir2)" file_metadata1="$(file_metadata1)" file_metadata2="$(file_metadata2)" data1="$(data1)" data2="$(data2)" panel="$(panel)" pca1="$(pca1)" pca2="$(pca2)" merging1="$(merging1)" merging2="$(merging2)"
endef

$(eval $(call Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_rule,))


###############################################################################################################
# Merging data 23 and 29 for panel 1 CD4
###############################################################################################################

data_dir_merged := CK_2016-06-merged_23_29/01_CD4

data_dir1 := CK_2016-06-23_01_CD4_mergingNEW2
data_dir2 := CK_2016-06-29_01_CD4_merging2
file_metadata1 := $(METADATA)/metadata_23_01.xlsx
file_metadata2 := $(METADATA)/metadata_29_01.xlsx

data1 := 23CD4
data2 := 29CD4
pca1 := pca1
pca2 := pca1

panel := 01CD4

merging1 := merging5
merging2 := merging5


.PHONY: Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_goal

define Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_rule
Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_goal:
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_8_merging_2datasets_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir_merged="$(data_dir_merged)" data_dir1="$(data_dir1)" data_dir2="$(data_dir2)" file_metadata1="$(file_metadata1)" file_metadata2="$(file_metadata2)" data1="$(data1)" data2="$(data2)" panel="$(panel)" pca1="$(pca1)" pca2="$(pca2)" merging1="$(merging1)" merging2="$(merging2)"
endef

$(eval $(call Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_rule,))



###############################################################################################################
# Merging data 23 and 29 for panel 1 CD8
###############################################################################################################

data_dir_merged := CK_2016-06-merged_23_29/01_CD8

data_dir1 := CK_2016-06-23_01_CD8_mergingNEW2
data_dir2 := CK_2016-06-29_01_CD8_merging2
file_metadata1 := $(METADATA)/metadata_23_01.xlsx
file_metadata2 := $(METADATA)/metadata_29_01.xlsx

data1 := 23CD8
data2 := 29CD8
pca1 := pca1
pca2 := pca1

panel := 01CD8

merging1 := merging6
merging2 := merging6


.PHONY: Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_goal

define Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_rule
Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_goal:
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_8_merging_2datasets_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir_merged="$(data_dir_merged)" data_dir1="$(data_dir1)" data_dir2="$(data_dir2)" file_metadata1="$(file_metadata1)" file_metadata2="$(file_metadata2)" data1="$(data1)" data2="$(data2)" panel="$(panel)" pca1="$(pca1)" pca2="$(pca2)" merging1="$(merging1)" merging2="$(merging2)"
endef

$(eval $(call Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_rule,))



###############################################################################################################
# Merging data 23 and 29 for panel 1 CD4 --- no HD10
###############################################################################################################

data_dir_merged := CK_2016-06-merged_23_29/01_CD4_noHD10

data_dir1 := CK_2016-06-23_01_CD4_mergingNEW2
data_dir2 := CK_2016-06-29_01_CD4_merging2
file_metadata1 := $(METADATA)/metadata_23_01.xlsx
file_metadata2 := $(METADATA)/metadata_29_01_noHD10.xlsx

data1 := 23CD4
data2 := 29CD4
pca1 := pca1
pca2 := pca1

panel := 01CD4

merging1 := merging5
merging2 := merging5


.PHONY: Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_goal_noHD10

define Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_rule
Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_goal_noHD10:
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_8_merging_2datasets_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir_merged="$(data_dir_merged)" data_dir1="$(data_dir1)" data_dir2="$(data_dir2)" file_metadata1="$(file_metadata1)" file_metadata2="$(file_metadata2)" data1="$(data1)" data2="$(data2)" panel="$(panel)" pca1="$(pca1)" pca2="$(pca2)" merging1="$(merging1)" merging2="$(merging2)"
endef

$(eval $(call Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_rule,))


###############################################################################################################
# Merging data 23 and 29 for panel 1 CD8 --- no HD10
###############################################################################################################

data_dir_merged := CK_2016-06-merged_23_29/01_CD8_noHD10

data_dir1 := CK_2016-06-23_01_CD8_mergingNEW2
data_dir2 := CK_2016-06-29_01_CD8_merging2
file_metadata1 := $(METADATA)/metadata_23_01.xlsx
file_metadata2 := $(METADATA)/metadata_29_01_noHD10.xlsx

data1 := 23CD8
data2 := 29CD8
pca1 := pca1
pca2 := pca1

panel := 01CD8

merging1 := merging6
merging2 := merging6


.PHONY: Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_goal_noHD10

define Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_rule
Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_goal_noHD10:
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_8_merging_2datasets_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir_merged="$(data_dir_merged)" data_dir1="$(data_dir1)" data_dir2="$(data_dir2)" file_metadata1="$(file_metadata1)" file_metadata2="$(file_metadata2)" data1="$(data1)" data2="$(data2)" panel="$(panel)" pca1="$(pca1)" pca2="$(pca2)" merging1="$(merging1)" merging2="$(merging2)"
endef

$(eval $(call Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_rule,))



###############################################################################################################
# Merging data 23 and 29 for panel 3
###############################################################################################################


data_dir_merged := CK_2016-06-merged_23_29/03

data_dir1 := CK_2016-06-23_03
data_dir2 := CK_2016-06-29_03
file_metadata1 := $(METADATA)/metadata_23_03.xlsx
file_metadata2 := $(METADATA)/metadata_29_03.xlsx

data1 := 23
data2 := 29
pca1 := pca1
pca2 := pca1

panel := 03

merging1 := merging5
merging2 := merging2


.PHONY: Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_goal

define Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_rule
Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_goal:
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_8_merging_2datasets_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir_merged="$(data_dir_merged)" data_dir1="$(data_dir1)" data_dir2="$(data_dir2)" file_metadata1="$(file_metadata1)" file_metadata2="$(file_metadata2)" data1="$(data1)" data2="$(data2)" panel="$(panel)" pca1="$(pca1)" pca2="$(pca2)" merging1="$(merging1)" merging2="$(merging2)"
endef

$(eval $(call Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_rule,))



###############################################################################################################
# Merging data 23all and 29all3 for panel3_v2
###############################################################################################################

data_dir_merged := CK_2016-06-merged_23_29/03all

data_dir1 := CK_2016-06-23_03all
data_dir2 := CK_2016-06-29_03all3
file_metadata1 := $(METADATA)/metadata_23_03all.xlsx
file_metadata2 := $(METADATA)/metadata_29_03all3.xlsx

data1 := 23all
data2 := 29all3
pca1 := pca1
pca2 := pca1

panel := 03v2

merging1 := merging5
merging2 := merging1

.PHONY: Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_goal

define Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_rule
Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_goal:
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_8_merging_2datasets_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir_merged="$(data_dir_merged)" data_dir1="$(data_dir1)" data_dir2="$(data_dir2)" file_metadata1="$(file_metadata1)" file_metadata2="$(file_metadata2)" data1="$(data1)" data2="$(data2)" panel="$(panel)" pca1="$(pca1)" pca2="$(pca2)" merging1="$(merging1)" merging2="$(merging2)"
endef

$(eval $(call Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_rule,))



###############################################################################################################
# Merging data 23 and 29 for panel 2v2 CD4 cytokines
###############################################################################################################

## arguments for the specific data set
data_dir_merged := CK_2016-06-merged_23_29/02v2_CD4

extract_marker := all
extract_cluster := all naive CM

data_dir1 := CK_2016-06-23_02_CD4_merging3
data_dir2 := CK_2016-06-29_02_CD4_merging3
file_metadata1 := $(METADATA)/metadata_23_02.xlsx
file_metadata2 := $(METADATA)/metadata_29_02.xlsx

data1 := 23CD4
data2 := 29CD4

cytokines := 02CD4v2

som_dim := 7
nmetaclusts := 49

data_dir_back1 := CK_2016-06-23_02_CD4_merging3
data_dir_back2 := CK_2016-06-29_02_CD4_merging3
data_back1 := 23CD4
data_back2 := 29CD4
panel_back1 := 02CD4v2
panel_back2 := 02CD4v2
pca_back1 := pca0
pca_back2 := pca0
merging_back1 := merging1
merging_back2 := merging1


.PHONY: Analysis_block_8_merging_2datasets_cytokine_$(cytokines)_goal
Analysis_block_8_merging_2datasets_cytokine_$(cytokines)_goal: $(foreach i,$(extract_cluster),Analysis_block_8_merging_2datasets_cytokine_$(cytokines)_goal$(i))


define Analysis_block_8_merging_2datasets_cytokine_$(cytokines)_rule
Analysis_block_8_merging_2datasets_cytokine_$(cytokines)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_8_merging_2datasets_cytokine.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir_merged="$(data_dir_merged)" data_dir1="$(data_dir1)_$(1)_$(merging_back1)_$(extract_marker)" data_dir2="$(data_dir2)_$(1)_$(merging_back2)_$(extract_marker)" file_metadata1="$(file_metadata1)" file_metadata2="$(file_metadata2)" data1="$(data1)$(1)$(extract_marker)" data2="$(data2)$(1)$(extract_marker)" cytokines="$(cytokines)" som_dim="$(som_dim)" nmetaclusts="$(nmetaclusts)" \
	data_dir_back1="$(data_dir_back1)" data_dir_back2="$(data_dir_back2)" data_back1="$(data_back1)" data_back2="$(data_back2)" panel_back1="$(panel_back1)" panel_back2="$(panel_back2)" pca_back1="$(pca_back1)" pca_back2="$(pca_back2)" merging_back1="$(merging_back1)" merging_back2="$(merging_back2)"
endef
$(foreach i,$(extract_cluster),$(eval $(call Analysis_block_8_merging_2datasets_cytokine_$(cytokines)_rule,$(i))))





###############################################################################################################
# Merging data 23 and 29 for panel 2v2 CD8 cytokines
###############################################################################################################


## arguments for the specific data set
data_dir_merged := CK_2016-06-merged_23_29/02v2_CD8

extract_marker := all
extract_cluster := all naive CM

data_dir1 := CK_2016-06-23_02_CD8_merging3
data_dir2 := CK_2016-06-29_02_CD8_merging3
file_metadata1 := $(METADATA)/metadata_23_02.xlsx
file_metadata2 := $(METADATA)/metadata_29_02.xlsx

data1 := 23CD8
data2 := 29CD8

cytokines := 02CD8v2

som_dim := 7
nmetaclusts := 49

data_dir_back1 := CK_2016-06-23_02_CD8_merging3
data_dir_back2 := CK_2016-06-29_02_CD8_merging3
data_back1 := 23CD8
data_back2 := 29CD8
panel_back1 := 02CD8v2
panel_back2 := 02CD8v2
pca_back1 := pca0
pca_back2 := pca0
merging_back1 := merging1
merging_back2 := merging1


.PHONY: Analysis_block_8_merging_2datasets_cytokine_$(cytokines)_goal
Analysis_block_8_merging_2datasets_cytokine_$(cytokines)_goal: $(foreach i,$(extract_cluster),Analysis_block_8_merging_2datasets_cytokine_$(cytokines)_goal$(i))


define Analysis_block_8_merging_2datasets_cytokine_$(cytokines)_rule
Analysis_block_8_merging_2datasets_cytokine_$(cytokines)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_8_merging_2datasets_cytokine.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir_merged="$(data_dir_merged)" data_dir1="$(data_dir1)_$(1)_$(merging_back1)_$(extract_marker)" data_dir2="$(data_dir2)_$(1)_$(merging_back2)_$(extract_marker)" file_metadata1="$(file_metadata1)" file_metadata2="$(file_metadata2)" data1="$(data1)$(1)$(extract_marker)" data2="$(data2)$(1)$(extract_marker)" cytokines="$(cytokines)" som_dim="$(som_dim)" nmetaclusts="$(nmetaclusts)" \
	data_dir_back1="$(data_dir_back1)" data_dir_back2="$(data_dir_back2)" data_back1="$(data_back1)" data_back2="$(data_back2)" panel_back1="$(panel_back1)" panel_back2="$(panel_back2)" pca_back1="$(pca_back1)" pca_back2="$(pca_back2)" merging_back1="$(merging_back1)" merging_back2="$(merging_back2)"
endef
$(foreach i,$(extract_cluster),$(eval $(call Analysis_block_8_merging_2datasets_cytokine_$(cytokines)_rule,$(i))))



###############################################################################################################
# Merging data 23 and 29 for panel 2 CD4 cytokines
###############################################################################################################

## arguments for the specific data set
data_dir_merged := CK_2016-06-merged_23_29/02_CD4

data_dir1 := CK_2016-06-23_02_CD4_merging2_Tmem_merging2_CD69
data_dir2 := CK_2016-06-29_02_CD4_merging_Tmem_merging3_CD69
file_metadata1 := $(METADATA)/metadata_23_02.xlsx
file_metadata2 := $(METADATA)/metadata_29_02.xlsx

data1 := 23CD4TmemCD69
data2 := 29CD4TmemCD69

cytokines := 02CD4

som_dim := 7
nmetaclusts := 49

data_dir_back1 := CK_2016-06-23_02_CD4_merging2
data_dir_back2 := CK_2016-06-29_02_CD4_merging
data_back1 := 23CD4
data_back2 := 29CD4
panel_back1 := 02CD4
panel_back2 := 02CD4
pca_back1 := pca1
pca_back2 := pca1
merging_back1 := merging2
merging_back2 := merging3


.PHONY: Analysis_block_8_merging_2datasets_cytokine_$(data1)_$(data2)_$(cytokines)_goal

define Analysis_block_8_merging_2datasets_cytokine_$(data1)_$(data2)_$(cytokines)_rule
Analysis_block_8_merging_2datasets_cytokine_$(data1)_$(data2)_$(cytokines)_goal:
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_8_merging_2datasets_cytokine.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir_merged="$(data_dir_merged)" data_dir1="$(data_dir1)" data_dir2="$(data_dir2)" file_metadata1="$(file_metadata1)" file_metadata2="$(file_metadata2)" data1="$(data1)" data2="$(data2)" cytokines="$(cytokines)" som_dim="$(som_dim)" nmetaclusts="$(nmetaclusts)" \
	data_dir_back1="$(data_dir_back1)" data_dir_back2="$(data_dir_back2)" data_back1="$(data_back1)" data_back2="$(data_back2)" panel_back1="$(panel_back1)" panel_back2="$(panel_back2)" pca_back1="$(pca_back1)" pca_back2="$(pca_back2)" merging_back1="$(merging_back1)" merging_back2="$(merging_back2)"
endef

$(eval $(call Analysis_block_8_merging_2datasets_cytokine_$(data1)_$(data2)_$(cytokines)_rule,))




###############################################################################################################
# Merging data 23 and 29 for panel 2 CD8 cytokines
###############################################################################################################


## arguments for the specific data set
data_dir_merged := CK_2016-06-merged_23_29/02_CD8

data_dir1 := CK_2016-06-23_02_CD8_merging2_Tmem_merging2_CD69
data_dir2 := CK_2016-06-29_02_CD8_merging_Tmem_merging3_CD69
file_metadata1 := $(METADATA)/metadata_23_02.xlsx
file_metadata2 := $(METADATA)/metadata_29_02.xlsx

data1 := 23CD8TmemCD69
data2 := 29CD8TmemCD69

cytokines := 02CD8

som_dim := 7
nmetaclusts := 49

data_dir_back1 := CK_2016-06-23_02_CD8_merging2
data_dir_back2 := CK_2016-06-29_02_CD8_merging
data_back1 := 23CD8
data_back2 := 29CD8
panel_back1 := 02CD8
panel_back2 := 02CD8
pca_back1 := pca1
pca_back2 := pca1
merging_back1 := merging2
merging_back2 := merging3

.PHONY: Analysis_block_8_merging_2datasets_cytokine_$(data1)_$(data2)_$(cytokines)_goal

define Analysis_block_8_merging_2datasets_cytokine_$(data1)_$(data2)_$(cytokines)_rule
Analysis_block_8_merging_2datasets_cytokine_$(data1)_$(data2)_$(cytokines)_goal:
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_8_merging_2datasets_cytokine.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir_merged="$(data_dir_merged)" data_dir1="$(data_dir1)" data_dir2="$(data_dir2)" file_metadata1="$(file_metadata1)" file_metadata2="$(file_metadata2)" data1="$(data1)" data2="$(data2)" cytokines="$(cytokines)" som_dim="$(som_dim)" nmetaclusts="$(nmetaclusts)" \
	data_dir_back1="$(data_dir_back1)" data_dir_back2="$(data_dir_back2)" data_back1="$(data_back1)" data_back2="$(data_back2)" panel_back1="$(panel_back1)" panel_back2="$(panel_back2)" pca_back1="$(pca_back1)" pca_back2="$(pca_back2)" merging_back1="$(merging_back1)" merging_back2="$(merging_back2)"
endef

$(eval $(call Analysis_block_8_merging_2datasets_cytokine_$(data1)_$(data2)_$(cytokines)_rule,))





###############################################################################################################
# Merging data 23 and 29 for panel 2v2 Tregs cytokines
###############################################################################################################

## arguments for the specific data set
data_dir_merged := CK_2016-06-merged_23_29/02v2_Tregs

extract_marker := all
extract_cluster := Tregs

data_dir1 := CK_2016-06-23_02
data_dir2 := CK_2016-06-29_02
file_metadata1 := $(METADATA)/metadata_23_02.xlsx
file_metadata2 := $(METADATA)/metadata_29_02.xlsx

data1 := 23
data2 := 29

cytokines := 02Tregsv2

som_dim := 3
nmetaclusts := 9

data_dir_back1 := CK_2016-06-23_02
data_dir_back2 := CK_2016-06-29_02
data_back1 := 23
data_back2 := 29
panel_back1 := 02v2
panel_back2 := 02v2
pca_back1 := pca0
pca_back2 := pca0
merging_back1 := merging3
merging_back2 := merging3


.PHONY: Analysis_block_8_merging_2datasets_cytokine_$(cytokines)_goal
Analysis_block_8_merging_2datasets_cytokine_$(cytokines)_goal: $(foreach i,$(extract_cluster),Analysis_block_8_merging_2datasets_cytokine_$(cytokines)_goal$(i))


define Analysis_block_8_merging_2datasets_cytokine_$(cytokines)_rule
Analysis_block_8_merging_2datasets_cytokine_$(cytokines)_goal$(1):
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_8_merging_2datasets_cytokine.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir_merged="$(data_dir_merged)" data_dir1="$(data_dir1)_$(1)_$(merging_back1)_$(extract_marker)" data_dir2="$(data_dir2)_$(1)_$(merging_back2)_$(extract_marker)" file_metadata1="$(file_metadata1)" file_metadata2="$(file_metadata2)" data1="$(data1)$(1)$(extract_marker)" data2="$(data2)$(1)$(extract_marker)" cytokines="$(cytokines)" som_dim="$(som_dim)" nmetaclusts="$(nmetaclusts)" \
	data_dir_back1="$(data_dir_back1)" data_dir_back2="$(data_dir_back2)" data_back1="$(data_back1)" data_back2="$(data_back2)" panel_back1="$(panel_back1)" panel_back2="$(panel_back2)" pca_back1="$(pca_back1)" pca_back2="$(pca_back2)" merging_back1="$(merging_back1)" merging_back2="$(merging_back2)"
endef
$(foreach i,$(extract_cluster),$(eval $(call Analysis_block_8_merging_2datasets_cytokine_$(cytokines)_rule,$(i))))














#
