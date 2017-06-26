
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
	Analysis_block_8_merging_2datasets_cytokine_23CD4TmemCD69_29CD4TmemCD69_02CD4_goal \
	Analysis_block_8_merging_2datasets_cytokine_23CD8TmemCD69_29CD8TmemCD69_02CD8_goal \
	Analysis_block_8_merging_2datasets_main_03_merging4_merging2_goal \
	Analysis_block_8_merging_2datasets_main_03v2_merging4_merging1_goal



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

merging1 := merging4
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

merging1 := merging4
merging2 := merging1

.PHONY: Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_goal

define Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_rule
Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_goal:
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_8_merging_2datasets_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir_merged="$(data_dir_merged)" data_dir1="$(data_dir1)" data_dir2="$(data_dir2)" file_metadata1="$(file_metadata1)" file_metadata2="$(file_metadata2)" data1="$(data1)" data2="$(data2)" panel="$(panel)" pca1="$(pca1)" pca2="$(pca2)" merging1="$(merging1)" merging2="$(merging2)"
endef

$(eval $(call Analysis_block_8_merging_2datasets_main_$(panel)_$(merging1)_$(merging2)_rule,))



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


.PHONY: Analysis_block_8_merging_2datasets_cytokine_$(data1)_$(data2)_$(cytokines)_goal

define Analysis_block_8_merging_2datasets_cytokine_$(data1)_$(data2)_$(cytokines)_rule
Analysis_block_8_merging_2datasets_cytokine_$(data1)_$(data2)_$(cytokines)_goal:
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_8_merging_2datasets_cytokine.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir_merged="$(data_dir_merged)" data_dir1="$(data_dir1)" data_dir2="$(data_dir2)" file_metadata1="$(file_metadata1)" file_metadata2="$(file_metadata2)" data1="$(data1)" data2="$(data2)" cytokines="$(cytokines)" som_dim="$(som_dim)" nmetaclusts="$(nmetaclusts)"
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


.PHONY: Analysis_block_8_merging_2datasets_cytokine_$(data1)_$(data2)_$(cytokines)_goal

define Analysis_block_8_merging_2datasets_cytokine_$(data1)_$(data2)_$(cytokines)_rule
Analysis_block_8_merging_2datasets_cytokine_$(data1)_$(data2)_$(cytokines)_goal:
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_8_merging_2datasets_cytokine.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir_merged="$(data_dir_merged)" data_dir1="$(data_dir1)" data_dir2="$(data_dir2)" file_metadata1="$(file_metadata1)" file_metadata2="$(file_metadata2)" data1="$(data1)" data2="$(data2)" cytokines="$(cytokines)" som_dim="$(som_dim)" nmetaclusts="$(nmetaclusts)"
endef

$(eval $(call Analysis_block_8_merging_2datasets_cytokine_$(data1)_$(data2)_$(cytokines)_rule,))





















#
