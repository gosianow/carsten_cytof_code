
MAKEARGS :=

R := R CMD BATCH --no-restore --no-save

RWD_MAIN := ../carsten_cytof/PD1_project
RCODE := .

METADATA := $(RWD_MAIN)/CK_metadata
PANELS := $(RWD_MAIN)/CK_panels


## Define the default rule (makefiles are usually written so that the first target is for compiling the entire program)
.PHONY: all
all: Analysis_block_8_merging_2datasets_main_01_goal1 Analysis_block_8_merging_2datasets_main_01_goal2 \
	Analysis_block_8_merging_2datasets_main_01CD4_goal Analysis_block_8_merging_2datasets_main_01CD8_goal

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


.PHONY: Analysis_block_8_merging_2datasets_main_$(panel)_goal1

Analysis_block_8_merging_2datasets_main_$(panel)_goal1:
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_8_merging_2datasets_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir_merged="$(data_dir_merged)" data_dir1="$(data_dir1)" data_dir2="$(data_dir2)" file_metadata1="$(file_metadata1)" file_metadata2="$(file_metadata2)" data1="$(data1)" data2="$(data2)" panel="$(panel)" pca1="$(pca1)" pca2="$(pca2)" merging1="$(merging1)" merging2="$(merging2)"


merging1 := merging7
merging2 := merging4

.PHONY: Analysis_block_8_merging_2datasets_main_$(panel)_goal2

Analysis_block_8_merging_2datasets_main_$(panel)_goal2:
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_8_merging_2datasets_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir_merged="$(data_dir_merged)" data_dir1="$(data_dir1)" data_dir2="$(data_dir2)" file_metadata1="$(file_metadata1)" file_metadata2="$(file_metadata2)" data1="$(data1)" data2="$(data2)" panel="$(panel)" pca1="$(pca1)" pca2="$(pca2)" merging1="$(merging1)" merging2="$(merging2)"


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

merging1 := merging6
merging2 := merging6


.PHONY: Analysis_block_8_merging_2datasets_main_$(panel)_goal

Analysis_block_8_merging_2datasets_main_$(panel)_goal:
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_8_merging_2datasets_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir_merged="$(data_dir_merged)" data_dir1="$(data_dir1)" data_dir2="$(data_dir2)" file_metadata1="$(file_metadata1)" file_metadata2="$(file_metadata2)" data1="$(data1)" data2="$(data2)" panel="$(panel)" pca1="$(pca1)" pca2="$(pca2)" merging1="$(merging1)" merging2="$(merging2)"



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

merging1 := merging7
merging2 := merging7


.PHONY: Analysis_block_8_merging_2datasets_main_$(panel)_goal


Analysis_block_8_merging_2datasets_main_$(panel)_goal:
	echo "\n>> make"
	make $(MAKEARGS) -f Analysis_block_8_merging_2datasets_main.mk R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)" data_dir_merged="$(data_dir_merged)" data_dir1="$(data_dir1)" data_dir2="$(data_dir2)" file_metadata1="$(file_metadata1)" file_metadata2="$(file_metadata2)" data1="$(data1)" data2="$(data2)" panel="$(panel)" pca1="$(pca1)" pca2="$(pca2)" merging1="$(merging1)" merging2="$(merging2)"



###############################################################################################################
# Merging data 23 and 29 for panel 1 CD4 --- no HD10
###############################################################################################################



###############################################################################################################
# Merging data 23 and 29 for panel 1 CD8 --- no HD10
###############################################################################################################




###############################################################################################################
# Merging data 23 and 29 for panel 3
###############################################################################################################




###############################################################################################################
# Merging data 23all and 29all3 for panel3_v2
###############################################################################################################









###############################################################################################################
# Merging data 23 and 29 for panel 2 CD4 cytokines
###############################################################################################################





###############################################################################################################
# Merging data 23 and 29 for panel 2 CD8 cytokines
###############################################################################################################























#
