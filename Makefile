
MAKEARGS :=

R := R CMD BATCH --no-restore --no-save

RWD_MAIN := ../carsten_cytof/PD1_project
RCODE := .


## Define the default rule (makefiles are usually written so that the first target is for compiling the entire program)
.PHONY: all
all: ck_pipeline23 ck_pipeline29 ck_pipeline_merging_2datasets 


### Make sure no intermediate files are deleted
.SECONDARY:


###############################################################################################################
# Run the main ck pipelines
###############################################################################################################

### Analysis of data 23

.PHONY: ck_pipeline23

ck_pipeline23:
	echo "\n>> make ck_pipeline23"
	make -f ck_pipeline23.mk MAKEARGS="$(MAKEARGS)" R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)"


### Analysis of data 29

.PHONY: ck_pipeline29

ck_pipeline29:
	echo "\n>> make ck_pipeline29"
	make -f ck_pipeline29.mk MAKEARGS="$(MAKEARGS)" R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)"


### Merging data 23 and data 29


.PHONY: ck_pipeline_merging_2datasets

ck_pipeline_merging_2datasets:
	echo "\n>> make ck_pipeline_merging_2datasets"
	make -f ck_pipeline_merging_2datasets.mk MAKEARGS="$(MAKEARGS)" R="$(R)" RWD_MAIN="$(RWD_MAIN)" RCODE="$(RCODE)"



#
