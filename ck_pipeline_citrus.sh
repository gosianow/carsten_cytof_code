#!/bin/bash

###############################################################################################################
## Define paths to software and reference files

RWD_MAIN=/Users/gosia/Dropbox/UZH/carsten_cytof
RCODE=/Users/gosia/Dropbox/UZH/carsten_cytof_code
METADATA=$RWD_MAIN/CK_metadata
PANELS=$RWD_MAIN/CK_panels

###############################################################################################################
# Run citrus on CK_2016-06-23_01
###############################################################################################################

# ---------------------------------------
### First run that manually

## Copy the fcs files from CK_2016-06-23_01/010_cleanfcs/ info CK_2016-06-23_01/0_citrus/

## In R:
# library("citrus")
# citrus.checkFileParameterConsistencyUI()

## Generate  CK_panels/panel1.csv file from CK_panels/panel1.xlsx. panel1.csv is a data frame with two columns without column names. The first column corresponds to fcs_colname elements and the second to Antigen elements from panel1.xlsx for which transform equals to 1.

## In R:
# library("cytofCore")
# cytofCore.updatePanel(fcsFolder = "/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01/0_citrus", templateFile = "/Users/gosia/Dropbox/UZH/carsten_cytof/CK_panels/panel1.csv")

## Split base and tx samples from relabeled/ into relabeled_base/ and relabeled_tx/ and remove HD samples

# ---------------------------------------

data_dir="CK_2016-06-23_01"

file_panel="panel1.xlsx"
file_metadata="metadata_23_01.xlsx"

prefix_data="23_"
prefix_panel="01_"
prefix_pca="pca1_"

RWD=$RWD_MAIN/${data_dir}
ROUT=$RWD/Rout
mkdir -p $ROUT

R CMD BATCH --no-save --no-restore "--args rwd='$RWD' citrus_outdir='citrusOutput_1000'  path_clustering_observables='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' path_metadata='${METADATA}/${file_metadata}' fileSampleSize=1000" $RCODE/0_runCitrus.R $ROUT/0_runCitrus.Rout
tail $ROUT/0_runCitrus.Rout


R CMD BATCH --no-save --no-restore "--args rwd='$RWD' citrus_outdir='citrusOutput_2000'  path_clustering_observables='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' path_metadata='${METADATA}/${file_metadata}' fileSampleSize=2000" $RCODE/0_runCitrus.R $ROUT/0_runCitrus.Rout
tail $ROUT/0_runCitrus.Rout

###############################################################################################################
# Run citrus on CK_2016-06-23_03
###############################################################################################################

# ---------------------------------------
### First run that manually

## Copy the fcs files from CK_2016-06-23_03/010_cleanfcs/ info CK_2016-06-23_03/0_citrus/

## In R:
# library("citrus")
# citrus.checkFileParameterConsistencyUI()

## Generate  CK_panels/panel3.csv file from CK_panels/panel3.xlsx. panel3.csv is a data frame with two columns without column names. The first column corresponds to fcs_colname elements and the second to Antigen elements from panel3.xlsx for which transform equals to 1.

## In R:
# library("cytofCore")
# cytofCore.updatePanel(fcsFolder = "/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_03/0_citrus", templateFile = "/Users/gosia/Dropbox/UZH/carsten_cytof/CK_panels/panel3.csv")

## Split base and tx samples from relabeled/ into relabeled_base/ and relabeled_tx/ and remove HD samples

# ---------------------------------------

data_dir="CK_2016-06-23_03"

file_panel="panel3.xlsx"
file_metadata="metadata_23_03.xlsx"

prefix_data="23_"
prefix_panel="03_"
prefix_pca="pca1_"

RWD=$RWD_MAIN/${data_dir}
ROUT=$RWD/Rout
mkdir -p $ROUT

R CMD BATCH --no-save --no-restore "--args rwd='$RWD' citrus_outdir='citrusOutput_1000'  path_clustering_observables='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' path_metadata='${METADATA}/${file_metadata}' fileSampleSize=1000" $RCODE/0_runCitrus.R $ROUT/0_runCitrus.Rout
tail $ROUT/0_runCitrus.Rout


R CMD BATCH --no-save --no-restore "--args rwd='$RWD' citrus_outdir='citrusOutput_2000'  path_clustering_observables='030_heatmaps/${prefix_data}${prefix_panel}${prefix_pca}clustering_observables.xls' path_metadata='${METADATA}/${file_metadata}' fileSampleSize=2000" $RCODE/0_runCitrus.R $ROUT/0_runCitrus.Rout
tail $ROUT/0_runCitrus.Rout





#
