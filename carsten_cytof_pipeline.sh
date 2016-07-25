#!/bin/bash
## Define paths to software and reference files

RWD_MAIN=/Users/gosia/Dropbox/UZH/carsten_cytof
RCODE=$RWD_MAIN

##############################################################################
# Analysis of CK_2016-06-23_01 data
##############################################################################

RWD=$RWD_MAIN/CK_2016-06-23_01
ROUT=$RWD/Rout
mkdir -p $ROUT

R CMD BATCH --no-save --no-restore "--args rwd='$RWD' save_cd4_cd8=TRUE outdir_cd4='$RWD_MAIN/CK_2016-06-23_01_CD4/010_cleanfcs' outdir_cd8='$RWD_MAIN/CK_2016-06-23_01_CD8/010_cleanfcs' keep_cd4='CD4' keep_cd8='CD8'" $RCODE/carsten_cytof_pipeline.R $ROUT/pipeline.Rout

tail $ROUT/pipeline.Rout





##############################################################################
# Analysis of CK_2016-06-23_02 data
##############################################################################

RWD=$RWD_MAIN/CK_2016-06-23_02
ROUT=$RWD/Rout
mkdir -p $ROUT

R CMD BATCH --no-save --no-restore "--args rwd='$RWD' save_cd4_cd8=TRUE outdir_cd4='$RWD_MAIN/CK_2016-06-23_02_CD4/010_cleanfcs' outdir_cd8='$RWD_MAIN/CK_2016-06-23_02_CD8/010_cleanfcs' keep_cd4='CD4_T_cells' keep_cd8='CD8_T_cells'" $RCODE/carsten_cytof_pipeline.R $ROUT/pipeline.Rout

tail $ROUT/pipeline.Rout





##############################################################################
# Analysis of CK_2016-06-29_01 data
##############################################################################

RWD=$RWD_MAIN/CK_2016-06-29_01
ROUT=$RWD/Rout
mkdir -p $ROUT

R CMD BATCH --no-save --no-restore "--args rwd='$RWD' save_cd4_cd8=TRUE outdir_cd4='$RWD_MAIN/CK_2016-06-29_01_CD4/010_cleanfcs' outdir_cd8='$RWD_MAIN/CK_2016-06-29_01_CD8/010_cleanfcs' keep_cd4='CD4_T_cells' keep_cd8='CD8_T_cells'" $RCODE/carsten_cytof_pipeline.R $ROUT/pipeline.Rout

tail $ROUT/pipeline.Rout





##############################################################################
# Analysis of CK_2016-06-29_02 data
##############################################################################

RWD=$RWD_MAIN/CK_2016-06-29_02
ROUT=$RWD/Rout
mkdir -p $ROUT

R CMD BATCH --no-save --no-restore "--args rwd='$RWD' save_cd4_cd8=FALSE outdir_cd4='$RWD_MAIN/CK_2016-06-29_02_CD4/010_cleanfcs' outdir_cd8='$RWD_MAIN/CK_2016-06-29_02_CD8/010_cleanfcs' keep_cd4='CD4_T_cells' keep_cd8='CD8_T_cells'" $RCODE/carsten_cytof_pipeline.R $ROUT/pipeline.Rout

tail $ROUT/pipeline.Rout



















