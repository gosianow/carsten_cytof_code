#!/bin/bash

###############################################################################################################
## Define paths to software and reference files

RWD_MAIN=../../carsten_cytof/PD1_project
RCODE=../../carsten_cytof_code
METADATA=$RWD_MAIN/CK_metadata
PANELS=$RWD_MAIN/CK_panels


###############################################################################################################


DATA=23
PANEL=1
data_dir="CK_2016-06-23_01"
file_metadata="metadata_23_01.xlsx"
file_panel="panel1.xlsx"
file_panel_export="panel1_export.xlsx"
panels_prefix="panel${PANEL}"
data_name="Data${DATA}"
panel_name="Panel${PANEL}"


RWD=$RWD_MAIN/Flow_repository/ck_repo_data3/CK_fcs_files/${data_dir}
ROUT=$RWD/Rout
mkdir -p $ROUT

R CMD BATCH --no-save --no-restore "--args rwd='${RWD}' path_metadata='${METADATA}/${file_metadata}' path_panel_export='$RWD_MAIN/CK_helpfiles/${file_panel_export}' panels_path='${PANELS}' panels_prefix='${panels_prefix}' data_name='${data_name}' panel_name='${panel_name}'" $RCODE/quick_analysis/rename_fcs.R $ROUT/rename_fcs.Rout
tail $ROUT/rename_fcs.Rout


DATA=23
PANEL=2
data_dir="CK_2016-06-23_02"
file_metadata="metadata_23_02.xlsx"
file_panel="panel2.xlsx"
file_panel_export="panel2_export.xlsx"
panels_prefix="panel${PANEL}"
data_name="Data${DATA}"
panel_name="Panel${PANEL}"


RWD=$RWD_MAIN/Flow_repository/ck_repo_data3/CK_fcs_files/${data_dir}
ROUT=$RWD/Rout
mkdir -p $ROUT

R CMD BATCH --no-save --no-restore "--args rwd='${RWD}' path_metadata='${METADATA}/${file_metadata}' path_panel_export='$RWD_MAIN/CK_helpfiles/${file_panel_export}' panels_path='${PANELS}' panels_prefix='${panels_prefix}' data_name='${data_name}' panel_name='${panel_name}'" $RCODE/quick_analysis/rename_fcs.R $ROUT/rename_fcs.Rout
tail $ROUT/rename_fcs.Rout


DATA=23
PANEL=3
data_dir="CK_2016-06-23_03"
file_metadata="metadata_23_03.xlsx"
file_panel="panel3.xlsx"
file_panel_export="panel3_export.xlsx"
panels_prefix="panel${PANEL}"
data_name="Data${DATA}"
panel_name="Panel${PANEL}"


RWD=$RWD_MAIN/Flow_repository/ck_repo_data3/CK_fcs_files/${data_dir}
ROUT=$RWD/Rout
mkdir -p $ROUT

R CMD BATCH --no-save --no-restore "--args rwd='${RWD}' path_metadata='${METADATA}/${file_metadata}' path_panel_export='$RWD_MAIN/CK_helpfiles/${file_panel_export}' panels_path='${PANELS}' panels_prefix='${panels_prefix}' data_name='${data_name}' panel_name='${panel_name}'" $RCODE/quick_analysis/rename_fcs.R $ROUT/rename_fcs.Rout
tail $ROUT/rename_fcs.Rout



DATA=29
PANEL=1
data_dir="CK_2016-06-29_01"
file_metadata="metadata_29_01.xlsx"
file_panel="panel1.xlsx"
file_panel_export="panel1_export.xlsx"
panels_prefix="panel${PANEL}"
data_name="Data${DATA}"
panel_name="Panel${PANEL}"


RWD=$RWD_MAIN/Flow_repository/ck_repo_data3/CK_fcs_files/${data_dir}
ROUT=$RWD/Rout
mkdir -p $ROUT

R CMD BATCH --no-save --no-restore "--args rwd='${RWD}' path_metadata='${METADATA}/${file_metadata}' path_panel_export='$RWD_MAIN/CK_helpfiles/${file_panel_export}' panels_path='${PANELS}' panels_prefix='${panels_prefix}' data_name='${data_name}' panel_name='${panel_name}'" $RCODE/quick_analysis/rename_fcs.R $ROUT/rename_fcs.Rout
tail $ROUT/rename_fcs.Rout


DATA=29
PANEL=2
data_dir="CK_2016-06-29_02"
file_metadata="metadata_29_02.xlsx"
file_panel="panel2.xlsx"
file_panel_export="panel2_export.xlsx"
panels_prefix="panel${PANEL}"
data_name="Data${DATA}"
panel_name="Panel${PANEL}"


RWD=$RWD_MAIN/Flow_repository/ck_repo_data3/CK_fcs_files/${data_dir}
ROUT=$RWD/Rout
mkdir -p $ROUT

R CMD BATCH --no-save --no-restore "--args rwd='${RWD}' path_metadata='${METADATA}/${file_metadata}' path_panel_export='$RWD_MAIN/CK_helpfiles/${file_panel_export}' panels_path='${PANELS}' panels_prefix='${panels_prefix}' data_name='${data_name}' panel_name='${panel_name}'" $RCODE/quick_analysis/rename_fcs.R $ROUT/rename_fcs.Rout
tail $ROUT/rename_fcs.Rout


DATA=29
PANEL=3
data_dir="CK_2016-06-29_03"
file_metadata="metadata_29_03.xlsx"
file_panel="panel3.xlsx"
file_panel_export="panel3_export.xlsx"
panels_prefix="panel${PANEL}"
data_name="Data${DATA}"
panel_name="Panel${PANEL}"


RWD=$RWD_MAIN/Flow_repository/ck_repo_data3/CK_fcs_files/${data_dir}
ROUT=$RWD/Rout
mkdir -p $ROUT

R CMD BATCH --no-save --no-restore "--args rwd='${RWD}' path_metadata='${METADATA}/${file_metadata}' path_panel_export='$RWD_MAIN/CK_helpfiles/${file_panel_export}' panels_path='${PANELS}' panels_prefix='${panels_prefix}' data_name='${data_name}' panel_name='${panel_name}'" $RCODE/quick_analysis/rename_fcs.R $ROUT/rename_fcs.Rout
tail $ROUT/rename_fcs.Rout

















#
