##############################################################################

# BioC 3.3
# Created 24 Aug 2016
# Updated 8 Mar 2018


##############################################################################
Sys.time()
##############################################################################

# Load packages
library(flowCore)
library(gdata)
library(tools)
library(xlsx)

##############################################################################
# Test arguments
##############################################################################


# rwd='/Users/nowickm2/Desktop/tmp/carsten_cytof/PD1_project/Flow_repository/ck_repo_data3/CK_fcs_files/CK_2016-06-23_01'
# path_metadata='/Users/nowickm2/Desktop/tmp/carsten_cytof/PD1_project/CK_metadata/metadata_23_01.xlsx'
# path_panel_export='/Users/nowickm2/Desktop/tmp/carsten_cytof/PD1_project/CK_helpfiles/panel1_export.xlsx'
# panels_path="/Users/nowickm2/Desktop/tmp/carsten_cytof/PD1_project/CK_panels"
# panels_prefix="panel1"
# data_name='Data23'
# panel_name='Panel1'


rwd='../../carsten_cytof/PD1_project/Flow_repository/ck_repo_data3/CK_fcs_files/CK_2016-06-29_01'
path_metadata='../../carsten_cytof/PD1_project/CK_metadata/metadata_29_01.xlsx'
path_panel_export='../../carsten_cytof/PD1_project/CK_helpfiles/panel1_export.xlsx'
panels_path="../../carsten_cytof/PD1_project/CK_panels"
panels_prefix="panel1"
data_name='Data29'
panel_name='Panel1'

##############################################################################
# Read in the arguments
##############################################################################

rm(list = ls())

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

cat(paste0(args, collapse = "\n"), fill = TRUE)


##############################################################################

# ------------------------------------------------------------
# define directories
# ------------------------------------------------------------

fcsDir <- paste0(rwd, "/010_cleanfcs")
fcsDir_new <- "../../carsten_cytof/PD1_project/Flow_repository/ck_repo_data3/CK_fcs_files"
outdir_metadata <- "../../carsten_cytof/PD1_project/Flow_repository/ck_repo_data3/CK_metadata"
outdir_panels <- "../../carsten_cytof/PD1_project/Flow_repository/ck_repo_data3/CK_panels"


# ------------------------------------------------------------


if(!file.exists(outdir_metadata)) 
  dir.create(outdir_metadata, recursive = TRUE)

    
if(!file.exists(outdir_panels)) 
  dir.create(outdir_panels, recursive = TRUE)


# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

# read metadata
md <- read.xls(path_metadata, stringsAsFactors=FALSE)

# ------------------------------------------------------------
# Load panel
# ------------------------------------------------------------

# read panel to export
panel <- read.xls(path_panel_export, stringsAsFactors=FALSE)

export_to_fcs <- panel$fcs_colname[as.logical(panel$export_to_fcs)]

# ------------------------------------------------------------
# Export fcs files
# ------------------------------------------------------------


new_names <- paste0(data_name, "_", panel_name, "_", md$shortname, "_", gsub("patient", "Patient", md$patient_id), ".fcs")

  
for(i in 1:nrow(md)){
  # i = 1
  
  f <- file.path(fcsDir, md$filename[i])
  
  fcs <- read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
  
  stopifnot(all(export_to_fcs %in% colnames(fcs)))

  export_to_fcs2 <- colnames(fcs)
  
  pData(fcs@parameters)
  
  # Change the file name 
  print(fcs@description$FIL)
  print(fcs@description$GUID)
  print(fcs@description$ORIGINALGUID)
  
  fcs@description$FIL <- new_names[i]
  fcs@description$GUID <- new_names[i]
  fcs@description$ORIGINALGUID <- new_names[i]
  
  
  write.FCS(fcs, file.path(fcsDir_new, new_names[i]))

  fcs_new2 <- read.FCS(file.path(fcsDir_new, new_names[i]), transformation = FALSE, truncate_max_range = FALSE)

  print(fcs_new2@description$FIL)
  print(fcs_new2@description$GUID)
  print(fcs_new2@description$ORIGINALGUID)
  
  all(exprs(fcs) == exprs(fcs_new2))
  colnames(fcs) == colnames(fcs_new2)

}


# ------------------------------------------------------------
# Change metadata
# ------------------------------------------------------------

md_new <- data.frame(filename_old = md$filename, md, stringsAsFactors = FALSE)

md_new$filename <- new_names

write.xlsx(md_new, file.path(outdir_metadata, basename(path_metadata)), row.names = FALSE)



# ------------------------------------------------------------
# Change panels
# ------------------------------------------------------------


panel_files <- list.files(path = panels_path, pattern = paste0(panels_prefix, ".*xlsx"), full.names = TRUE)

for(i in 1:length(panel_files)){
  # i = 1
  panel <- read.xls(panel_files[i], stringsAsFactors=FALSE)
  
  out <- panel[panel$fcs_colname %in% export_to_fcs2, , drop = FALSE]
  
  write.xlsx(out, file.path(outdir_panels, basename(panel_files[i])), row.names = FALSE)
  
}













sessionInfo()

