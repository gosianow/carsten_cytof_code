##############################################################################
## <<01_data_normalization.R>>

# BioC 3.3
# Created 24 Aug 2016


##############################################################################
Sys.time()
##############################################################################

# Load packages
library(flowCore)
library(gdata)
library(Repitools)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2)


##############################################################################
# Test arguments
##############################################################################

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
# data_prefix='23_01_'
# data_outdir='010_data'
# path_panel='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_panels/panel1.xlsx'
# path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.xlsx'

##############################################################################
# Read in the arguments
##############################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(args)

##############################################################################

setwd(rwd)

prefix <- data_prefix

# ------------------------------------------------------------
# define directories
# ------------------------------------------------------------

fcsDir <- "010_cleanfcs"
dataDir <- data_outdir
if( !file.exists(dataDir) ) dir.create(dataDir)


# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

# read metadata
md <- read.xls(path_metadata, stringsAsFactors=FALSE)

# ------------------------------------------------------------
# Load fcs files
# ------------------------------------------------------------

# define FCS file names
f <- file.path(fcsDir, md$filename)
names(f) <- md$shortname


# read raw FCS files in
fcs <- lapply(f, read.FCS)

fcs_colnames <- colnames(fcs[[1]])

## Create sample info
samp <- rep(names(fcs), sapply(fcs, nrow))

# ------------------------------------------------------------
# Load panel
# ------------------------------------------------------------

# read panel, pick which columns to use
panel <- read.xls(path_panel, stringsAsFactors=FALSE)

if(!all(c("fcs_colname", "Isotope", "Antigen", "transform") %in% colnames(panel)))
  stop("Wrong columns in panel!!!")


# cols - get fcs columns that are in the panel with transform = 1
cols <- which(fcs_colnames %in% panel$fcs_colname[panel$transform==1])

# get the isotope and antigen for fcs markers
m <- match(fcs_colnames, panel$fcs_colname)

fcs_panel <- data.frame(fcs_colname = fcs_colnames, Isotope = panel$Isotope[m], Antigen = panel$Antigen[m], stringsAsFactors = FALSE)



# --------------------------------------------------------------------------
# Get the marker expression
# --------------------------------------------------------------------------


# arc-sin-h the columns specific 
fcsT <- lapply(fcs, function(u) {
  e <- exprs(u)
  e[,cols] <- asinh( e[,cols] / 5 )
  exprs(u) <- e
  u
})


es <- lapply(fcsT, function(u) {
  exprs(u)[, cols]
})

e <- do.call("rbind",es)

# normalize to 0-1
el <- e
rng <- apply(el,2,quantile,p=c(.01,.99))
for(i in 1:ncol(el)) {
  el[,i] <- (el[,i]-rng[1,i])/(rng[2,i]-rng[1,i])
}
el[el<0] <- 0
el[el>1] <- 1


### Save the data

e_out <- data.frame(cell_id = 1:nrow(e), sample_id = samp, e, check.names = FALSE, stringsAsFactors = FALSE)

saveRDS(e_out, file.path(dataDir, paste0(prefix, "expr_raw.rds")))



el_out <- data.frame(cell_id = 1:nrow(el), sample_id = samp, el, check.names = FALSE, stringsAsFactors = FALSE)

saveRDS(el_out, file.path(dataDir, paste0(prefix, "expr_norm.rds")))







sessionInfo()













################################
### 01_data_normalization.R done!
################################