##############################################################################
## <<02_fcs_saving.R>>

# BioC 3.3
# Created 7 Aug 2016
# Updated 14 Oct 2016

##############################################################################
Sys.time()
##############################################################################

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

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-29_02_CD4_merging'
save_prefix='29CD4_02CD4_'
save_outdir='060_dumpfcs'
path_panel='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_panels/panel2CD4.xlsx'
path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_29_02.xlsx'

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

setwd(rwd)


# ------------------------------------------------------------
# define directories
# ------------------------------------------------------------

fcsDir <- "010_cleanfcs"
outdir <- save_outdir
if( !file.exists(outdir) ) dir.create(outdir)


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
fcs_colnames

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




# fcsT01 = fcsT with replaced expression which is normalized to 0-1
samp <- factor(samp, levels = names(fcsT))

el_spl <- split(data.frame(el), samp)


fcsT01 <- lapply(seq(length(fcsT)), function(i){
  # i = 1
  u <- fcsT[[i]]
  exprs(u)[,cols] <- as.matrix(el_spl[[i]])
  u
})



# ------------------------------------------------------------
# Save the new fcs files
# ------------------------------------------------------------

### arcsineh normalized
dummy <- lapply(seq(length(fcsT)), function(i){
  # i = 1
  
  fn <- file.path(outdir, paste0(save_prefix, "arcsineh_", basename(f[i])))
  fcs_out <- fcsT[[i]]
  write.FCS(fcs_out, fn)
  
})


### arcsineh + 0-1 normalized
dummy <- lapply(seq(length(fcsT01)), function(i){
  # i = 1
  
  fn <- file.path(outdir, paste0(save_prefix, "arcsineh01_", basename(f[i])))
  fcs_out <- fcsT01[[i]]
  
  write.FCS(fcs_out, fn)
  
})










################################
### 02_fcs_saving done!
################################