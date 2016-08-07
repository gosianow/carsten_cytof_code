##############################################################################
## <<02_fcs_saving.R>>

# BioC 3.3
# Created 7 Aug 2016

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

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_02_CD4_merging2'
# save_prefix='pnlCD4_'
# save_dir='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_02_CD4_merging2/070_dumpfcs'
# path_panel='panel_CD4.xlsx'


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

dir.create(save_dir, recursive = TRUE)


# ------------------------------------------------------------
# define directories
# ------------------------------------------------------------

fcsDir <- "010_cleanfcs"; if( !file.exists(fcsDir) ) dir.create(fcsDir)
dmpDir <- "070_dumpfcs"; if( !file.exists(dmpDir) ) dir.create(dmpDir)


# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

# read metadata
md <- read.xls("metadata.xlsx",stringsAsFactors=FALSE)


# define FCS file names
f <- file.path(fcsDir, md$filename)
names(f) <- md$shortname


# read raw FCS files in
fcs <- lapply(f, read.FCS)


# read panel, pick which columns to use
panel <- read.xls(path_panel,stringsAsFactors=FALSE)

# get isotope mass of columns in fcs files.. to match against the panel
panel_mass <- as.numeric(gsub("[[:alpha:]]", "", colnames(fcs[[1]])))

# cols - get fcs columns that are in the panel with transform = 1
cols <- which(panel_mass %in% panel$Isotope[panel$transform==1])


# Antigen - get the antigen name
m <- match(panel_mass, panel$Isotope)
fcs_panel <- data.frame(colnames = colnames(fcs[[1]]), Isotope = panel_mass, cols = panel_mass %in% panel$Isotope[panel$transform==1], Antigen = panel$Antigen[m], stringsAsFactors = FALSE)
fcs_panel$Antigen[is.na(fcs_panel$Antigen)] <- ""


# arc-sin-h the columns specific 
fcsT <- lapply(fcs, function(u) {
  e <- exprs(u)
  e[,cols] <- asinh( e[,cols] / 5 )
  exprs(u) <- e
  u
})


# ------------------------------------------------------------
# Normalization to 0-1
# ------------------------------------------------------------


# re-extract the data as list of tables
es <- lapply(fcsT, function(u) {
  exprs(u)[,cols]
})

e <- do.call("rbind", es)


# normalize to 0-1
el <- e
rng <- apply(el,2,quantile,p=c(.01,.99))
for(i in 1:ncol(el)) {
  el[,i] <- (el[,i]-rng[1,i])/(rng[2,i]-rng[1,i])
}
el[el<0] <- 0
el[el>1] <- 1


## replace expression in fcsT with 0-1 normalized expr.
samp <- rep( names(fcsT), sapply(fcsT, nrow) )
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
  
  fn <- file.path(save_dir, paste0(save_prefix, "arcsineh_", basename(f[i])))
  fcs_out <- fcsT[[i]][, cols]
  write.FCS(fcs_out, fn)
  
})


### arcsineh + 0-1 normalized
dummy <- lapply(seq(length(fcsT01)), function(i){
  # i = 1
  
  fn <- file.path(save_dir, paste0(save_prefix, "arcsineh01_", basename(f[i])))
  fcs_out <- fcsT01[[i]][, cols]
  write.FCS(fcs_out, fn)
  
})










################################
### 02_fcs_saving done!
################################