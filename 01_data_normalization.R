
Sys.time()

# Load packages
library(flowCore)
library(gdata)
library(ggplot2)
library(reshape2)

##############################################################################
# Test arguments
##############################################################################

args <- NULL

dir_fcs='CK_2016-06-23_01/010_cleanfcs'
path_panel='CK_panels/panel1.xlsx'
path_metadata='CK_metadata/metadata_23_01.xlsx'
prefix='23_01_'
outdir='CK_2016-06-23_01/010_data'

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


if(!file.exists(outdir)) 
  dir.create(outdir, recursive = TRUE)


# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

# read metadata
md <- read.xls(path_metadata, stringsAsFactors=FALSE)

## colors per sample
color_values <- md$color
names(color_values) <- md$shortname

# ------------------------------------------------------------
# Load fcs files
# ------------------------------------------------------------

# define FCS file names
f <- file.path(dir_fcs, md$filename)
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
cols <- which(fcs_colnames %in% panel$fcs_colname[panel$transform == 1])

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

saveRDS(e_out, file.path(outdir, paste0(prefix, "expr_raw.rds")))



el_out <- data.frame(cell_id = 1:nrow(el), sample_id = samp, el, check.names = FALSE, stringsAsFactors = FALSE)

saveRDS(el_out, file.path(outdir, paste0(prefix, "expr_norm.rds")))





sessionInfo()


