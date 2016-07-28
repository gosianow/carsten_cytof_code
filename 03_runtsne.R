##############################################################################
## <<03_runtsne.R>>

# BioC 3.3
# Created 28 July 2016
# Updated 

##############################################################################
Sys.time()
##############################################################################

# Load packages
library(flowCore)
library(gdata)
library(FlowSOM)
library(Repitools)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2)
library(coop)
library(pheatmap)
library(RColorBrewer)
library(ggdendro)
library(Rtsne)


##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
path_clustering_observables='pca1_cl20_clustering_observables.xls'
tsne_prefix='pca1_cl20_'
tsne_pmin=1500

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

prefix <- tsne_prefix


# ------------------------------------------------------------
# define directories
# ------------------------------------------------------------

fcsDir <- "010_cleanfcs"; if( !file.exists(fcsDir) ) dir.create(fcsDir)
pcaDir <- "020_pcascores"; if( !file.exists(pcaDir) ) dir.create(pcaDir)
hmDir <- "030_heatmaps"; if( !file.exists(hmDir) ) dir.create(hmDir)
sneDir <- "040_tsnemaps"; if( !file.exists(sneDir) ) dir.create(sneDir)


# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

# read metadata
md <- read.xls("metadata.xlsx",stringsAsFactors=FALSE)

# read panel, pick which columns to use
panel <- read.xls("panel.xlsx",stringsAsFactors=FALSE)


# define FCS file names
f <- file.path(fcsDir, md$filename)
names(f) <- md$shortname


# read raw FCS files in
fcs <- lapply(f, read.FCS)


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
# Load more data
# ------------------------------------------------------------

clust_observ <- read.table(file.path(hmDir, path_clustering_observables), header = TRUE, sep = "\t", as.is = TRUE)
clust_observ <- clust_observ[, 1]

# ------------------------------------------------------------

### Indeces of observables used for clustering 

scols <- which(fcs_panel$colnames %in% clust_observ)


# ---------------------------------------
# tSNE analyses
# ---------------------------------------


# re-extract the data as list of tables
es <- lapply(fcsT, function(u) {
  exprs(u)[,scols]
})

e <- do.call("rbind",es)
colnames(e) <- fcs_panel$Antigen[scols]



# normalize to 0-1
el <- e
rng <- apply(el,2,quantile,p=c(.01,.99))
for(i in 1:ncol(el)) {
  el[,i] <- (el[,i]-rng[1,i])/(rng[2,i]-rng[1,i])
}
el[el<0] <- 0
el[el>1] <- 1



# find duplicates
dups <- duplicated(el)  
w <- which(!dups)


### Data subsampling
samp <- rep( names(fcsT), sapply(fcsT,nrow) )

# create indices by sample
inds <- split(1:nrow(el), samp) 

# per-sample, how many cells to downsample
ts <- table(samp)
ns <- pmin(ts, tsne_pmin)  


# get subsampled indices
subs <- mapply(function(u,v) {
  set.seed(1234)
  s <- sample(u, ns[v], replace = FALSE)
  intersect(s,w)
}, inds, names(inds))


inds2keep <- c(unlist(subs))

el_sub <- el[inds2keep, ]



### Run tSNE
set.seed(1234)
rtsne_out <- Rtsne(el_sub, pca = FALSE, verbose = TRUE)

# load(file.path(sneDir, paste0(prefix, "rtsne_out.rda")))


### Save rtsne results

rtsne_data <- data.frame(cell_index = inds2keep, sample_name = samp[inds2keep], el_sub)

write.table(rtsne_data, file.path(sneDir, paste0(prefix, "rtsne_data.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


save(rtsne_out, file = file.path(sneDir, paste0(prefix, "rtsne_out.rda")))












sessionInfo()










