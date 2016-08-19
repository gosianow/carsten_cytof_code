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

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
# path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.xlsx'
# tsne_prefix='23_01_pca1_'
# path_clustering_observables='23_01_pca1_clustering_observables.xls'
# tsne_pmin=1500


##############################################################################
# Read in the arguments
##############################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(args)
print(rwd)
print(path_clustering_observables)
print(tsne_prefix)
print(tsne_pmin)

##############################################################################

setwd(rwd)
rand_seed <- 1234

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
md <- read.xls(path_metadata, stringsAsFactors=FALSE)

# define FCS file names
f <- file.path(fcsDir, md$filename)
names(f) <- md$shortname

# read raw FCS files in
fcs <- lapply(f, read.FCS)


# ------------------------------------------------------------
# Load more data
# ------------------------------------------------------------

if(!grepl("/", path_clustering_observables)){
  clust_observ <- read.table(file.path(hmDir, path_clustering_observables), header = TRUE, sep = "\t", as.is = TRUE)
}else{
  clust_observ <- read.table(path_clustering_observables, header = TRUE, sep = "\t", as.is = TRUE)
}

clust_observ <- clust_observ[, 1]

# -------------------------------------

# selected columns for clustering  

scols <- which(colnames(fcs[[1]]) %in% clust_observ)

# arc-sin-h the columns specific 
fcsT <- lapply(fcs, function(u) {
  e <- exprs(u)
  e[ , scols] <- asinh( e[ , scols] / 5 )
  exprs(u) <- e
  u
})


# ---------------------------------------
# tSNE analyses
# ---------------------------------------


# re-extract the data as list of tables
es <- lapply(fcsT, function(u) {
  exprs(u)[ , scols]
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
  set.seed(rand_seed)
  s <- sample(u, ns[v], replace = FALSE)
  intersect(s,w)
}, inds, names(inds))


inds2keep <- c(unlist(subs))





### Run tSNE on normalized data

# el_sub <- el[inds2keep, ]
# 
# set.seed(rand_seed)
# rtsne_out <- Rtsne(el_sub, pca = FALSE, verbose = TRUE)
# 
# 
# # Save rtsne results
# 
# rtsne_data <- data.frame(cell_index = inds2keep, sample_name = samp[inds2keep], el_sub)
# 
# write.table(rtsne_data, file.path(sneDir, paste0(prefix, "rtsne_data_norm.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
# 
# 
# save(rtsne_out, file = file.path(sneDir, paste0(prefix, "rtsne_out_norm.rda")))






### Run tSNE on raw data

e_sub <- e[inds2keep, ]

set.seed(rand_seed)
rtsne_out <- Rtsne(e_sub, pca = FALSE, verbose = TRUE)


# Save rtsne results

rtsne_data <- data.frame(cell_index = inds2keep, sample_name = samp[inds2keep], e_sub)

write.table(rtsne_data, file.path(sneDir, paste0(prefix, "rtsne_data_raw.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


save(rtsne_out, file = file.path(sneDir, paste0(prefix, "rtsne_out_raw.rda")))








sessionInfo()























################################
### 03_runtsne done!
################################