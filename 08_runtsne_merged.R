##############################################################################
## <<08_runtsne_merged.R>>

# BioC 3.3
# Created 16 Oct 2016
# Updated 16 Oct 2016

##############################################################################
Sys.time()
##############################################################################

# Load packages
library(gdata)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2)
library(Rtsne)

##############################################################################
# Test arguments
##############################################################################

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-merged_23_29/01'
# tsne_prefix='23m6_29m4_'
# tsne_outdir='08_tsnemaps_merged'
# 
# path_data=c('/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01/010_data/23_01_expr_raw.rds','/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-29_01/010_data/29_01_expr_raw.rds')
# path_metadata=c('/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.xlsx','/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_29_01.xlsx')
# path_clustering_observables=c('/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01/030_heatmaps/23_01_pca1_clustering_observables.xls','/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-29_01/030_heatmaps/29_01_pca1_clustering_observables.xls')
# tsne_pmin=1500
# data_name=c('data23','data29')

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
rand_seed <- 1234

prefix <- tsne_prefix
outdir <- tsne_outdir

if(!file.exists(outdir)) 
  dir.create(outdir)


# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

# read expression data

expr <- lapply(1:length(data_name), function(i){
  
  if(grepl(".txt", path_data[i])){
    expr <- read.table(path_data[i], header = TRUE, sep = "\t", as.is = TRUE)
  }
  if(grepl(".rds", path_data[i])){
    expr <- readRDS(path_data[i])
  }
  
  expr
})

expr <- rbind.fill(expr)

### cell_id is a paste or cell_id and sample_id because cell ids are the same for data 23 and 29 
cell_id <- paste0(expr[, "cell_id"], "-", expr[, "sample_id"])
samp <- expr[, "sample_id"]
fcs_colnames <- colnames(expr)[!grepl("cell_id|sample_id", colnames(expr))]
e <- expr[, fcs_colnames]


# ------------------------------------------------------------
# Load metadata - is not really needed
# ------------------------------------------------------------

md <- lapply(1:length(data_name), function(i){
  
  path <- path_metadata[i]
  md <- read.xls(path, stringsAsFactors=FALSE)
  md$data <- data_name[i]
  md
  
})

md <- rbind.fill(md)

# ------------------------------------------------------------
# Load clustering_observables
# ------------------------------------------------------------
# Use a union of observalbles used in each dataset


clustering_observables <- lapply(1:length(data_name), function(i){
  
  clustering_observables <- read.table(path_clustering_observables[i], header = TRUE, sep = "\t", as.is = TRUE)
  return(clustering_observables)
  
})

clustering_observables <- Reduce(function(...) merge(..., by = c("mass", "marker"), all=TRUE, sort = FALSE), clustering_observables)

clustering_observables$clustering_observable <- apply(clustering_observables[, grep("clustering_observable", colnames(clustering_observables))], 1, any)

clust_observ <- clustering_observables[clustering_observables$clustering_observable, "mass"]

# -------------------------------------

# Indeces of observables used for clustering 
scols <- which(fcs_colnames %in% clust_observ)


# ---------------------------------------
# tSNE analyses
# ---------------------------------------

et <- e[, scols]

### find duplicates
dups <- duplicated(et)  
w <- which(!dups)


### Data subsampling
# create indices by sample
inds <- split(1:length(samp), samp) 

# per-sample, how many cells to downsample
ts <- table(samp)
ns <- pmin(ts, tsne_pmin)  

# get subsampled indices
subs <- mapply(function(u,v) {
  set.seed(rand_seed)
  s <- sample(u, ns[v], replace = FALSE)
  intersect(s,w)
}, inds, names(inds))


cells2keep <- c(unlist(subs))



### Run tSNE

et_sub <- et[cells2keep, ]

set.seed(rand_seed)
rtsne_out <- Rtsne(et_sub, pca = FALSE, verbose = TRUE)


# Save rtsne results

rtsne_data <- data.frame(cell_index = cell_id[cells2keep], sample_name = samp[cells2keep], et_sub)

write.table(rtsne_data, file.path(outdir, paste0(prefix, "rtsne_data.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


save(rtsne_out, file = file.path(outdir, paste0(prefix, "rtsne_out.rda")))








sessionInfo()























################################
### 08_runtsne_merged.R done!
################################