##############################################################################
## <<03_runtsne.R>>

# BioC 3.3
# Created 28 July 2016
# Updated 24 Aug 2016

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

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
# tsne_prefix='23_01_pca1_raw_'
# tsne_outdir='040_tsnemaps'
# path_data='010_data/23_01_expr_raw.rds'
# path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.xlsx'
# path_clustering_observables='030_heatmaps/23_01_pca1_clustering_observables.xls'
# tsne_pmin=1500

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
if(grepl(".txt", path_data)){
  expr <- read.table(path_data, header = TRUE, sep = "\t", as.is = TRUE)
}
if(grepl(".rds", path_data)){
  expr <- readRDS(path_data)
}

cell_id <- expr[, "cell_id"]
samp <- expr[, "sample_id"]
fcs_colnames <- colnames(expr)[!grepl("cell_id|sample_id", colnames(expr))]
e <- expr[, fcs_colnames]


# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

md <- read.xls(path_metadata, stringsAsFactors=FALSE)

# ------------------------------------------------------------
# Load clustering_observables
# ------------------------------------------------------------

clustering_observables <- read.table(path_clustering_observables, header = TRUE, sep = "\t", as.is = TRUE)

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


et_sub <- et[cells2keep, ]



### Run tSNE
set.seed(rand_seed)
rtsne_out <- Rtsne(et_sub, perplexity = 30, pca = FALSE, max_iter = 1000, verbose = TRUE)


# Save rtsne results

rtsne_data <- data.frame(cell_index = cells2keep, sample_name = samp[cells2keep], et_sub)

write.table(rtsne_data, file.path(outdir, paste0(prefix, "rtsne_data.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


save(rtsne_out, file = file.path(outdir, paste0(prefix, "rtsne_out.rda")))










sessionInfo()























################################
### 03_runtsne done!
################################