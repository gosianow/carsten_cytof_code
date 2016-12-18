##############################################################################
## <<03_run_dimension_reduction.R>>

# BioC 3.3
# Created 2 Nov 2016
# Updated 2 Nov 2016

##############################################################################
Sys.time()
##############################################################################

# Load packages
library(gdata)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2)
library(cytofkit)

##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
dr_prefix='23_01_pca1_raw_'
dr_outdir='040_dimension_reduction'
path_data='010_data/23_01_expr_raw.rds'
path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.xlsx'
path_clustering_observables='030_heatmaps/23_01_pca1_clustering_observables.xls'
pmin=1500

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
rand_seed <- 1234

prefix <- dr_prefix
outdir <- dr_outdir

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
ns <- pmin(ts, pmin)  

# get subsampled indices
subs <- mapply(function(u,v) {
  set.seed(rand_seed)
  s <- sample(u, ns[v], replace = FALSE)
  intersect(s,w)
}, inds, names(inds))


cells2keep <- c(unlist(subs))


et_sub <- et[cells2keep, ]

dim(et_sub)


dr_isomap <- cytofkit::cytof_dimReduction(et_sub, method="isomap")

dr_diffusion <- cytofkit::cytof_dimReduction(et_sub, method="diffusionmap")

dr_pca <- cytofkit::cytof_dimReduction(et_sub, method="pca")




# Save the results

dr_isomap <- data.frame(dr_isomap)
colnames(dr_isomap) <- c("dim1", "dim2")

isomap_data <- data.frame(cell_index = cells2keep, sample_name = samp[cells2keep], dr_isomap)

write.table(isomap_data, file.path(outdir, paste0(prefix, "isomap_data.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


dr_diffusion <- data.frame(dr_diffusion)
colnames(dr_diffusion) <- c("dim1", "dim2")

diffusion_data <- data.frame(cell_index = cells2keep, sample_name = samp[cells2keep], dr_diffusion)

write.table(diffusion_data, file.path(outdir, paste0(prefix, "diffusion_data.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


dr_pca <- data.frame(dr_pca)
colnames(dr_pca) <- c("dim1", "dim2")

pca_data <- data.frame(cell_index = cells2keep, sample_name = samp[cells2keep], dr_pca)

write.table(pca_data, file.path(outdir, paste0(prefix, "pca_data.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



























################################
### 03_run_dimension_reduction.R done!
################################