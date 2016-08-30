##############################################################################
## <<02_flowsom.R>>

# BioC 3.3
# Created 27 July 2016
# Updated 24 Aug 2016

##############################################################################
Sys.time()
##############################################################################

# Load packages
library(gdata)
library(FlowSOM)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2)
library(ConsensusClusterPlus)


##############################################################################
# Test arguments
##############################################################################

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
# flowsom_prefix='23_01_pca1_cl20_'
# flowsom_outdir='030_heatmaps'
# path_data='010_data/23_01_expr_raw.rds'
# path_clustering_observables='030_heatmaps/23_01_pca1_clustering_observables.xls'
# rand_seed_consensus=123
# nmetaclusts=20

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

prefix <- flowsom_prefix
outdir <- flowsom_outdir

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
# Load clustering_observables
# ------------------------------------------------------------


clustering_observables <- read.table(path_clustering_observables, header = TRUE, sep = "\t", as.is = TRUE)

clust_observ <- clustering_observables[clustering_observables$clustering_observable, "mass"]

# selected columns for clustering 

scols <- which(fcs_colnames %in% clust_observ)

ef <- as.matrix(e[, scols])



# -------------------------------------
# run FlowSOM and make heatmaps
# -------------------------------------

# SOM
set.seed(rand_seed)
fsom <- FlowSOM::SOM(ef)


# consensus clustering that is reproducible with seed
data <- fsom$codes
k <- nmetaclusts

pdf(file.path(outdir, paste0(prefix, "ConsensusClusterPlus.pdf")), width = 7, height = 7)

results <- ConsensusClusterPlus::ConsensusClusterPlus(t(data),
  maxK = k, reps = 100, pItem = 0.9, pFeature = 1, title = tempdir(),
  plot = NULL, verbose = FALSE, clusterAlg = "hc", distance = "euclidean", seed = rand_seed_consensus)

dev.off()


# get cluster ids
fsom_mc <- results[[k]]$consensusClass
clust <- fsom_mc[fsom$mapping[,1]]


### For plotting tracking plot with my colors
# pdf(file.path(outdir, paste0(prefix, "ConsensusClusterPlus_bar.pdf")), width = 14, height = 4)
# barplot(rep(1, length(results[[k]]$clrs[[1]])), col = results[[k]]$clrs[[1]][results[[k]]$consensusTree$order])
# dev.off()
# 
# results[[k]]$consensusClass[results[[k]]$consensusTree$order]



### Save clustering results

save(fsom, file = file.path(outdir, paste0(prefix, "fsom.rda")))
save(fsom_mc, file = file.path(outdir, paste0(prefix, "fsom_mc.rda")))


clust_out <- data.frame(cluster = clust, cell_id = cell_id, sample_id = samp, stringsAsFactors = FALSE)
write.table(clust_out, file = file.path(outdir, paste0(prefix, "clustering.xls")), row.names=FALSE, quote=FALSE, sep="\t")


# get cluster frequencies
freq_clust <- table(clust)

# make data frame with labels
labels <- data.frame(cluster = sort(unique(clust)), label = sort(unique(clust)), counts = as.numeric(freq_clust))
write.table(labels, file = file.path(outdir, paste0(prefix, "clustering_labels.xls")), row.names=FALSE, quote=FALSE, sep="\t")


sessionInfo()













################################
### 02_flowsom done!
################################