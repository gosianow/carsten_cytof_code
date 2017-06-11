

Sys.time()

# Load packages
library(FlowSOM)
library(ConsensusClusterPlus)

##############################################################################
# Test arguments
##############################################################################

prefix='23_01_pca1_cl20_'
outdir='CK_2016-06-23_01/030_heatmaps'
path_data='CK_2016-06-23_01/010_data/23_01_expr_raw.rds'
path_clustering_observables='CK_2016-06-23_01/030_heatmaps/23_01_pca1_clustering_observables.xls'
rand_seed_consensus=123
nmetaclusts=20

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


rand_seed <- 1234
linkage <- "average"


if(!file.exists(outdir)) 
  dir.create(outdir, recursive = TRUE)


# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

expr <- readRDS(path_data)

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


# --------------------------------------------------------------------------
# run FlowSOM 
# --------------------------------------------------------------------------


# -------------------------------------
# SOM
# -------------------------------------

set.seed(rand_seed)
fsom <- FlowSOM::SOM(ef, rlen = 10)


# -------------------------------------
# metaClustering_consensus 
# -------------------------------------
# consensus clustering that is reproducible with seed

### Sometimes not all the codes are used in mapping
codes <- fsom$codes
k <- nmetaclusts

pdf(file.path(outdir, paste0(prefix, "ConsensusClusterPlus.pdf")), width = 7, height = 7)

fccp <- ConsensusClusterPlus::ConsensusClusterPlus(t(codes),
  maxK = k, reps = 100, pItem = 0.9, pFeature = 1, title = tempdir(),
  plot = NULL, verbose = FALSE, clusterAlg = "hc", innerLinkage = linkage, finalLinkage = "average", distance = "euclidean", seed = rand_seed_consensus)

dev.off()


### Get cluster ids
fsom_mc <- fccp[[k]]$consensusClass

clust <- fsom_mc[fsom$mapping[,1]]


# -------------------------------------
# Save clustering results
# -------------------------------------

saveRDS(fsom, file = file.path(outdir, paste0(prefix, "fsom.rds")))
saveRDS(fccp, file = file.path(outdir, paste0(prefix, "fccp.rds")))


clust_out <- data.frame(cluster = clust, cell_id = cell_id, sample_id = samp, stringsAsFactors = FALSE)
write.table(clust_out, file = file.path(outdir, paste0(prefix, "clustering.xls")), row.names=FALSE, quote=FALSE, sep="\t")


# -------------------------------------
# Save clustering labels
# -------------------------------------

# get cluster frequencies
freq_clust <- table(clust)

# make data frame with labels
labels <- data.frame(cluster = sort(unique(clust)), label = sort(unique(clust)), counts = as.numeric(freq_clust))
labels$proportions <- round(labels$counts/sum(labels$counts) * 100, 2)


write.table(labels, file = file.path(outdir, paste0(prefix, "clustering_labels.xls")), row.names=FALSE, quote=FALSE, sep="\t")




























sessionInfo()



