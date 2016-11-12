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
library(igraph)
library(RColorBrewer)
library(pheatmap)
library(cytofkit)


##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
flowsom_prefix='23_01_pca1_cl20_'
flowsom_outdir='030_heatmaps'
path_data='010_data/23_01_expr_raw.rds'
path_clustering_observables='030_heatmaps/23_01_pca1_clustering_observables.xls'
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

setwd(rwd)
rand_seed <- 1234

prefix <- flowsom_prefix
outdir <- flowsom_outdir

if(!file.exists(outdir)) 
  dir.create(outdir, recursive = TRUE)

if(all(!grepl("linkage=", args))){
  linkage <- "average"
}


# ------------------------------------------------------------
# Colors for clusters
# ------------------------------------------------------------

# ggplot palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=60 , c=100)[1:n]
}

# ------------------------------ 
# color blind palette

colors_muted <- c("#DC050C", "#E8601C", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6", "#F1932D", "#F6C141", "#F7EE55", "#4EB265", "#90C987", "#CAEDAB")
color_ramp <- c(colors_muted, gg_color_hue(max(1, nmetaclusts - length(colors_muted))))

colors_clusters <- color_ramp[1:nmetaclusts]
names(colors_clusters) <- 1:nmetaclusts


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
data <- fsom$codes
k <- nmetaclusts

pdf(file.path(outdir, paste0(prefix, "ConsensusClusterPlus.pdf")), width = 7, height = 7)

results <- ConsensusClusterPlus::ConsensusClusterPlus(t(data),
  maxK = k, reps = 100, pItem = 0.9, pFeature = 1, title = tempdir(),
  plot = NULL, verbose = FALSE, clusterAlg = "hc", innerLinkage = linkage, finalLinkage = "average", distance = "euclidean", seed = rand_seed_consensus)

dev.off()


# get cluster ids
fsom_mc <- results[[k]]$consensusClass
fsom_mc_tree <- results[[k]]$consensusTree

clust <- fsom_mc[fsom$mapping[,1]]


### Save clustering results


save(fsom, file = file.path(outdir, paste0(prefix, "fsom.rda")))
save(fsom_mc, file = file.path(outdir, paste0(prefix, "fsom_mc.rda")))


clust_out <- data.frame(cluster = clust, cell_id = cell_id, sample_id = samp, stringsAsFactors = FALSE)
write.table(clust_out, file = file.path(outdir, paste0(prefix, "clustering.xls")), row.names=FALSE, quote=FALSE, sep="\t")


# get cluster frequencies
freq_clust <- table(clust)

# make data frame with labels
labels <- data.frame(cluster = sort(unique(clust)), label = sort(unique(clust)), counts = as.numeric(freq_clust))
labels$proportions <- round(labels$counts/sum(labels$counts) * 100, 2)


write.table(labels, file = file.path(outdir, paste0(prefix, "clustering_labels.xls")), row.names=FALSE, quote=FALSE, sep="\t")


### Code sizes

code_sizes <- table(fsom$mapping[, 1]) 

# Sometimes not all the codes have mapped cells so they will have size 0
code_sizes_full <- rep(0, nrow(data))
names(code_sizes_full) <- 1:nrow(data)

code_sizes_full[names(code_sizes)] <- as.numeric(code_sizes)


### Plot codes as a heatmap

rownames(data) <- 1:nrow(data)

color <- rev(colorRampPalette(brewer.pal(n = 8, name = "Spectral"))(100))
labels_col <- clustering_observables[clustering_observables$clustering_observable, "marker"]

cluster_rows <- fsom_mc_tree
labels_row <- paste0(rownames(data), "  ( ", format(code_sizes_full, big.mark=",", scientific=FALSE), " )")

annotation_row <- data.frame(cluster = factor(fsom_mc))
rownames(annotation_row) <- rownames(data)
annotation_colors <- list(cluster = colors_clusters)  

pheatmap(data[, clustering_observables[clustering_observables$clustering_observable, "mass"]], color = color, cellwidth = 24, cellheight = 12, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, display_numbers = TRUE, number_color = "black", fontsize_number = 6, fontsize_row = 7, fontsize_col = 14, fontsize = 12, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "codes_pheatmap.pdf")))




### Save codes


codes <- data.frame(code_id = rownames(data), cluster = fsom_mc, size = code_sizes_full, data)

write.table(codes, file = file.path(outdir, paste0(prefix, "codes.xls")), row.names=FALSE, quote=FALSE, sep="\t")




























sessionInfo()














################################
### 02_flowsom done!
################################