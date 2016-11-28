##############################################################################
## <<06_cytokines_bimatrix_clustering.R>>

# BioC 3.3
# Created 27 July 2016
# Updated 3 Nov 2016

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
library(pheatmap)
library(RColorBrewer)

##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_02_CD4_merging2'
clust_prefix='23CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_cl40_'
clust_outdir='060_cytokines_bimatrix_flowsom/01_clustering'
path_data='060_cytokines_bimatrix/01_clustering/23CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_bimatrix.txt'
path_clustering_observables='060_cytokines_bimatrix/01_clustering/23CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_clustering_observables.xls'
som_dim=10
nmetaclusts=40

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-29_02_CD4_merging'
clust_prefix='29CD4_02CD4_pca1_merging3_EM_CM_cytCM_raw2_cl25_'
clust_outdir='060_cytokines_bimatrix/01_clustering'
path_data='060_cytokines_bimatrix/01_clustering/29CD4_02CD4_pca1_merging3_EM_CM_cytCM_raw2_bimatrix.txt'
path_clustering_observables='060_cytokines_bimatrix/01_clustering/29CD4_02CD4_pca1_merging3_EM_CM_cytCM_raw2_clustering_observables.xls'
som_dim=5
nmetaclusts=25

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

prefix <- clust_prefix
outdir <- clust_outdir

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
fsom <- FlowSOM::SOM(ef, xdim = som_dim, ydim = som_dim)

data <- fsom$codes
rownames(data) <- 1:nrow(data)



if(som_dim^2 > nmetaclusts){
  
  # -------------------------------------
  # metaClustering_consensus 
  # -------------------------------------
  # consensus clustering that is reproducible with seed
  
  ### Sometimes not all the codes are used in mapping
  k <- nmetaclusts
  
  pdf(file.path(outdir, paste0(prefix, "ConsensusClusterPlus.pdf")), width = 7, height = 7)
  
  results <- ConsensusClusterPlus::ConsensusClusterPlus(t(data),
    maxK = k, reps = 100, pItem = 0.9, pFeature = 1, title = tempdir(),
    plot = NULL, verbose = FALSE, clusterAlg = "hc", innerLinkage = linkage, finalLinkage = "average", distance = "euclidean", seed = rand_seed)
  
  dev.off()
  
  # get cluster ids
  fsom_mc <- results[[k]]$consensusClass
  fsom_mc_tree <- results[[k]]$consensusTree
  
  clust <- fsom_mc[fsom$mapping[,1]]
  
  # elements for the heatmap
  annotation_row <- data.frame(cluster = factor(fsom_mc))
  rownames(annotation_row) <- rownames(data)
  annotation_colors <- list(cluster = colors_clusters)
  
  ### Save clustering results
  save(fsom, file = file.path(outdir, paste0(prefix, "fsom.rda")))
  save(fsom_mc, file = file.path(outdir, paste0(prefix, "fsom_mc.rda")))
  
}else{
  
  clust <- fsom$mapping[,1]
  
  fsom_mc <- 1:nrow(data)
  fsom_mc_tree <- hclust(dist(data), method = linkage)
  
  # elements for the heatmap
  annotation_row <- NA
  annotation_colors = NA
  
}



### Save clustering results

clust_out <- data.frame(cluster = clust, cell_id = cell_id, sample_id = samp, stringsAsFactors = FALSE)
write.table(clust_out, file = file.path(outdir, paste0(prefix, "clustering.xls")), row.names=FALSE, quote=FALSE, sep="\t")


### Code sizes

code_sizes <- table(fsom$mapping[, 1]) 

# Sometimes not all the codes have mapped cells so they will have size 0
code_sizes_full <- rep(0, nrow(data))
names(code_sizes_full) <- 1:nrow(data)

code_sizes_full[names(code_sizes)] <- as.numeric(code_sizes)



### Plot codes as a heatmap

color <- rev(colorRampPalette(brewer.pal(n = 8, name = "Spectral"))(100))
labels_col <- clustering_observables[clustering_observables$clustering_observable, "marker"]

cluster_rows <- fsom_mc_tree

labels_row <- paste0(rownames(data), "  ( ", format(code_sizes_full, big.mark=",", scientific=FALSE), " )")


pheatmap(data[, clustering_observables[clustering_observables$clustering_observable, "mass"]], color = color, cellwidth = 24, cellheight = 12, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, display_numbers = TRUE, number_color = "black", fontsize_number = 6, fontsize_row = 7, fontsize_col = 14, fontsize = 12, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "codes_pheatmap.pdf")))


### Save codes

codes <- data.frame(code_id = rownames(data), cluster = fsom_mc, size = code_sizes_full, data)

write.table(codes, file = file.path(outdir, paste0(prefix, "codes.xls")), row.names=FALSE, quote=FALSE, sep="\t")



### Save data frames with labels

# get cluster frequencies
freq_clust <- table(clust)
prop_clust <- round(freq_clust/sum(freq_clust)*100, 2)

# get the names of cytokines that are positive in a given cluster

a <- aggregate(ef, by = list(clust), FUN = "mean")

data_new_labels <- a[, clustering_observables[clustering_observables$clustering_observable, "mass"]]
colnames(data_new_labels) <- clustering_observables[clustering_observables$clustering_observable, "marker"]

cluster_binames <- apply(data_new_labels, 1, function(x){
  
  paste0(rev(colnames(data_new_labels)[x > 0.75]), collapse = "/")
  
})


cluster <- sort(unique(clust))

label <- paste0(cluster, " - ", cluster_binames, " (", prop_clust, ")")

labels <- data.frame(cluster = cluster, label = label, stringsAsFactors = FALSE)


labels <- data.frame(cluster = cluster, label = label, counts = as.numeric(freq_clust))
labels$proportions <- round(labels$counts/sum(labels$counts) * 100, 2)


write.table(labels, file = file.path(outdir, paste0(prefix, "clustering_labels.xls")), row.names=FALSE, quote=FALSE, sep="\t")






### Plot heatmap of generated clusters


color <- rev(colorRampPalette(brewer.pal(n = 8, name = "Spectral"))(100))
labels_col <- colnames(data_new_labels)

cluster_rows <- hclust(dist(data_new_labels), method = linkage)

cluster <- sort(unique(clust))

labels_row <- paste0(cluster, " - ", cluster_binames, " (", prop_clust, ")")

annotation_row <- data.frame(cluster = factor(cluster))
rownames(annotation_row) <- rownames(data_new_labels)

annotation_colors <- list(cluster = colors_clusters)


pheatmap(data_new_labels, color = color, cellwidth = 20, cellheight = 20, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, display_numbers = TRUE, number_color = "black", fontsize_number = 6, fontsize_row = 10, fontsize_col = 10, fontsize = 10, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "clusters_pheatmap.pdf")))


# Error in annotation_colors[[colnames(annotation)[i]]] :
# subscript out of bounds















sessionInfo()














################################
### 06_cytokines_bimatrix_clustering done!
################################