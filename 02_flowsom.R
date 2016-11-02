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
  dir.create(outdir)

if(all(!grepl("linkage=", args))){
  linkage <- "average"
}


# ------------------------------------------------------------
# Colors for MST maps
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
fsom <- FlowSOM::SOM(ef)



# -------------------------------------
# metaClustering_consensus 
# -------------------------------------


# consensus clustering that is reproducible with seed

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


# Save clustering results


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



### Plot codes as a heatmap

rownames(data) <- 1:nrow(data)
data_order <- match(colnames(data), clustering_observables$mass)

color <- rev(colorRampPalette(brewer.pal(n = 8, name = "Spectral"))(100))
labels_col <- clustering_observables$marker[1:length(data_order)]
annotation_colors <- list(cluster = colors_clusters)  

cluster_rows <- fsom_mc_tree
labels_row <- 1:nrow(data)
annotation_row <- data.frame(cluster = factor(fsom_mc))
rownames(annotation_row) <- 1:nrow(data)


pheatmap(data[, data_order], color = color, cellwidth = 24, cellheight = 12, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, display_numbers = TRUE, number_color = "black", fontsize_number = 6, fontsize_row = 7, fontsize_col = 14, fontsize = 12, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "codes_pheatmap.pdf")))




# -------------------------------------
# BuildMST
# -------------------------------------


adjacency <- stats::dist(data, method = "euclidean")

fullGraph <- igraph::graph.adjacency(as.matrix(adjacency), mode = "undirected", weighted = TRUE)

MST_graph <- igraph::minimum.spanning.tree(fullGraph)

MST_l <- igraph::layout.kamada.kawai(MST_graph) 


layout <- MST_l
lty <- 1
code_sizes <- table(fsom$mapping[,1]) 
vertex_sizes <- code_sizes/max(code_sizes) * 10
vertex_sizes[vertex_sizes < 2] <- 2

### Plot MST

vertex_colors <- colors_clusters[as.character(fsom_mc)]

pdf(file.path(outdir, paste0(prefix, "codes_mst.pdf")), width = 7, height = 7)

igraph::plot.igraph(MST_graph, layout = layout, vertex.size = vertex_sizes, vertex.label = NA, vertex.label.cex = 0.5, vertex.color = vertex_colors, edge.lty = lty)

dev.off()


# -------------------------------------
# Dimention reduction 
# -------------------------------------

dr <- list()

set.seed(rand_seed)
dr[["isomap"]] <- cytofkit::cytof_dimReduction(data, method="isomap")
set.seed(rand_seed)
dr[["tsne"]] <- cytofkit::cytof_dimReduction(data, method="tsne", perplexity = 20, tsneSeed = rand_seed)
set.seed(rand_seed)
dr[["diffusion"]] <- cytofkit::cytof_dimReduction(data, method="diffusionmap")
set.seed(rand_seed)
dr[["pca"]] <- cytofkit::cytof_dimReduction(data, method="pca")

dr[["graph"]] <- MST_l

dr_methods <- names(dr)

for(i in 1:length(dr)){
  # i = 1
  
  dr_data <- data.frame(dr[[i]])
  colnames(dr_data) <- c("dim1", "dim2")
  
  dr_data$cluster <- factor(fsom_mc)
  dr_data$size <- as.numeric(vertex_sizes/10)
  
  ggdf <- dr_data
  
  ggp <- ggplot(ggdf,  aes(x = dim1, y = dim2, color = cluster, size = size)) +
    geom_point(alpha = 0.7) +
    labs(x = "dim1", y="dim2") + 
    theme_bw() +
    theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold"), legend.key = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_color_manual(values = colors_clusters[as.character(fsom_mc)]) +
    guides(colour = guide_legend(override.aes = list(size = 5)))
  
  pdf(file.path(outdir, paste0(prefix, "codes_", dr_methods[[i]], ".pdf")), width = 9, height = 7)                 
  print(ggp)
  dev.off()
  
}




# --------------------------------------------------------------------------

# --------------------------------------------------------------------------


# -------------------------------------
# hierarchical clustering 
# -------------------------------------


# fsom_hc <- hclust(dist(data), method = linkage)
# 
# fsom_ct <- cutree(fsom_hc, k)
# 
# clust_hc <- fsom_ct[fsom$mapping[,1]]
# 
# 
# # Save clustering results
# 
# clust_out <- data.frame(cluster = clust_hc, cell_id = cell_id, sample_id = samp, stringsAsFactors = FALSE)
# write.table(clust_out, file = file.path(outdir, paste0(prefix, "hc_clustering.xls")), row.names=FALSE, quote=FALSE, sep="\t")
# 
# 
# # get cluster frequencies
# freq_clust <- table(clust_hc)
# 
# # make data frame with labels
# labels <- data.frame(cluster = sort(unique(clust_hc)), label = sort(unique(clust_hc)), counts = as.numeric(freq_clust))
# labels$proportions <- round(labels$counts/sum(labels$counts) * 100, 2)
# 
# 
# write.table(labels, file = file.path(outdir, paste0(prefix, "hc_clustering_labels.xls")), row.names=FALSE, quote=FALSE, sep="\t")
# 
# 
# 
# ### Plot MST
# 
# vertex_colors <- colors_clusters[as.character(fsom_ct)]
# 
# pdf(file.path(outdir, paste0(prefix, "hc_codes_mst.pdf")), width = 7, height = 7)
# 
# igraph::plot.igraph(MST_graph, layout = layout, vertex.size = vertex_sizes, vertex.label = NA, vertex.label.cex = 0.5, vertex.color = vertex_colors, edge.lty = lty)
# 
# dev.off()
# 
# 
# 
# ### Plot codes
# rownames(data) <- 1:nrow(data)
# data_order <- match(colnames(data), clustering_observables$mass)
# 
# color <- rev(colorRampPalette(brewer.pal(n = 8, name = "Spectral"))(100))
# labels_col <- clustering_observables$marker[1:length(data_order)]
# labels_row <- 1:nrow(data)
# annotation_colors <- list(cluster = colors_clusters)  
# 
# cluster_rows <- fsom_hc
# annotation_row <- data.frame(cluster = factor(fsom_ct))
# rownames(annotation_row) <- 1:nrow(data)
# 
# 
# pheatmap(data[, data_order], color = color, cellwidth = 24, cellheight = 12, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, display_numbers = TRUE, number_color = "black", fontsize_number = 6, fontsize_row = 7, fontsize_col = 14, fontsize = 12, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "hc_codes_pheatmap.pdf")))































sessionInfo()














################################
### 02_flowsom done!
################################