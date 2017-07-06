

Sys.time()

# Load packages
library(FlowSOM)
library(ConsensusClusterPlus)
library(gdata)
library(pheatmap)
library(RColorBrewer)
library(Rtsne)
library(igraph)
library(ggplot2)
library(ggraph)

##############################################################################
# Test arguments
##############################################################################


prefix='23CD4_01CD4_pca1_cl5_merging5_glmer_binomial_interglht_'
outdir='../carsten_cytof/PD1_project/CK_2016-06-23_01_CD4_mergingNEW2/030_codes'
path_pvs='../carsten_cytof/PD1_project/CK_2016-06-23_01_CD4_mergingNEW2/050_frequencies_codes/23CD4_01CD4_pca1_cl5_frequencies_pvs_glmer_binomial_interglht_top10.xls'
path_fsom='../carsten_cytof/PD1_project/CK_2016-06-23_01_CD4_mergingNEW2/030_heatmaps/23CD4_01CD4_pca1_cl5_fsom.rds'
path_fccp='../carsten_cytof/PD1_project/CK_2016-06-23_01_CD4_mergingNEW2/030_heatmaps/23CD4_01CD4_pca1_cl5_fccp.rds'
path_cluster_merging='../carsten_cytof/PD1_project/CK_2016-06-23_01_CD4_mergingNEW2/010_helpfiles/23CD4_01CD4_pca1_cl5_cluster_merging5.xlsx'
FDR_cutoff='10'

# path_cluster_merging=NULL
args <- NULL

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


if(!file.exists(outdir)) 
  dir.create(outdir, recursive = TRUE)

suffix <- paste0("_top", FDR_cutoff)
FDR_cutoff <- as.numeric(paste0("0.", FDR_cutoff))
FDR_cutoff


# ------------------------------------------------------------
# Load pvalues
# ------------------------------------------------------------

pvs <- read.table(path_pvs, header = TRUE, sep = "\t", as.is = TRUE)

pvs <- pvs[order(pvs$cluster), ]

comparisons <- colnames(pvs)[grep("adjp_", colnames(pvs))]
comparisons


# ------------------------------------------------------------
# Load FlowSOM objects
# ------------------------------------------------------------

fsom <- readRDS(path_fsom)

fccp <- readRDS(path_fccp)

k <- nmetaclusts <- length(fccp)
k

codes <- fsom$codes
ncodes <- nrow(codes)
rownames(codes) <- 1:ncodes

fsom_mc <- fccp[[k]]$consensusClass

if(!length(fsom_mc) == ncodes)
  stop("Some of the codes have zero cells assigned!")

fsom_mc_tree <- fccp[[k]]$consensusTree

code_sizes <- as.numeric(table(fsom$mapping[, 1]))


# ------------------------------------------------------------
# Prepare a color annotation for heatmaps 
# ------------------------------------------------------------

# ggplot palette
gg_color_hue <- function(n){
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 60 , c=100)[1:n]
}

# color blind palette

colors_muted <- c("#DC050C", "#E8601C", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6", "#F1932D", "#F6C141", "#F7EE55", "#4EB265", "#90C987", "#CAEDAB")
color_ramp <- c(colors_muted, gg_color_hue(max(1, k - length(colors_muted))))

colors_clusters <- color_ramp[1:k]
names(colors_clusters) <- 1:k
colors_clusters


# ------------------------------ 
# Annotation for merging or for the original clusters 
# ------------------------------ 


annotation_row <- data.frame(cluster = factor(fsom_mc))
rownames(annotation_row) <- 1:ncodes

annotation_colors <- list(cluster = colors_clusters)
rows_order <- order(fsom_mc)


if(!is.null(path_cluster_merging)){
  
  ### Read in cluster merging file
  cm <- gdata::read.xls(path_cluster_merging)
  
  if(!all(c("old_cluster", "label", "new_cluster") %in% colnames(cm)))
    stop("Merging file must contain 'old_cluster', 'label' and 'new_cluster' columns!")
  
  ### Remove spaces in labels bcs they are problematic...
  cm$label <- factor(cm$label, labels = gsub(" ", "_", levels(cm$label))) 
  
  cm_unique <- unique(cm[, c("label", "new_cluster")])
  cm_unique <- cm_unique[order(cm_unique$new_cluster), ]
  
  ### Add merging to the annotation
  mm <- match(annotation_row$cluster, cm$old_cluster)
  annotation_row$cluster_merging <- cm$label[mm]
  annotation_row$cluster_merging <- factor(annotation_row$cluster_merging, levels = cm_unique$label)
  
  ### Add colors for merging
  color_ramp <- c(colors_muted, gg_color_hue(max(1, nlevels(cm_unique$label) - length(colors_muted))))
  
  colors_clusters_merging <- color_ramp[1:nlevels(cm_unique$label)]
  names(colors_clusters_merging) <- cm_unique$label
  
  annotation_colors[["cluster_merging"]] <- colors_clusters_merging
  
  rows_order <- order(annotation_row$cluster_merging, annotation_row$cluster)
  
  ### Drop the "drop" cluster
  rows_order <- rows_order[annotation_row$cluster_merging[rows_order] != "drop"]
  
}




# ------------------------------------------------------------
# Row clustering from the fccp object
# ------------------------------------------------------------

### This clustering is used in all the heatmaps
cluster_rows <- fsom_mc_tree

# ------------------------------------------------------------
# Heatmaps with adjusted p-values for the codes
# ------------------------------------------------------------


### Plot all the p-values

pvs_heat <- as.matrix(pvs[, comparisons, drop = FALSE])
head(pvs_heat)
rownames(pvs_heat) <- pvs$cluster


labels_row <- rownames(pvs_heat)
labels_col <- colnames(pvs_heat)


## With row clustering
pheatmap(pvs_heat[cluster_rows$order, , drop = FALSE], cellwidth = 40, cellheight = 8, border_color = NA, color = c("grey50", "grey70", "grey90"), breaks = c(0, 0.05, 0.1, 1), legend_breaks = c(0, 0.05, 0.1, 1), legend = FALSE, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row[cluster_rows$order], display_numbers = TRUE, number_format = "%.02e", number_color = "black", fontsize_row = 8, fontsize_col = 14, fontsize = 6, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "pheatmap_codes_pvs_row_clust_all.pdf")))


## No row clustering
pheatmap(pvs_heat[rows_order, , drop = FALSE], cellwidth = 40, cellheight = 8, border_color = NA, color = c("grey50", "grey70", "grey90"), breaks = c(0, 0.05, 0.1, 1), legend_breaks = c(0, 0.05, 0.1, 1), legend = FALSE, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row[rows_order], display_numbers = TRUE, number_format = "%.02e", number_color = "black", fontsize_row = 8, fontsize_col = 14, fontsize = 6, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "pheatmap_codes_pvs_no_clust_all.pdf")))



### Plot p-values per comparison

# for(i in 1:length(comparisons)){
#   # i = 1
#   
#   comparison <- comparisons[i]
#   print(comparison)
#   comparison_suffix <- paste0(gsub("adjp_", "", comparison))
#   
#   pvs_heat <- as.matrix(pvs[, comparison, drop = FALSE])
#   head(pvs_heat)
#   rownames(pvs_heat) <- pvs$cluster
#   
#   labels_row <- rownames(pvs_heat)
#   labels_col <- colnames(pvs_heat)
#   
#   ## With row clustering
#   pheatmap(pvs_heat[cluster_rows$order, , drop = FALSE], cellwidth = 40, cellheight = 8, border_color = NA, color = c("grey50", "grey70", "grey90"), breaks = c(0, 0.05, 0.1, 1), legend_breaks = c(0, 0.05, 0.1, 1), legend = FALSE, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row[cluster_rows$order], display_numbers = FALSE, number_format = "%.02e", number_color = "black", fontsize_row = 8, fontsize_col = 14, fontsize = 6, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "pheatmap_codes_pvs_row_clust_",comparison_suffix ,".pdf")))
#   
#   
#   ## No row clustering
#   pheatmap(pvs_heat[rows_order, , drop = FALSE], cellwidth = 40, cellheight = 8, border_color = NA, color = c("grey50", "grey70", "grey90"), breaks = c(0, 0.05, 0.1, 1), legend_breaks = c(0, 0.05, 0.1, 1), legend = FALSE, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row[rows_order], display_numbers = FALSE, number_format = "%.02e", number_color = "black", fontsize_row = 8, fontsize_col = 14, fontsize = 6, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "pheatmap_codes_pvs_no_clust_", comparison_suffix, ".pdf")))
#   
#   
# }



# ------------------------------------------------------------
# t-SNE of codes with adjusted p-values
# ------------------------------------------------------------

rand_seed <- 1234

### Run tSNE
set.seed(rand_seed)
rtsne_out <- Rtsne(codes, perplexity = 20, pca = FALSE, max_iter = 1000, verbose = TRUE)

pvs_discrete <- lapply(1:length(comparisons), function(i){
  cut(pvs[, comparisons[i]], c(0, 0.05, 0.1, 1))
})

pvs_discrete <- data.frame(pvs_discrete)
colnames(pvs_discrete) <- comparisons

colors_pvs <- c("grey20", "grey60", "grey100")
names(colors_pvs) <- levels(pvs_discrete[, 1])


ggdf <- data.frame(codes, dim1 = rtsne_out$Y[, 1], dim2 = rtsne_out$Y[, 2], pvs_discrete, annotation_row, code_sizes = code_sizes)


if(is.null(path_cluster_merging)){
  clustering <- "cluster"
  colors <- colors_clusters
}else{
  colors <- colors_clusters_merging
  clustering <- "cluster_merging"
  
  ### Drop the "drop" cluster
  codes2drop <- ggdf$cluster_merging == "drop"
  ggdf <- ggdf[!codes2drop, ]
  
}


for(i in 1:length(comparisons)){
  # i = 1
  
  comparison <- comparisons[i]
  print(comparison)
  comparison_suffix <- paste0(gsub("adjp_", "", comparison))
  
  ggp <- ggplot(ggdf, aes_string(x = "dim1", y = "dim2", color = clustering, size = "code_sizes", fill = comparison)) +
    geom_point(alpha = 0.9, shape = 21, stroke = 2) +
    labs(x = "t-SNE1", y = "t-SNE2") + 
    theme_bw() +
    theme(axis.text = element_text(size = 12), 
      axis.title  = element_text(size = 15),
      strip.text = element_text(size = 15, hjust = 0),
      strip.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_rect(size = 0.3, linetype = "solid", color = "black")) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors_pvs) +
    guides(colour = guide_legend(override.aes = list(size = 5)))
  
  pdf(file.path(outdir, paste0(prefix, "tsne_codes_pvs_", comparison_suffix, ".pdf")), width = 9, height = 7)       
  print(ggp)
  dev.off()
  
  
}








# ------------------------------------------------------------
# MST of codes with adjusted p-values
# ------------------------------------------------------------


adjacency <- stats::dist(ggdf[, colnames(codes)], method = "euclidean")

fullGraph <- igraph::graph.adjacency(as.matrix(adjacency), mode = "undirected", weighted = TRUE)

MST_graph <- igraph::minimum.spanning.tree(fullGraph)

MST_l <- igraph::layout.kamada.kawai(MST_graph) 


# layout <- MST_l
# lty <- 1
# 
# vertex_sizes <- ggdf$code_sizes/max(ggdf$code_sizes) * 10
# # vertex_sizes[vertex_sizes < 2] <- 2
# 
# 
# if(is.null(path_cluster_merging)){
#   vertex_colors <- colors_clusters[ggdf$cluster]
# }else{
#   vertex_colors <- colors_clusters_merging[ggdf$cluster_merging]
# }
# 
# mark.col <- adjustcolor(vertex_colors, alpha = 0.5)
# names(mark.col) <- paste0(names(vertex_colors), 1:length(vertex_colors))
# 
# mark.groups <- lapply(1:length(vertex_colors), function(x){
#   x
# })
# 
# names(mark.groups) <- names(mark.col)
# 
# 
# 
# ### Plot MST colored by p-values
# 
# for(i in 1:length(comparisons)){
#   # i = 1
#   
#   comparison <- comparisons[i]
#   print(comparison)
#   comparison_suffix <- paste0(gsub("adjp_", "", comparison))
#   
#   vertex_fill <- colors_pvs[ggdf[, comparison]]
#   
#   pdf(file.path(outdir, paste0(prefix, "mst_codes_pvs_", comparison_suffix, ".pdf")), width = 7, height = 7)
#   
#   igraph::plot.igraph(MST_graph, layout = layout, vertex.size = vertex_sizes, vertex.label = NA, vertex.label.cex = 0.5, vertex.color = vertex_fill, edge.lty = lty, mark.groups = mark.groups, mark.border = NA, mark.col = mark.col, mark.shape = 1)
#   
#   dev.off()
#   
# }


### Plot MST colored by cell population


# pdf(file.path(outdir, paste0(prefix, "mst_codes_pvs_clusters.pdf")), width = 7, height = 7)
# 
# igraph::plot.igraph(MST_graph, layout = layout, vertex.size = vertex_sizes, vertex.label = NA, vertex.label.cex = 0.5, vertex.color = vertex_colors, edge.lty = lty)
# 
# dev.off()




### Plotting with ggraph


V(MST_graph)$code_sizes <- ggdf$code_sizes
V(MST_graph)$clustering <- as.character(ggdf[, clustering])



for(i in 1:length(comparisons)){
  # i = 1
  
  comparison <- comparisons[i]
  print(comparison)
  comparison_suffix <- paste0(gsub("adjp_", "", comparison))
  
  V(MST_graph)$comparison <- as.character(ggdf[, comparison])

  ggp <- ggraph(MST_graph, layout = 'igraph', algorithm = 'kk') + 
    geom_edge_link() + 
    geom_node_point(aes(size = code_sizes, fill = comparison, color = clustering), shape = 21, stroke = 2) +
    theme_bw() +
    theme(axis.text = element_blank(), 
      axis.title  = element_blank(),
      axis.ticks  = element_blank(),
      strip.text = element_text(size = 15, hjust = 0),
      strip.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank()) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors_pvs) +
    guides(colour = guide_legend(override.aes = list(size = 5)))
  
  
  pdf(file.path(outdir, paste0(prefix, "ggraph_codes_pvs_", comparison_suffix, ".pdf")), width = 9, height = 7)       
  print(ggp)
  dev.off()
  
  
}











sessionInfo()


