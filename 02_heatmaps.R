##############################################################################
## <<02_heatmaps.R>>

# BioC 3.3
# Created 27 July 2016
# Updated 24 Aug 2016

##############################################################################
Sys.time()
##############################################################################

# Load packages
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
library(ConsensusClusterPlus)

##############################################################################
# Test arguments
##############################################################################

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
# heatmap_prefix='23_01_pca1_cl20_raw_'
# heatmap_outdir='030_heatmaps'
# path_data='010_data/23_01_expr_raw.rds'
# path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.xlsx'
# path_clustering_observables='030_heatmaps/23_01_pca1_clustering_observables.xls'
# path_clustering='030_heatmaps/23_01_pca1_cl20_clustering.xls'
# path_clustering_labels='030_heatmaps/23_01_pca1_cl20_clustering_labels.xls'
# path_marker_selection='030_heatmaps/23_01_pca1_cl20_marker_selection.txt'
# aggregate_fun='median'
# pheatmap_palette='YlGnBu' # RdYlBu
# pheatmap_palette_rev=FALSE
# pheatmap_scale=TRUE


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

prefix <- heatmap_prefix
outdir <- heatmap_outdir

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
# Load clustering data
# ------------------------------------------------------------

# clustering
clustering <- read.table(path_clustering, header = TRUE, sep = "\t", as.is = TRUE)

clust <- clustering[, "cluster"]
names(clust) <- clustering[, "cell_id"]

# clustering labels
labels <- read.table(path_clustering_labels, header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))
rownames(labels) <- labels$cluster

# clustering observables
clustering_observables <- read.table(path_clustering_observables, header = TRUE, sep = "\t", as.is = TRUE)
rownames(clustering_observables) <- clustering_observables$mass

clust_observ <- clustering_observables[clustering_observables$clustering_observable, "mass"]

# ------------------------------------------------------------
# load marker selection for plotting on the heatmaps
marker_selection <- NULL

if(file.exists(file.path(path_marker_selection))){
  
  marker_selection <- read.table(file.path(path_marker_selection), header = TRUE, sep = "\t", as.is = TRUE)
  marker_selection <- marker_selection[, 1]
  
  if(!all(marker_selection %in% clustering_observables$marker))
    stop("Marker selection is wrong")
  
}


# ------------------------------------------------------------
# get the isotope and antigen for fcs markers

m <- match(fcs_colnames, clustering_observables$mass)

fcs_panel <- data.frame(fcs_colname = fcs_colnames, Isotope = clustering_observables$mass[m], Antigen = clustering_observables$marker[m], stringsAsFactors = FALSE)


# ------------------------------------------------------------

# Indeces of observables used for clustering 
scols <- which(fcs_colnames %in% clust_observ)

# Indeces of other observables
xcols <- which(!fcs_colnames %in% clust_observ)


# ordered by decreasing pca score
if("avg_score" %in% colnames(clustering_observables)){
  scols <- scols[order(clustering_observables[fcs_colnames[scols], "avg_score"], decreasing = TRUE)]
  xcols <- xcols[order(clustering_observables[fcs_colnames[xcols], "avg_score"], decreasing = TRUE)]
}


# ------------------------------------------------------------
# Keep expression and clustering results for the cells that are common in both
# ------------------------------------------------------------

common_cells <- intersect(clustering[, "cell_id"], expr[, "cell_id"])

samp <- expr[expr[, "cell_id"] %in% common_cells, "sample_id"]
e <- expr[expr[, "cell_id"] %in% common_cells, fcs_colnames]
clust <- clustering[clustering[, "cell_id"] %in% common_cells, "cluster"]


# ------------------------------------------------------------
# Get the median expression
# ------------------------------------------------------------

colnames(e) <- fcs_panel$Antigen

a <- aggregate(e, by = list(clust), FUN = aggregate_fun)

mlab <- match(a$Group.1, labels$cluster)
rownames(a) <- labels$label[mlab]

colnames(a)[1] <- "cluster"


# get cluster frequencies
freq_clust <- table(clust)


### Save cluster frequencies and the median expression

clusters_out <- data.frame(cluster = names(freq_clust), label = labels[names(freq_clust), "label"], counts = as.numeric(freq_clust), frequencies = as.numeric(freq_clust)/sum(freq_clust), a[, fcs_panel$Antigen[c(scols, xcols)]])

write.table(clusters_out, file.path(outdir, paste0(prefix, "clusters.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)




# ------------------------------------------------------------
# pheatmaps of median expression
# ------------------------------------------------------------

### Cluster clustering is based on all markers here
expr <- as.matrix(a[, fcs_panel$Antigen[c(scols, xcols)]])
rownames(expr) <- labels$label

cluster_rows <- hclust(dist(expr), method = "average")

labels_row <- paste0(as.character(labels$label), " (", round(as.numeric(freq_clust)/sum(freq_clust)*100, 2), "%)")
labels_col <- colnames(expr)

if(pheatmap_palette_rev){
  color <- colorRampPalette(rev(brewer.pal(n = 8, name = pheatmap_palette)))(100)
}else{
  color <- colorRampPalette(brewer.pal(n = 8, name = pheatmap_palette))(100)
}


pheatmap(expr, color = color, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, display_numbers = TRUE, number_color = "black", fontsize_number = 7, gaps_col = length(scols), fontsize_row = 10, fontsize_col = 10, fontsize = 7, filename = file.path(outdir, paste0(prefix, "pheatmap_all_row_clust.pdf")), width = 10, height = 7)

pheatmap(expr, color = color, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, display_numbers = TRUE, number_color = "black", fontsize_number = 7, gaps_col = length(scols), fontsize_row = 10, fontsize_col = 10, fontsize = 7, filename = file.path(outdir, paste0(prefix, "pheatmap_all.pdf")), width = 10, height = 7)


## Plot only the selected markers
if(!is.null(marker_selection)){
  expr_sub <- expr[, marker_selection, drop = FALSE]
  labels_col_sub <- colnames(expr_sub)
  
  pheatmap(expr_sub, color = color, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col_sub, labels_row = labels_row, fontsize_number = 7, fontsize_row = 10, fontsize_col = 10, fontsize = 7, filename = file.path(outdir, paste0(prefix, "pheatmap_s1_row_clust.pdf")), width = 10, height = 7)
  
  pheatmap(expr_sub, color = color, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col_sub, labels_row = labels_row, fontsize_number = 7, fontsize_row = 10, fontsize_col = 10, fontsize = 7, filename = file.path(outdir, paste0(prefix, "pheatmap_s1.pdf")), width = 10, height = 7)
  
}



# ------------------------------------------------------------
# pheatmaps of median expression scalled by marker (column)
# ------------------------------------------------------------

if(pheatmap_scale){
  
  ### Cluster clustering is based on all markers here
  expr <- as.matrix(a[, fcs_panel$Antigen[c(scols, xcols)]])
  rownames(expr) <- labels$label
  
  
  expr_scaled <- apply(expr, 2, function(x){
    (x-min(x))/(max(x)-min(x))
  })
  
  cluster_rows <- hclust(dist(expr_scaled), method = "average")
  
  labels_row <- paste0(as.character(labels$label), " (", round(as.numeric(freq_clust)/sum(freq_clust)*100, 2), "%)")
  labels_col <- colnames(expr_scaled)
  
  color <- colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlGn")))(100)
  
  pheatmap(expr_scaled, color = color, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, display_numbers = TRUE, number_color = "black", fontsize_number = 7, gaps_col = length(scols), fontsize_row = 10, fontsize_col = 10, fontsize = 7, filename = file.path(outdir, paste0(prefix, "pheatmap_all_row_clust_scale.pdf")), width = 10, height = 7)
  
  pheatmap(expr_scaled, color = color, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, display_numbers = TRUE, number_color = "black", fontsize_number = 7, gaps_col = length(scols), fontsize_row = 10, fontsize_col = 10, fontsize = 7, filename = file.path(outdir, paste0(prefix, "pheatmap_all_scale.pdf")), width = 10, height = 7)
  
  
  ## Plot only the selected markers
  if(!is.null(marker_selection)){
    expr_sub <- expr_scaled[, marker_selection, drop = FALSE]
    labels_col_sub <- colnames(expr_sub)
    
    pheatmap(expr_sub, color = color, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col_sub, labels_row = labels_row, fontsize_number = 7, fontsize_row = 10, fontsize_col = 10, fontsize = 7, filename = file.path(outdir, paste0(prefix, "pheatmap_s1_row_clust_scale.pdf")), width = 10, height = 7)
    
    pheatmap(expr_sub, color = color, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col_sub, labels_row = labels_row, fontsize_number = 7, fontsize_row = 10, fontsize_col = 10, fontsize = 7, filename = file.path(outdir, paste0(prefix, "pheatmap_s1_scale.pdf")), width = 10, height = 7)
    
  }
  
  
}



# ------------------------------------------------------------
# Plot expression of markers in each cluster - distros plots
# ------------------------------------------------------------


plotting_wrapper <- function(e, suffix){
  
  df <- data.frame(e, clust = factor(clust))
  dfm <- melt(df, id.vars = c("clust"))
  
  dfm$clust <- factor(dfm$clust, labels = labels$label)
  
  # ggp <- ggplot(dfm, aes(x=value)) + 
  #   geom_density(adjust=3, fill = "black", alpha = 0.3) + 
  #   facet_grid(clust ~ variable, scales = "free") +
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # 
  # pdf(file.path(outdir, paste0(prefix, "distros_raw_in.pdf")), w = ncol(e), h = nrow(labels))
  # print(ggp)
  # dev.off()
  
  ## with free scales using facet_wrap 
  dfm$variable_clust <- interaction(dfm$variable, dfm$clust, lex.order = FALSE)
  
  ggp <- ggplot(dfm, aes(x=value)) +
    geom_density(adjust=1, fill = "black", alpha = 0.3) +
    facet_wrap(~ variable_clust, nrow = nlevels(dfm$clust), scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text = element_text(size = 6), strip.text = element_text(size = 9))
  
  pdf(file.path(outdir, paste0(prefix, "distrosfree", suffix, ".pdf")), w = ncol(e)*3/2, h = nrow(labels)*3/2)
  print(ggp)
  dev.off()
  
  return(NULL)
  
}



# ## Expression, included observables
# 
# plotting_wrapper(e = e[, fcs_panel$Antigen[scols]], suffix = "_in")
# 
# 
# if(length(xcols) > 0){
# 
#   # Expression, excluded observables
# 
#   plotting_wrapper(e = e[, fcs_panel$Antigen[xcols]], suffix = "_ex")
# 
# }


# ------------------------------------------------------------
# Plot expression of markers for all data and strat. per sample
# ------------------------------------------------------------


plotting_wrapper2 <- function(e, suffix){
  
  df <- data.frame(samp = samp, e)
  dfm <- melt(df, id.var = "samp")
  
  # add group info
  mm <- match(dfm$samp, md$shortname)
  dfm$group <- factor(md$condition[mm])
  
  ggp <- ggplot(dfm, aes(x=value)) + 
    geom_density(adjust = 1, fill = "black", alpha = 0.3) + 
    facet_wrap(~ variable, nrow = 3, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  pdf(file.path(outdir, paste0(prefix, "distrosmer", suffix,".pdf")), w = ncol(e), h = 10)
  print(ggp)
  dev.off()
  
  ## create colors per group not per sample
  colors <- unique(md[, c("condition", "color")])
  
  mm <- match(levels(dfm$samp), md$shortname)
  groups <- factor(md$condition[mm])
  color_values <- colorRampPalette(brewer.pal(12,"Paired"))(12)[c(1,3,5,2,4,6)]
  color_values <- color_values[as.numeric(groups)]
  names(color_values) <- levels(dfm$samp)
  
  ggp <- ggplot(dfm, aes(x=value, color = samp)) + 
    geom_density(adjust = 1) + 
    facet_wrap(~ variable, nrow = 3, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title = element_blank(), legend.position = "bottom") +
    guides(color = guide_legend(nrow = 2)) +
    scale_color_manual(values = color_values)
  
  pdf(file.path(outdir, paste0(prefix, "distrosgrp", suffix,".pdf")), w = ncol(e), h = 10)
  print(ggp)
  dev.off()
  
  
  return(NULL)
  
}




# ## Expression, included observables
# 
# plotting_wrapper2(e = e[, fcs_panel$Antigen[scols]], suffix = "_in")
# 
# if(length(xcols) > 0){
# 
#   # Expression, excluded observables
# 
#   plotting_wrapper2(e = e[, fcs_panel$Antigen[xcols]], suffix = "_ex")
# 
# }











sessionInfo()













################################
### 02_heatmaps done!
################################