##############################################################################
## <<02_heatmaps.R>>

# BioC 3.3
# Created 27 July 2016
# Updated 18 Aug 2016

##############################################################################
Sys.time()
##############################################################################

# Load packages
library(flowCore)
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
# path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.xlsx'
# heatmap_prefix='23_01_pca1_cl20_'
# path_clustering_observables='23_01_pca1_clustering_observables.xls'
# path_clustering='23_01_pca1_cl20_clustering.xls'
# path_clustering_labels='23_01_pca1_cl20_clustering_labels.xls'


# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/FACS_data'
# path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/FACS_data/metadata_facs.xlsx'
# heatmap_prefix='facs_pca1_cl20_'
# path_clustering_observables='facs_pca1_clustering_observables.xls'
# path_clustering='facs_pca1_cl20_clustering.xls'
# path_clustering_labels='facs_pca1_cl20_clustering_labels.xls'


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


# ------------------------------------------------------------
# define directories
# ------------------------------------------------------------

fcsDir <- "010_cleanfcs"; if( !file.exists(fcsDir) ) dir.create(fcsDir)
hmDir <- "030_heatmaps"; if( !file.exists(hmDir) ) dir.create(hmDir)


# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

# read metadata
md <- read.xls(path_metadata, stringsAsFactors=FALSE)

# define FCS file names
f <- file.path(fcsDir, md$filename)
names(f) <- md$shortname

# read raw FCS files in
fcs <- lapply(f, read.FCS)

fcs_colnames <- colnames(fcs[[1]])

# ------------------------------------------------------------
# Load more data
# ------------------------------------------------------------

# clustering
clust <- read.table(file.path(hmDir, path_clustering), header = TRUE, sep = "\t", as.is = TRUE)
clust <- clust[, 1]

# clustering labels
labels <- read.table(file.path(hmDir, path_clustering_labels), header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))
rownames(labels) <- labels$cluster


# clustering_observables
if(!grepl("/", path_clustering_observables)){
  clustering_observables <- read.table(file.path(hmDir, path_clustering_observables), header = TRUE, sep = "\t", as.is = TRUE)
}else{
  clustering_observables <- read.table(path_clustering_observables, header = TRUE, sep = "\t", as.is = TRUE)
}

rownames(clustering_observables) <- clustering_observables$mass

clust_observ <- clustering_observables[clustering_observables$clustering_observable, "mass"]


# -------------------------------------
# get the isotope and antigen for fcs markers

m <- match(fcs_colnames, clustering_observables$mass)

fcs_panel <- data.frame(colnames = fcs_colnames, Isotope = clustering_observables$mass[m], Antigen = clustering_observables$marker[m], stringsAsFactors = FALSE)


# -------------------------------------

cols <- which(fcs_colnames %in% clustering_observables$mass)

# arc-sin-h the columns specific 
fcsT <- lapply(fcs, function(u) {
  e <- exprs(u)
  e[ , cols] <- asinh( e[ , cols] / 5 )
  exprs(u) <- e
  u
})



# ------------------------------------------------------------

### Indeces of observables used for clustering 

scols <- which(fcs_colnames %in% clust_observ)

# ordered by decreasing pca score
scols <- scols[order(clustering_observables[fcs_colnames[scols], "avg_score"], decreasing = TRUE)]


### Indeces of other observables

xcols <- setdiff(cols, scols)

# ordered by decreasing pca score
xcols <- xcols[order(clustering_observables[fcs_colnames[xcols], "avg_score"], decreasing = TRUE)]




# ------------------------------------------------------------
# Get expression
# ------------------------------------------------------------

### Get marker expression

es <- lapply(fcsT, function(u) {
  exprs(u)[, cols]
})

e <- do.call("rbind",es)
colnames(e) <- fcs_panel$Antigen[cols]

a <- aggregate( e, by=list(clust), FUN=median)

mlab <- match(a$Group.1, labels$cluster)
rownames(a) <- labels$label[mlab]

colnames(a)[1] <- "cluster"

# normalize to 0-1
el <- e
rng <- apply(el,2,quantile,p=c(.01,.99))
for(i in 1:ncol(el)) {
  el[,i] <- (el[,i]-rng[1,i])/(rng[2,i]-rng[1,i])
}
el[el<0] <- 0
el[el>1] <- 1

al <- aggregate( el, by=list(clust), FUN=median)

mlab <- match(al$Group.1, labels$cluster)
rownames(al) <- labels$label[mlab]

colnames(al)[1] <- "cluster"


### Save cluster frequencies and the median expression

# get cluster frequencies
freq_clust <- table(clust)


# raw
clusters_out <- data.frame(cluster = names(freq_clust), label = labels[names(freq_clust), "label"], counts = as.numeric(freq_clust), frequencies = as.numeric(freq_clust)/sum(freq_clust), a[, fcs_panel$Antigen[c(scols, xcols)]])

write.table(clusters_out, file.path(hmDir, paste0(prefix, "clusters_raw.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# normalized 
clusters_out <- data.frame(cluster = names(freq_clust), label = labels[names(freq_clust), "label"], counts = as.numeric(freq_clust), frequencies = as.numeric(freq_clust)/sum(freq_clust), al[, fcs_panel$Antigen[c(scols, xcols)]])

write.table(clusters_out, file.path(hmDir, paste0(prefix, "clusters_norm.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


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
  # pdf(file.path(hmDir, paste0(prefix, "distros_raw_in.pdf")), w = ncol(e), h = nrow(labels))
  # print(ggp)
  # dev.off()
  
  ## with free scales using facet_wrap 
  dfm$variable_clust <- interaction(dfm$variable, dfm$clust, lex.order = FALSE)
  
  ggp <- ggplot(dfm, aes(x=value)) +
    geom_density(adjust=1, fill = "black", alpha = 0.3) +
    facet_wrap(~ variable_clust, nrow = nlevels(dfm$clust), scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text = element_text(size = 6), strip.text = element_text(size = 9))
  
  pdf(file.path(hmDir, paste0(prefix, "distrosfree", suffix, ".pdf")), w = ncol(e)*3/2, h = nrow(labels)*3/2)
  print(ggp)
  dev.off()
  
  return(NULL)
  
}


## Raw expression, included observables

plotting_wrapper(e = e[, fcs_panel$Antigen[scols]], suffix = "_raw_in")

# Normalized expression, included observables

plotting_wrapper(e = el[, fcs_panel$Antigen[scols]], suffix = "_norm_in")


if(length(xcols) > 0){
  
  # Normalized expression, excluded observables
  
  plotting_wrapper(e = e[, fcs_panel$Antigen[xcols]], suffix = "_norm_ex")
  
  # Raw expression, excluded observables
  
  plotting_wrapper(e = el[, fcs_panel$Antigen[xcols]], suffix = "_raw_ex")
  
}


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
  
  pdf(file.path(hmDir, paste0(prefix, "distrosmer", suffix,".pdf")), w = ncol(e), h = 10)
  print(ggp)
  dev.off()
  
  ## create colors per group not per sample
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
  
  pdf(file.path(hmDir, paste0(prefix, "distrosgrp", suffix,".pdf")), w = ncol(e), h = 10)
  print(ggp)
  dev.off()
  
  
  return(NULL)
  
}




## Raw expression, included observables

# plotting_wrapper2(e = e[, fcs_panel$Antigen[scols]], suffix = "_raw_in")
# 
# 
# # Normalized expression, included observables
# 
# plotting_wrapper2(e = el[, fcs_panel$Antigen[scols]], suffix = "_norm_in")
# 
# 
# if(length(xcols) > 0){
#   
#   # Raw expression, excluded observables
#   
#   plotting_wrapper2(e = e[, fcs_panel$Antigen[xcols]], suffix = "_raw_ex")
#   
#   
#   # Normalized expression, excluded observables
#   
#   plotting_wrapper2(e = el[, fcs_panel$Antigen[xcols]], suffix = "_norm_ex")
#   
# }


# ------------------------------------------------------------
# pheatmaps for raw data
# ------------------------------------------------------------

### Cluster clustering is based on all markers here

expr <- as.matrix(a[, fcs_panel$Antigen[c(scols, xcols)]])
rownames(expr) <- labels$label

cluster_rows <- hclust(dist(expr), method = "average")


labels_row <- paste0(as.character(labels$label), " (", round(as.numeric(freq_clust)/sum(freq_clust)*100, 2), "%)")
labels_col <- colnames(expr)

pheatmap(expr, color = colorRampPalette(brewer.pal(n = 8, name = "YlGnBu"))(100), cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, display_numbers = TRUE, number_color = "black", fontsize_number = 7, gaps_col = length(scols), fontsize_row = 10, fontsize_col = 10, fontsize = 7, filename = file.path(hmDir, paste0(prefix, "pheatmap_raw_row_clust.pdf")), width = 10, height = 7)


pheatmap(expr, color = colorRampPalette(brewer.pal(n = 8, name = "YlGnBu"))(100), cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, display_numbers = TRUE, number_color = "black", fontsize_number = 7, gaps_col = length(scols), fontsize_row = 10, fontsize_col = 10, fontsize = 7, filename = file.path(hmDir, paste0(prefix, "pheatmap_raw.pdf")), width = 10, height = 7)



# ------------------------------------------------------------
# pheatmaps for normalized data
# ------------------------------------------------------------

### Cluster clustering is based on all markers here

expr <- as.matrix(al[, fcs_panel$Antigen[c(scols, xcols)]])
rownames(expr) <- labels$label

cluster_rows <- hclust(dist(expr), method = "average")


labels_row <- paste0(as.character(labels$label), " (", round(as.numeric(freq_clust)/sum(freq_clust)*100, 2), "%)")
labels_col <- colnames(expr)

pheatmap(expr, color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(100), cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, breaks = seq(from = 0, to = 1, length.out = 101), legend_breaks = seq(from = 0, to = 1, by = 0.2), display_numbers = TRUE, number_color = "black", fontsize_number = 7, gaps_col = length(scols), fontsize_row = 10, fontsize_col = 10, fontsize = 7, filename = file.path(hmDir, paste0(prefix, "pheatmap_norm_row_clust.pdf")), width = 10, height = 7)


pheatmap(expr, color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(100), cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, breaks = seq(from = 0, to = 1, length.out = 101), legend_breaks = seq(from = 0, to = 1, by = 0.2), display_numbers = TRUE, number_color = "black", fontsize_number = 7, gaps_col = length(scols), fontsize_row = 10, fontsize_col = 10, fontsize = 7, filename = file.path(hmDir, paste0(prefix, "pheatmap_norm.pdf")), width = 10, height = 7)





# ------------------------------------------------------------
# ggplot tile
# ------------------------------------------------------------

# ggdf <- data.frame(cluster = labels$label, expr)

# ggdf$cluster <- factor(ggdf$cluster, levels = rev(levels(labels$label)))

# ggdfm <- melt(ggdf, id.var = "cluster", variable.name = "observable")
# ggdfm$observable <- factor(ggdfm$observable, levels = colnames(ggdf)[-1])


# ggp <- ggplot(ggdfm, aes(x = observable, y = cluster, fill = value)) + 
#   geom_tile() + 
#   geom_text(aes(label = round(value, 2)), color = "black", size = 2.5) + 
#   scale_x_discrete(expand = c(0, 0)) + 
#   scale_y_discrete(expand = c(0, 0)) + 
#   xlab("") + 
#   ylab("") + 
#   theme_bw() +
#   scale_fill_gradientn("", colours = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(100), limits=c(0, 1)) +
#   theme(panel.background = element_rect(fill = NA, colour = NA), axis.ticks = element_blank(), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size = 14)) 


# pdf(file.path(hmDir, paste0(prefix, "ggtile.pdf")), width = 10, height = 7)
# print(ggp)
# dev.off()



# ### Row order from clustering

# ggdfm$cluster <- factor(ggdfm$cluster, levels = rev(levels(labels$label)[cluster_rows$order]))

# ggp <- ggplot(ggdfm, aes(x = observable, y = cluster, fill = value)) + 
#   geom_tile() + 
#   geom_text(aes(label = round(value, 2)), color = "black", size = 2.5) + 
#   scale_x_discrete(expand = c(0, 0)) + 
#   scale_y_discrete(expand = c(0, 0)) + 
#   xlab("") + 
#   ylab("") + 
#   theme_bw() +
#   scale_fill_gradientn("", colours = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(100), limits=c(0, 1)) +
#   theme(panel.background = element_rect(fill = NA, colour = NA), axis.ticks = element_blank(), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size = 14)) 


# pdf(file.path(hmDir, paste0(prefix, "ggtile_row_clust.pdf")), width = 10, height = 7)
# print(ggp)
# dev.off()



# ------------------------------------------------------------
# ggdendogram
# ------------------------------------------------------------


# ggp <- ggdendrogram(cluster_rows, rotate = FALSE, size = 2)


# pdf(file.path(hmDir, paste0(prefix, "ggdendro_row_clust.pdf")), width = 7, height = 3)
# print(ggp)
# dev.off()


# dhc <- as.dendrogram(cluster_rows)
# ddata <- dendro_data(dhc, type = "rectangle")
# 
# 
# ggp <- ggplot(segment(ddata)) +
#   geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
#   geom_text(data = label(ddata), aes(x=x, y=y, label=label, hjust=0), size=3) +
#   coord_flip() +
#   scale_y_reverse(expand=c(0.2, 0)) +
#   theme_dendro()
# 
# 
# 
# pdf(file.path(hmDir, paste0(prefix, "ggdendro_row_clust.pdf")), width = 2, height = 7)
# print(ggp)
# dev.off()

# geom_segment(data=segment(hcdata), aes(x=x, y=y, xend=xend, yend=yend)) +
#   geom_text(data=label(hcdata), aes(x=x, y=y, label=label, hjust=0), size=3) +
#   coord_flip() + 
#   scale_y_reverse(expand=c(0.2, 0))
# 




# ------------------------------------------------------------
# ConsensusClusterPlus
# ------------------------------------------------------------

# rand_seed <- 1234

# pdf(file.path(hmDir, paste0(prefix, "ConsensusClusterPlus.pdf")), width = 7, height = 10)

# results <- ConsensusClusterPlus(d = t(as.matrix(expr)), maxK = 10, reps = 100, pItem = 0.9, pFeature = 1, plot = NULL, verbose = FALSE, clusterAlg = "hc", innerLinkage="average", finalLinkage="average", distance = "euclidean", seed = rand_seed)

# dev.off()


# pdf(file.path(hmDir, paste0(prefix, "ConsensusClusterPlus_pheatmap.pdf")), width = 10, height = 7)

# for(k in 2:10){
#   # k <- 5

#   fm <- results[[k]]$ml
#   cluster_rows_c <- hclust(as.dist( 1 - fm ), method = "average")

#   pheatmap(mat = expr, color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(100), cluster_cols = FALSE, cluster_rows = cluster_rows_c, labels_col = labels_col, labels_row = labels_row, border_color = NA, breaks = seq(from = 0, to = 1, length.out = 101), legend_breaks = seq(from = 0, to = 1, by = 0.2), display_numbers = TRUE, number_color = "black", fontsize_number = 7, gaps_col = length(scols), fontsize_row = 10, fontsize_col = 10, fontsize = 7, filename = NA, main = paste0("Consensus cluster number k = ", k))

# }

# dev.off()



# pdf(file.path(hmDir, paste0(prefix, "ConsensusClusterPlus_ggdendro.pdf")), width = 7, height = 3)

# for(k in 2:10){
#   # k <- 5

#   fm <- results[[k]]$ml
#   cluster_rows_c <- hclust(as.dist( 1 - fm ), method = "average")

#   ggp <- ggdendrogram(cluster_rows_c, rotate = FALSE, size = 2) +
#   ggtitle(paste0("Consensus cluster number k = ", k))

#   print(ggp)

# }

# dev.off()





sessionInfo()













################################
### 02_heatmaps done!
################################