##############################################################################
## <<02_heatmaps.R>>

# BioC 3.3
# Created 27 July 2016
# Updated 26 Jan 2017

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

rwd='/home/Shared/data/cytof/carsten_cytof/CK_2016-06-23_01'
heatmap_prefix='23_01_pca1_merging6_raw_'
heatmap_outdir='030_heatmaps'
path_data='010_data/23_01_expr_raw.rds'
path_metadata='/home/Shared/data/cytof/carsten_cytof/CK_metadata/metadata_23_01.xlsx'
path_clustering_observables='030_heatmaps/23_01_pca1_clustering_observables.xls'
path_clustering='030_heatmaps/23_01_pca1_merging6_clustering.xls'
path_clustering_labels='030_heatmaps/23_01_pca1_merging6_clustering_labels.xls'
path_marker_selection='23_01_pca1_merging6_marker_selection.txt'
aggregate_fun='median'
pheatmap_palette='YlGnBu'
pheatmap_palette_rev=FALSE
pheatmap_scale=TRUE
extra_path_data='010_data/23_01_expr_norm.rds'
extra_pheatmap_palette='Greys'
extra_pheatmap_palette_rev=FALSE
extra_suffix='_norm'

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

prefix <- heatmap_prefix
outdir <- heatmap_outdir

if(!file.exists(outdir)) 
  dir.create(outdir, recursive = TRUE)



if(all(!grepl("linkage=", args))){
  linkage <- "average"
}

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
# Load cluster merging data
# ------------------------------------------------------------
### Used to plot an annotation on the heatmap of how the mering looks like 

# elements for the heatmap
annotation_row <- NA
annotation_colors = NA
rows_order <- 1:nrow(labels)

if(any(grepl("path_cluster_merging=", args))){
  
  # read cluster merging file
  cm <- read.xls(path_cluster_merging)
  
  if(!all(c("old_cluster", "label", "new_cluster") %in% colnames(cm)))
    stop("Merging file must contain 'old_cluster', 'label' and 'new_cluster' columns!")
  
  # remove spaces in labels bcs they are problematic...
  cm$label <- factor(cm$label, labels = gsub(" ", "_", levels(cm$label))) 
  
  cm_unique <- unique(cm[, c("label", "new_cluster")])
  cm_unique <- cm_unique[order(cm_unique$new_cluster), ]
  
  cm$new_label <- factor(cm$label, levels = cm_unique$label)
  
  cmm <- merge(labels, cm[, c("old_cluster", "new_cluster", "new_label"), drop = FALSE], by.x = "cluster", by.y = "old_cluster", all.x = TRUE)
  
  
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
  color_ramp <- c(colors_muted, gg_color_hue(max(1, nlevels(cmm$new_label) - length(colors_muted))))
  
  colors_clusters <- color_ramp[1:nlevels(cmm$new_label)]
  names(colors_clusters) <- levels(cmm$new_label)
  
  
  # elements for the heatmap
  annotation_row <- data.frame(cluster = cmm$new_label)
  rownames(annotation_row) <- rownames(cmm$label)
  
  annotation_colors <- list(cluster = colors_clusters)
  
  rows_order <- order(cmm$new_label)
  
}


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


samp_org <- samp
e_org <- e
clust_org <- clust


subset_samp <- list()
subset_samp[[1]] <- unique(samp)
subset_samp[[2]] <- unique(samp)[grep("_HD", unique(samp))]


suffix_subset <- c("", "_HD")

if(any(grepl("merging_suffix=", args))){
  
  suffix_subset <- paste0(c("", "_HD"), merging_suffix)
  
}


### Plot heatmaps based on all the data or HD samples only
for(ii in 1:length(subset_samp)){
  # ii = 1
  
  cells2keep <- samp_org %in% subset_samp[[ii]]
  
  samp <- samp_org[cells2keep]
  e <- e_org[cells2keep, , drop = FALSE]
  clust <- clust_org[cells2keep]
  
  suffix <- suffix_subset[ii]
  
  # ------------------------------------------------------------
  # Get the median expression
  # ------------------------------------------------------------
  
  colnames(e) <- fcs_panel$Antigen
  
  a <- aggregate(e, by = list(clust), FUN = aggregate_fun)
  
  # get cluster frequencies
  freq_clust <- table(clust)
  
  
  ### Save cluster frequencies and the median expression
  
  clusters_out <- data.frame(cluster = names(freq_clust), label = labels[names(freq_clust), "label"], counts = as.numeric(freq_clust), frequencies = as.numeric(freq_clust)/sum(freq_clust), a[, fcs_panel$Antigen[c(scols, xcols)]])
  
  write.table(clusters_out, file.path(outdir, paste0(prefix, "clusters", suffix, ".xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  
  # ------------------------------------------------------------
  # pheatmaps of median expression
  # ------------------------------------------------------------
  
  ### This clustering is based on the markers that were used for the main clustering, and it is used in all the heatmaps
  expr <- as.matrix(a[, fcs_panel$Antigen[scols]])
  rownames(expr) <- labels$label
  
  cluster_rows <- hclust(dist(expr), method = linkage)
  
  
  ### Use all markers for plotting
  expr <- as.matrix(a[, fcs_panel$Antigen[c(scols, xcols)]])
  rownames(expr) <- labels$label
  
  labels_row <- paste0(as.character(labels$label), " (", round(as.numeric(freq_clust)/sum(freq_clust)*100, 2), "%)")
  labels_col <- colnames(expr)
  
  if(pheatmap_palette_rev){
    color <- colorRampPalette(rev(brewer.pal(n = 8, name = pheatmap_palette)))(100)
  }else{
    color <- colorRampPalette(brewer.pal(n = 8, name = pheatmap_palette))(100)
  }
  
  ## With row clustering
  pheatmap(expr, color = color, cellwidth = 24, cellheight = 24, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, display_numbers = TRUE, number_color = "black", fontsize_number = 8, gaps_col = length(scols), fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "pheatmap_all_row_clust", suffix, ".pdf")))
  
  ## No row clustering
  pheatmap(expr[rows_order, ], color = color, cellwidth = 24, cellheight = 24, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row[rows_order], display_numbers = TRUE, number_color = "black", fontsize_number = 8, gaps_col = length(scols), fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "pheatmap_all", suffix, ".pdf")))
  
  
  ## Plot only the selected markers
  if(!is.null(marker_selection)){
    expr_sub <- expr[, marker_selection, drop = FALSE]
    labels_col_sub <- colnames(expr_sub)
    
    pheatmap(expr_sub, color = color, cellwidth = 24, cellheight = 24, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col_sub, labels_row = labels_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "pheatmap_s1_row_clust", suffix, ".pdf")))
    
    pheatmap(expr_sub[rows_order, ], color = color, cellwidth = 24, cellheight = 24, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col_sub, labels_row = labels_row[rows_order], fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "pheatmap_s1", suffix, ".pdf")))
    
  }
  
  
  
  # ------------------------------------------------------------
  # pheatmaps of median expression scalled by marker (column) 
  # ------------------------------------------------------------
  
  if(pheatmap_scale){
    
    scalling_type <- "s01" 
    
    switch(scalling_type, 
      snorm = {
        ## scalled to mean = 0, sd = 1
        expr_scaled <- apply(expr, 2, function(x){(x-mean(x))/sd(x)})
        th <- 2.5
        expr_scaled[expr_scaled > th] <- th
        expr_scaled[expr_scaled < -th] <- -th
        breaks = seq(from = -th, to = th, length.out = 101)
        legend_breaks = seq(from = -round(th), to = round(th), by = 1)
      },
      
      s01 = {
        ## scalled to 01
        expr_scaled <- apply(expr, 2, function(x){(x-min(x))/(max(x)-min(x))})
        breaks = seq(from = 0, to = 1, length.out = 101)
        legend_breaks = seq(from = 0, to = 1, by = 0.25)
        
      }
    )
    
    
    color <- colorRampPalette(brewer.pal(n = 8, name = "Greys"))(120)[11:110]
    
    ## With row clustering
    pheatmap(expr_scaled, color = color, cellwidth = 24, cellheight = 24, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, breaks = breaks, legend_breaks = legend_breaks, display_numbers = TRUE, number_color = "black", fontsize_number = 8, gaps_col = length(scols), fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "pheatmap_all_row_clust_scale", suffix, ".pdf")))
    
    ## No row clustering
    pheatmap(expr_scaled[rows_order, ], color = color, cellwidth = 24, cellheight = 24, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row[rows_order], breaks = breaks, legend_breaks = legend_breaks, display_numbers = TRUE, number_color = "black", fontsize_number = 8, gaps_col = length(scols), fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "pheatmap_all_scale", suffix, ".pdf")))
    
    
    ## Plot only the selected markers
    if(!is.null(marker_selection)){
      expr_sub <- expr_scaled[, marker_selection, drop = FALSE]
      labels_col_sub <- colnames(expr_sub)
      
      pheatmap(expr_sub, color = color, cellwidth = 24, cellheight = 24, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col_sub, labels_row = labels_row, breaks = breaks, legend_breaks = legend_breaks, fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "pheatmap_s1_row_clust_scale", suffix, ".pdf")))
      
      pheatmap(expr_sub[rows_order, ], color = color, cellwidth = 24, cellheight = 24, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col_sub, labels_row = labels_row[rows_order], breaks = breaks, legend_breaks = legend_breaks, fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "pheatmap_s1_scale", suffix, ".pdf")))
      
    }
    
    
  }
  
  
  
  # ------------------------------------------------------------
  # Plot phetmaps based on the EXTRA data - works only if the extra arguments are specified!!!
  # Had to do this way because I want to plot the 01 normalized data, but I want to keep row clustering from the raw data
  # ------------------------------------------------------------
  
  
  if(any(grepl("extra_path_data=", args))){
    
    
    # ------------------------------------------------------------
    # Load extra data
    # ------------------------------------------------------------
    
    # read expression data
    if(grepl(".txt", extra_path_data)){
      expr_extra <- read.table(extra_path_data, header = TRUE, sep = "\t", as.is = TRUE)
    }
    if(grepl(".rds", extra_path_data)){
      expr_extra <- readRDS(extra_path_data)
    }
    
    
    # ------------------------------------------------------------
    # Keep expression and clustering results for the cells that are common in both
    # ------------------------------------------------------------
    
    common_cells <- intersect(clustering[, "cell_id"], expr_extra[, "cell_id"])
    
    samp.ex <- expr_extra[expr_extra[, "cell_id"] %in% common_cells, "sample_id"]
    e.ex <- expr_extra[expr_extra[, "cell_id"] %in% common_cells, fcs_colnames]
    clust.ex <- clustering[clustering[, "cell_id"] %in% common_cells, "cluster"]
    
    cells2keep <- samp.ex %in% subset_samp[[ii]]
    
    samp.ex <- samp.ex[cells2keep]
    e.ex <- e.ex[cells2keep, , drop = FALSE]
    clust.ex <- clust.ex[cells2keep]
    
    # ------------------------------------------------------------
    # Get the median expression
    # ------------------------------------------------------------
    
    colnames(e.ex) <- fcs_panel$Antigen
    
    a.ex <- aggregate(e.ex, by = list(clust.ex), FUN = aggregate_fun)
    
    # ------------------------------------------------------------
    # pheatmaps of median expression
    # ------------------------------------------------------------
    
    ### Use all markers for plotting
    expr <- as.matrix(a.ex[, fcs_panel$Antigen[c(scols, xcols)]])
    rownames(expr) <- labels$label
    
    labels_row <- paste0(as.character(labels$label), " (", round(as.numeric(freq_clust)/sum(freq_clust)*100, 2), "%)")
    labels_col <- colnames(expr)
    
    if(extra_pheatmap_palette_rev){
      color <- colorRampPalette(rev(brewer.pal(n = 8, name = extra_pheatmap_palette)))(130)[c(1, 32:130)]
    }else{
      color <- colorRampPalette(brewer.pal(n = 8, name = extra_pheatmap_palette))(130)[c(1, 32:130)]
    }
    
    ### Fixed legend range from 0 to 1 
    # breaks <- NA
    # legend_breaks <- NA
    breaks = seq(from = 0, to = 1, length.out = 101)
    legend_breaks = seq(from = 0, to = 1, by = 0.2)
    
    pheatmap(expr, color = color, cellwidth = 24, cellheight = 24, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, breaks = breaks, legend_breaks = legend_breaks, display_numbers = TRUE, number_color = "black", fontsize_number = 8, gaps_col = length(scols), fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "pheatmap_all_row_clust", extra_suffix, suffix, ".pdf")))
    
    pheatmap(expr[rows_order, ], color = color, cellwidth = 24, cellheight = 24, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row[rows_order], breaks = breaks, legend_breaks = legend_breaks, display_numbers = TRUE, number_color = "black", fontsize_number = 8, gaps_col = length(scols), fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "pheatmap_all", extra_suffix, suffix, ".pdf")))
    
    
    ## Plot only the selected markers
    if(!is.null(marker_selection)){
      expr_sub <- expr[, marker_selection, drop = FALSE]
      labels_col_sub <- colnames(expr_sub)
      
      pheatmap(expr_sub, color = color, cellwidth = 24, cellheight = 24, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col_sub, labels_row = labels_row, breaks = breaks, legend_breaks = legend_breaks, fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "pheatmap_s1_row_clust", extra_suffix, suffix, ".pdf")))
      
      pheatmap(expr_sub[rows_order, ], color = color, cellwidth = 24, cellheight = 24, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col_sub, labels_row = labels_row[rows_order], breaks = breaks, legend_breaks = legend_breaks, fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "pheatmap_s1", extra_suffix, suffix, ".pdf")))
      
    }
    
    
    
  }
  
  
  
  # ------------------------------------------------------------
  # Plot expression of markers in each cluster - distros plots
  # ------------------------------------------------------------
  
  
  plotting_wrapper <- function(e, suffix){
    
    df <- data.frame(e, clust = clust, check.names = FALSE)
    dfm <- melt(df, id.vars = c("clust"))
    
    dfm$clust <- factor(dfm$clust, levels = labels$cluster[cluster_rows$order], labels = paste0(as.character(labels$label), " (", round(as.numeric(freq_clust)/sum(freq_clust)*100, 2), "%)")[cluster_rows$order])
    
    ggp <- ggplot(dfm, aes(x=value, y = ..scaled..)) +
      geom_density(adjust=1, fill = "blue", alpha = 0.3) +
      facet_grid(clust ~ variable, scales = "free") +
      xlab("Marker expression") +
      ylab("Scaled density") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text = element_text(size = 8), strip.text = element_text(size = 9)) 
    
    pdf(file.path(outdir, paste0(prefix, "distros", suffix, ".pdf")), w = ncol(e)*1.8, h = nrow(labels)*1.6)
    print(ggp)
    dev.off()
    
    return(NULL)
    
  }
  
  
  # ## Expression, included observables
  # plotting_wrapper(e = e[, fcs_panel$Antigen[scols]], suffix = paste0("_in", suffix))
  # 
  # ## Expression, excluded observables
  # if(length(xcols) > 0){
  #   plotting_wrapper(e = e[, fcs_panel$Antigen[xcols]], suffix = paste0("_ex", suffix))
  # }
  
  
  
  
  
}









sessionInfo()













################################
### 02_heatmaps done!
################################