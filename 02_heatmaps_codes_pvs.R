

Sys.time()

# Load packages
library(FlowSOM)
library(ConsensusClusterPlus)
library(gdata)
library(ComplexHeatmap)
library(RColorBrewer)


##############################################################################
# Test arguments
##############################################################################


prefix='23_03_pca1_cl20_merging4_'
outdir='../carsten_cytof/PD1_project/CK_2016-06-23_03/030_codes'
path_data='../carsten_cytof/PD1_project/CK_2016-06-23_03/010_data/23_03_expr_raw.rds'
path_data_norm='../carsten_cytof/PD1_project/CK_2016-06-23_03/010_data/23_03_expr_norm.rds'
path_clustering_observables='../carsten_cytof/PD1_project/CK_2016-06-23_03/030_heatmaps/23_03_pca1_clustering_observables.xls'
path_codes_clustering='../carsten_cytof/PD1_project/CK_2016-06-23_03/030_heatmaps/23_03_pca1_cl20_codes_clustering.xls'
path_codes_clustering_labels='../carsten_cytof/PD1_project/CK_2016-06-23_03/030_heatmaps/23_03_pca1_cl20_codes_clustering_labels.xls'
path_marker_selection='../carsten_cytof/PD1_project/CK_2016-06-23_03/010_helpfiles/23_03_pca1_merging4_marker_selection_codes.txt'
path_cluster_merging='../carsten_cytof/PD1_project/CK_2016-06-23_03/010_helpfiles/23_03_pca1_cl20_cluster_merging4.xlsx'
path_fsom='../carsten_cytof/PD1_project/CK_2016-06-23_03/030_heatmaps/23_03_pca1_cl20_fsom.rds'
path_fccp='../carsten_cytof/PD1_project/CK_2016-06-23_03/030_heatmaps/23_03_pca1_cl20_fccp.rds'
path_pvs='../carsten_cytof/PD1_project/CK_2016-06-23_03/050_frequencies_codes/23_03_pca1_cl20_frequencies_pvs_glmer_binomial_interglht_top10.xls'
path_coeffs='../carsten_cytof/PD1_project/CK_2016-06-23_03/050_frequencies_codes/23_03_pca1_cl20_frequencies_coeffs_glmer_binomial_interglht_top10.xls'
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


linkage <- "average"

pheatmap_palette <- 'YlGnBu'
pheatmap_palette_rev <- FALSE

pheatmap_palette_norm <- 'RdYlBu'
pheatmap_palette_norm_rev <- TRUE


if(!any(grepl("aggregate_fun=", args))){
  aggregate_fun='median'
}

if(!any(grepl("scale=", args))){
  scale=TRUE
}

suffix <- paste0("_top", FDR_cutoff)
FDR_cutoff <- as.numeric(paste0("0.", FDR_cutoff))
FDR_cutoff


# ------------------------------------------------------------
# Load expression data
# ------------------------------------------------------------

expr <- readRDS(path_data)

cell_id <- expr[, "cell_id"]
samp <- expr[, "sample_id"]
fcs_colnames <- colnames(expr)[!grepl("cell_id|sample_id", colnames(expr))]

e <- expr[, fcs_colnames]

if(!is.null(path_data_norm)){
  expr_norm <- readRDS(path_data_norm)
  e_norm <- expr_norm[, fcs_colnames]
}

# ------------------------------------------------------------
# Load clustering data
# ------------------------------------------------------------

# clustering
clustering <- read.table(path_codes_clustering, header = TRUE, sep = "\t", as.is = TRUE)

clust <- clustering[, "cluster"]
names(clust) <- clustering[, "cell_id"]

# clustering labels
labels <- read.table(path_codes_clustering_labels, header = TRUE, sep = "\t", as.is = TRUE)

labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))
rownames(labels) <- labels$cluster
labels

# clustering observables
clustering_observables <- read.table(path_clustering_observables, header = TRUE, sep = "\t", as.is = TRUE)
rownames(clustering_observables) <- clustering_observables$mass
clustering_observables

clust_observ <- clustering_observables[clustering_observables$clustering_observable, "mass"]
clust_observ


# ------------------------------------------------------------
# Load FlowSOM objects
# ------------------------------------------------------------

fsom <- readRDS(path_fsom)

fccp <- readRDS(path_fccp)

k <- nmetaclusts <- length(fccp)

codes <- fsom$codes
ncodes <- nrow(codes)
rownames(codes) <- 1:ncodes

fsom_mc <- fccp[[k]]$consensusClass

if(!length(fsom_mc) == ncodes)
  stop("Some of the codes have zero cells!")

fsom_mc_tree <- fccp[[k]]$consensusTree


# ------------------------------------------------------------
# Load pvalues
# ------------------------------------------------------------

pvs <- read.table(path_pvs, header = TRUE, sep = "\t", as.is = TRUE)

pvs <- pvs[order(pvs$cluster), ]

comparisons <- colnames(pvs)[grep("adjp_", colnames(pvs))]
comparisons

# comparisons <- "adjp_NRvsR"

# ------------------------------------------------------------
# Load coeffs
# ------------------------------------------------------------

coeffs <- read.table(path_coeffs, header = TRUE, sep = "\t", as.is = TRUE)

coeffs <- coeffs[order(coeffs$cluster), ]


### Check

stopifnot(all(pvs$label == coeffs$label))


# ------------------------------------------------------------
# Prepare a color annotation for heatmaps 
# ------------------------------------------------------------

# ggplot palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=60 , c=100)[1:n]
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
  
  ### Reorder the annotations so that merging is more to the left in the Figure
  annotation_colors <- annotation_colors[2:1]
  annotation_row <- annotation_row[, 2:1, drop = FALSE]
  
}



# ------------------------------------------------------------
# Load marker selection for plotting on the heatmaps
# ------------------------------------------------------------

marker_selection <- NULL

if(!is.null(path_marker_selection)){
  if(file.exists(path_marker_selection)){
    
    marker_selection <- read.table(file.path(path_marker_selection), header = TRUE, sep = "\t", as.is = TRUE)
    marker_selection <- marker_selection[, 1]
    
    if(!all(marker_selection %in% clustering_observables$marker))
      stop("Marker selection is wrong")
    
  }
}

marker_selection

# ------------------------------------------------------------
# Marker information
# ------------------------------------------------------------

# Get the isotope and antigen for fcs markers

m <- match(fcs_colnames, clustering_observables$mass)

fcs_panel <- data.frame(fcs_colname = fcs_colnames, Isotope = clustering_observables$mass[m], Antigen = clustering_observables$marker[m], stringsAsFactors = FALSE)


# Indeces of observables used for clustering 
scols <- which(fcs_colnames %in% clust_observ)

# Indeces of other observables
xcols <- which(!fcs_colnames %in% clust_observ)


# Ordered by decreasing pca score
if("avg_score" %in% colnames(clustering_observables)){
  scols <- scols[order(clustering_observables[fcs_colnames[scols], "avg_score"], decreasing = TRUE)]
  xcols <- xcols[order(clustering_observables[fcs_colnames[xcols], "avg_score"], decreasing = TRUE)]
}

smarkers <- fcs_panel$Antigen[scols]
smarkers
xmarkers <- fcs_panel$Antigen[xcols]
xmarkers

# ------------------------------------------------------------
# Plotting
# ------------------------------------------------------------


clust <- clustering[, "cluster"]


# ------------------------------------------------------------
# Get the median expression
# ------------------------------------------------------------

colnames(e) <- fcs_panel$Antigen

a <- aggregate(e, by = list(clust), FUN = aggregate_fun)

# get cluster frequencies
freq_clust <- table(clust)


### Save cluster frequencies and the median expression

clusters_out <- data.frame(cluster = names(freq_clust), label = labels[names(freq_clust), "label"], counts = as.numeric(freq_clust), frequencies = as.numeric(freq_clust)/sum(freq_clust), a[, fcs_panel$Antigen[c(scols, xcols)]])

write.table(clusters_out, file.path(outdir, paste0(prefix, "codes_cluster_median_expression_raw.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# ------------------------------------------------------------
# Row clustering from the fccp object
# ------------------------------------------------------------

### This clustering is used in all the heatmaps
cluster_rows <- fsom_mc_tree


# ------------------------------------------------------------
# Heatmaps of raw median expression
# ------------------------------------------------------------

### Use all markers for plotting
expr <- as.matrix(a[, c(smarkers, xmarkers)])
rownames(expr) <- 1:ncodes


labels_row <- paste0(as.character(1:ncodes), " (", as.numeric(freq_clust), ")")
labels_col <- colnames(expr)


if(pheatmap_palette_rev){
  color <- colorRampPalette(rev(brewer.pal(n = 8, name = pheatmap_palette)))(100)
}else{
  color <- colorRampPalette(brewer.pal(n = 8, name = pheatmap_palette))(100)
}



pvs_cut_tmp <- cut(c(0.01, 0.05, 0.1, 1), c(0, 0.01, 0.05, 0.1, 1))
pvs_cut_tmp

pvs_cut_tmp <- factor(c(paste0("up", pvs_cut_tmp), paste0("down", pvs_cut_tmp)), levels = c(paste0("up", pvs_cut_tmp), paste0("down", pvs_cut_tmp)))
pvs_cut_tmp


colors_pvs <- c(colorRampPalette(c("#00008b", "#f5f5f5"), space = "Lab")(4), colorRampPalette(c("#dc143c", "#f5f5f5"), space = "Lab")(4))
names(colors_pvs) <- levels(pvs_cut_tmp)



comparisons <- list("adjp_NRvsR", "adjp_NRvsR_base", "adjp_NRvsR_tx", c("adjp_NRvsR_base", "adjp_NRvsR_tx"))
comparison_suffixs <- c("NRvsR", "NRvsR_base", "NRvsR_tx", "NRvsR_baseANDtx")


for(i in 1:length(comparisons)){
  # i = 2
  
  comparison <- comparisons[[i]]
  print(comparison)
  comparison_suffix <- comparison_suffixs[i]
  
  pvs_discrete <- cut(as.numeric(as.matrix(pvs[, comparison, drop = FALSE])), c(0, 0.01, 0.05, 0.1, 1))
  sign_discrete <- sign(as.numeric(as.matrix(coeffs[, gsub("adjp_" , "", comparison), drop = FALSE])))
  sign_discrete <- ifelse(sign_discrete == 1, "up", "down")
  
  pvs_discrete <- paste0(sign_discrete, pvs_discrete)
  
  pvs_heat <- matrix(pvs_discrete, ncol = length(comparison), byrow = FALSE)
  colnames(pvs_heat) <- comparison
  rownames(pvs_heat) <- pvs$cluster
  head(pvs_heat)
  
  legend_breaks <- seq(from = floor(min(expr)), to = ceiling(max(expr)), by = 1)
  
  ## With row clustering
  
  ha <-  HeatmapAnnotation(df = annotation_row, col = annotation_colors, which = "row", width = unit(1.5, "cm"))
  
  ha_text = rowAnnotation(text = row_anno_text(rownames(expr), gp = gpar(fontsize = 10)))
  
  ht1s <- Heatmap(expr[, smarkers, drop = FALSE], name = "in", col = color, cluster_columns = FALSE, cluster_rows = cluster_rows, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"), show_row_names = FALSE, row_dend_width = unit(20, "mm"))
  
  ht1x <- Heatmap(expr[, xmarkers, drop = FALSE], name = "out", col = color, cluster_columns = FALSE, cluster_rows = cluster_rows, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"), show_row_names = FALSE)
  
  ht2 <- Heatmap(pvs_heat, name = "apvs", col = colors_pvs, cluster_columns = FALSE, cluster_rows = cluster_rows, row_dend_reorder = FALSE, show_row_names = FALSE, width = unit(1, "cm"))
  
  pdf(file.path(outdir, paste0(prefix, "ComplexHeatmap_codes_all_row_clust_raw_", comparison_suffix ,".pdf")), width = 12, height = 14)
  draw(ha + ht1s + ht1x + ht2 + ha_text, row_dend_side = "left")
  dev.off()
  
  
  ## No row clustering
  ha <-  HeatmapAnnotation(df = annotation_row[rows_order, , drop = FALSE], col = annotation_colors, which = "row", width = unit(1.5, "cm"))
  
  ha_text = rowAnnotation(text = row_anno_text(rownames(expr[rows_order, , drop = FALSE]), gp = gpar(fontsize = 10)))
  
  ht1s <- Heatmap(expr[rows_order, smarkers, drop = FALSE], name = "in", col = color, cluster_columns = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"), show_row_names = FALSE, row_dend_width = unit(20, "mm"))
  
  ht1x <- Heatmap(expr[rows_order, xmarkers, drop = FALSE], name = "out", col = color, cluster_columns = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"), show_row_names = FALSE)
  
  ht2 <- Heatmap(pvs_heat[rows_order, , drop = FALSE], name = "apvs", col = colors_pvs, cluster_columns = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, show_row_names = FALSE, width = unit(1, "cm"))
  
  pdf(file.path(outdir, paste0(prefix, "ComplexHeatmap_codes_all_no_clust_raw_", comparison_suffix ,".pdf")), width = 12, height = 14)
  draw(ha + ht1s + ht1x + ht2 + ha_text, row_dend_side = "left")
  dev.off()
  
  
  ## No row clustering + plot only the significant codes
  codes_sign <- rowSums(pvs[, comparison, drop = FALSE] < FDR_cutoff) > 0
  
  if(any(codes_sign)){
    
    hh <- sum(codes_sign) / 4 + 2
    
    ha <-  HeatmapAnnotation(df = annotation_row[rows_order, , drop = FALSE][codes_sign[rows_order], , drop = FALSE], col = annotation_colors, which = "row", width = unit(1.5, "cm"))
    
    ha_text = rowAnnotation(text = row_anno_text(rownames(expr[rows_order, , drop = FALSE][codes_sign[rows_order], , drop = FALSE]), gp = gpar(fontsize = 10)))
    
    ht1s <- Heatmap(expr[rows_order, smarkers, drop = FALSE][codes_sign[rows_order], , drop = FALSE], name = "in", col = color, cluster_columns = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"), show_row_names = FALSE, row_dend_width = unit(20, "mm"))
    
    ht1x <- Heatmap(expr[rows_order, xmarkers, drop = FALSE][codes_sign[rows_order], , drop = FALSE], name = "out", col = color, cluster_columns = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"), show_row_names = FALSE)
    
    ht2 <- Heatmap(pvs_heat[rows_order, , drop = FALSE][codes_sign[rows_order], , drop = FALSE], name = "apvs", col = colors_pvs, cluster_columns = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, show_row_names = FALSE, width = unit(1, "cm"))
    
    pdf(file.path(outdir, paste0(prefix, "ComplexHeatmap_codes_all_no_clust_raw_", comparison_suffix, suffix ,".pdf")), width = 12, height = hh)
    draw(ha + ht1s + ht1x + ht2 + ha_text, row_dend_side = "left")
    dev.off()
    
  }
  
  
  
  ## Plot only the selected markers
  if(!is.null(marker_selection)){
    
    expr_sub <- expr[, marker_selection, drop = FALSE]
    
  }
  
  
  if(scale){
    
    # ------------------------------------------------------------
    # Heatmaps of raw median expression scalled by marker (column) 
    # ------------------------------------------------------------
    
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
    
    ha <-  HeatmapAnnotation(df = annotation_row, col = annotation_colors, which = "row", width = unit(1.5, "cm"))
    
    ha_text = rowAnnotation(text = row_anno_text(rownames(expr_scaled), gp = gpar(fontsize = 10)))
    
    ht1s <- Heatmap(expr_scaled[, smarkers, drop = FALSE], name = "in", col = color, cluster_columns = FALSE, cluster_rows = cluster_rows, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"), show_row_names = FALSE, row_dend_width = unit(20, "mm"))
    
    ht1x <- Heatmap(expr_scaled[, xmarkers, drop = FALSE], name = "out", col = color, cluster_columns = FALSE, cluster_rows = cluster_rows, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"), show_row_names = FALSE)
    
    ht2 <- Heatmap(pvs_heat, name = "apvs", col = colors_pvs, cluster_columns = FALSE, cluster_rows = cluster_rows, row_dend_reorder = FALSE, show_row_names = FALSE, width = unit(1, "cm"))
    
    pdf(file.path(outdir, paste0(prefix, "ComplexHeatmap_codes_all_row_clust_scale_", comparison_suffix ,".pdf")), width = 12, height = 14)
    draw(ha + ht1s + ht1x + ht2 + ha_text, row_dend_side = "left")
    dev.off()
    
    ## No row clustering
    ha <-  HeatmapAnnotation(df = annotation_row[rows_order, , drop = FALSE], col = annotation_colors, which = "row", width = unit(1.5, "cm"))
    
    ha_text = rowAnnotation(text = row_anno_text(rownames(expr_scaled[rows_order, , drop = FALSE]), gp = gpar(fontsize = 10)))
    
    ht1s <- Heatmap(expr_scaled[rows_order, smarkers, drop = FALSE], name = "in", col = color, cluster_columns = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"), show_row_names = FALSE, row_dend_width = unit(20, "mm"))
    
    ht1x <- Heatmap(expr_scaled[rows_order, xmarkers, drop = FALSE], name = "out", col = color, cluster_columns = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"), show_row_names = FALSE)
    
    ht2 <- Heatmap(pvs_heat[rows_order, , drop = FALSE], name = "apvs", col = colors_pvs, cluster_columns = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, show_row_names = FALSE, width = unit(1, "cm"))
    
    pdf(file.path(outdir, paste0(prefix, "ComplexHeatmap_codes_all_no_clust_scale_", comparison_suffix ,".pdf")), width = 12, height = 14)
    draw(ha + ht1s + ht1x + ht2 + ha_text, row_dend_side = "left")
    dev.off()
    
    
    
    ## No row clustering + plot only the significant codes
    codes_sign <- rowSums(pvs[, comparison, drop = FALSE] < FDR_cutoff) > 0
    
    if(any(codes_sign)){
      
      hh <- sum(codes_sign) / 4 + 2
      
      ha <-  HeatmapAnnotation(df = annotation_row[rows_order, , drop = FALSE][codes_sign[rows_order], , drop = FALSE], col = annotation_colors, which = "row", width = unit(1.5, "cm"))
      
      ha_text = rowAnnotation(text = row_anno_text(rownames(expr_scaled[rows_order, , drop = FALSE][codes_sign[rows_order], , drop = FALSE]), gp = gpar(fontsize = 10)))
      
      ht1s <- Heatmap(expr_scaled[rows_order, smarkers, drop = FALSE][codes_sign[rows_order], , drop = FALSE], name = "in", col = color, cluster_columns = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"), show_row_names = FALSE, row_dend_width = unit(20, "mm"))
      
      ht1x <- Heatmap(expr_scaled[rows_order, xmarkers, drop = FALSE][codes_sign[rows_order], , drop = FALSE], name = "out", col = color, cluster_columns = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"), show_row_names = FALSE)
      
      ht2 <- Heatmap(pvs_heat[rows_order, , drop = FALSE][codes_sign[rows_order], , drop = FALSE], name = "apvs", col = colors_pvs, cluster_columns = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, show_row_names = FALSE, width = unit(1, "cm"))
      
      pdf(file.path(outdir, paste0(prefix, "ComplexHeatmap_codes_all_no_clust_scale_", comparison_suffix , suffix,".pdf")), width = 12, height = hh)
      draw(ha + ht1s + ht1x + ht2 + ha_text, row_dend_side = "left")
      dev.off()
      
    }
    
    
    ## Plot only the selected markers
    if(!is.null(marker_selection)){
      
      ## No row clustering
      ha <-  HeatmapAnnotation(df = annotation_row[rows_order, , drop = FALSE], col = annotation_colors, which = "row", width = unit(1.5, "cm"))
      
      ha_text = rowAnnotation(text = row_anno_text(rownames(expr_scaled[rows_order, , drop = FALSE]), gp = gpar(fontsize = 10)))
      
      ht1s <- Heatmap(expr_scaled[rows_order, intersect(marker_selection, smarkers), drop = FALSE], name = "in", col = color, cluster_columns = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"), show_row_names = FALSE, row_dend_width = unit(20, "mm"))
      
      ht1x <- Heatmap(expr_scaled[rows_order, intersect(marker_selection, xmarkers), drop = FALSE], name = "out", col = color, cluster_columns = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"), show_row_names = FALSE)
      
      ht2 <- Heatmap(pvs_heat[rows_order, , drop = FALSE], name = "apvs", col = colors_pvs, cluster_columns = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, show_row_names = FALSE, width = unit(1, "cm"))
      
      pdf(file.path(outdir, paste0(prefix, "ComplexHeatmap_codes_sel_no_clust_scale_", comparison_suffix ,".pdf")), width = 12, height = 14)
      draw(ha + ht1s + ht1x + ht2 + ha_text, row_dend_side = "left")
      dev.off()
      
      
    }
    
    
  }
  
  
  if(!is.null(path_data_norm)){
    
    # ------------------------------------------------------------
    # Heatmaps of norm median expression
    # Had to do this way because I want to plot the 01 normalized data, but I want to keep row clustering from the raw data
    # ------------------------------------------------------------
    
    
    # ------------------------------------------------------------
    # Get the median expression
    # ------------------------------------------------------------
    
    colnames(e_norm) <- fcs_panel$Antigen
    
    a_norm <- aggregate(e_norm, by = list(clust), FUN = aggregate_fun)
    
    # ------------------------------------------------------------
    # pheatmaps of median expression
    # ------------------------------------------------------------
    
    ### Use all markers for plotting
    expr_norm <- as.matrix(a_norm[, fcs_panel$Antigen[c(scols, xcols)]])
    rownames(expr_norm) <- 1:ncodes
    
    labels_row <- paste0(as.character(1:ncodes), " (", as.numeric(freq_clust), ")")
    labels_col <- colnames(expr_norm)
    
    
    if(pheatmap_palette_norm_rev){
      color <- colorRampPalette(rev(brewer.pal(n = 8, name = pheatmap_palette_norm)))(101)
    }else{
      color <- colorRampPalette(brewer.pal(n = 8, name = pheatmap_palette_norm))(101)
    }
    
    ### Fixed legend range from 0 to 1 
    breaks = seq(from = 0, to = 1, length.out = 101)
    legend_breaks = seq(from = 0, to = 1, by = 0.2)
    
    
    ## With row clustering
    
    ha <-  HeatmapAnnotation(df = annotation_row, col = annotation_colors, which = "row", width = unit(1.5, "cm"))
    
    ha_text = rowAnnotation(text = row_anno_text(rownames(expr_norm), gp = gpar(fontsize = 10)))
    
    ht1s <- Heatmap(expr_norm[, smarkers, drop = FALSE], name = "in", col = color, cluster_columns = FALSE, cluster_rows = cluster_rows, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"), show_row_names = FALSE, row_dend_width = unit(20, "mm"))
    
    ht1x <- Heatmap(expr_norm[, xmarkers, drop = FALSE], name = "out", col = color, cluster_columns = FALSE, cluster_rows = cluster_rows, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"), show_row_names = FALSE)
    
    ht2 <- Heatmap(pvs_heat, name = "apvs", col = colors_pvs, cluster_columns = FALSE, cluster_rows = cluster_rows, row_dend_reorder = FALSE, show_row_names = FALSE, width = unit(1, "cm"))
    
    pdf(file.path(outdir, paste0(prefix, "ComplexHeatmap_codes_all_row_clust_norm_", comparison_suffix ,".pdf")), width = 12, height = 14)
    draw(ha + ht1s + ht1x + ht2 + ha_text, row_dend_side = "left")
    dev.off()
    
    ## No row clustering
    ha <-  HeatmapAnnotation(df = annotation_row[rows_order, , drop = FALSE], col = annotation_colors, which = "row", width = unit(1.5, "cm"))
    
    ha_text = rowAnnotation(text = row_anno_text(rownames(expr_norm[rows_order, , drop = FALSE]), gp = gpar(fontsize = 10)))
    
    ht1s <- Heatmap(expr_norm[rows_order, smarkers, drop = FALSE], name = "in", col = color, cluster_columns = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"), show_row_names = FALSE, row_dend_width = unit(20, "mm"))
    
    ht1x <- Heatmap(expr_norm[rows_order, xmarkers, drop = FALSE], name = "out", col = color, cluster_columns = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"), show_row_names = FALSE)
    
    ht2 <- Heatmap(pvs_heat[rows_order, , drop = FALSE], name = "apvs", col = colors_pvs, cluster_columns = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, show_row_names = FALSE, width = unit(1, "cm"))
    
    pdf(file.path(outdir, paste0(prefix, "ComplexHeatmap_codes_all_no_clust_norm_", comparison_suffix ,".pdf")), width = 12, height = 14)
    draw(ha + ht1s + ht1x + ht2 + ha_text, row_dend_side = "left")
    dev.off()
    
    
    ## No row clustering + plot only the significant codes
    codes_sign <- rowSums(pvs[, comparison, drop = FALSE] < FDR_cutoff) > 0
    
    if(any(codes_sign)){
      
      hh <- sum(codes_sign) / 4 + 2
      
      ha <-  HeatmapAnnotation(df = annotation_row[rows_order, , drop = FALSE][codes_sign[rows_order], , drop = FALSE], col = annotation_colors, which = "row", width = unit(1.5, "cm"))
      
      ha_text = rowAnnotation(text = row_anno_text(rownames(expr_norm[rows_order, , drop = FALSE][codes_sign[rows_order], , drop = FALSE]), gp = gpar(fontsize = 10)))
      
      ht1s <- Heatmap(expr_norm[rows_order, smarkers, drop = FALSE][codes_sign[rows_order], , drop = FALSE], name = "in", col = color, cluster_columns = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"), show_row_names = FALSE, row_dend_width = unit(20, "mm"))
      
      ht1x <- Heatmap(expr_norm[rows_order, xmarkers, drop = FALSE][codes_sign[rows_order], , drop = FALSE], name = "out", col = color, cluster_columns = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"), show_row_names = FALSE)
      
      ht2 <- Heatmap(pvs_heat[rows_order, , drop = FALSE][codes_sign[rows_order], , drop = FALSE], name = "apvs", col = colors_pvs, cluster_columns = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, show_row_names = FALSE, width = unit(1, "cm"))
      
      pdf(file.path(outdir, paste0(prefix, "ComplexHeatmap_codes_all_no_clust_norm_", comparison_suffix , suffix, ".pdf")), width = 12, height = hh)
      draw(ha + ht1s + ht1x + ht2 + ha_text, row_dend_side = "left")
      dev.off()
      
    }
    
    
    
    ## Plot only the selected markers
    if(!is.null(marker_selection)){
      
      expr_sub <- expr_norm[, marker_selection, drop = FALSE]
      
    }
    
    
  }
  
  
  # ------------------------------------------------------------
  # Heatmaps of the original codes
  # ------------------------------------------------------------
  
  
  ### Use the code markers for plotting
  expr_codes <- codes
  rownames(expr_codes) <- 1:ncodes
  
  mm <- match(colnames(codes), fcs_panel$fcs_colname)
  colnames(expr_codes) <- fcs_panel$Antigen[mm]
  expr_codes <- expr_codes[, fcs_panel$Antigen[scols]]  
  
  labels_row <- paste0(as.character(1:ncodes), " (", as.numeric(freq_clust), ")")
  labels_col <- colnames(expr_codes)
  
  color <- colorRampPalette(rev(brewer.pal(n = 8, name = "Spectral")))(101)
  
  
  ## With row clustering
  
  ha <-  HeatmapAnnotation(df = annotation_row, col = annotation_colors, which = "row", width = unit(1.5, "cm"))
  
  ha_text = rowAnnotation(text = row_anno_text(rownames(expr_codes), gp = gpar(fontsize = 10)))
  
  ht1s <- Heatmap(expr_codes, name = "in", col = color, cluster_columns = FALSE, cluster_rows = cluster_rows, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"), show_row_names = FALSE, row_dend_width = unit(20, "mm"))
  
  ht1x <- Heatmap(expr_codes, name = "out", col = color, cluster_columns = FALSE, cluster_rows = cluster_rows, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"), show_row_names = FALSE)
  
  ht2 <- Heatmap(pvs_heat, name = "apvs", col = colors_pvs, cluster_columns = FALSE, cluster_rows = cluster_rows, row_dend_reorder = FALSE, show_row_names = FALSE, width = unit(1, "cm"))
  
  pdf(file.path(outdir, paste0(prefix, "ComplexHeatmap_codes_codes_row_clust_raw_", comparison_suffix ,".pdf")), width = 12, height = 14)
  draw(ha + ht1s + ht1x + ht2 + ha_text, row_dend_side = "left")
  dev.off()
  
  ## No row clustering
  ha <-  HeatmapAnnotation(df = annotation_row[rows_order, , drop = FALSE], col = annotation_colors, which = "row", width = unit(1.5, "cm"))
  
  ha_text = rowAnnotation(text = row_anno_text(rownames(expr_codes[rows_order, , drop = FALSE]), gp = gpar(fontsize = 10)))
  
  ht1s <- Heatmap(expr_codes[rows_order, , drop = FALSE], name = "in", col = color, cluster_columns = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"), show_row_names = FALSE, row_dend_width = unit(20, "mm"))
  
  ht1x <- Heatmap(expr_codes[rows_order, , drop = FALSE], name = "out", col = color, cluster_columns = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"), show_row_names = FALSE)
  
  ht2 <- Heatmap(pvs_heat[rows_order, , drop = FALSE], name = "apvs", col = colors_pvs, cluster_columns = FALSE, cluster_rows = FALSE, row_dend_reorder = FALSE, show_row_names = FALSE, width = unit(1, "cm"))
  
  pdf(file.path(outdir, paste0(prefix, "ComplexHeatmap_codes_codes_no_clust_raw_", comparison_suffix ,".pdf")), width = 12, height = 14)
  draw(ha + ht1s + ht1x + ht2 + ha_text, row_dend_side = "left")
  dev.off()
  
  
}




















sessionInfo()


