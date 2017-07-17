

Sys.time()

# Load packages
library(gdata)
library(pheatmap)
library(RColorBrewer)


##############################################################################
# Test arguments
##############################################################################


prefix='23CD8allall_29CD8allall_02CD8v2_cl49_top10_glmer_binomial_interglht_'
outdir='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/02v2_CD8/090_cytokine_bimatrix_frequencies_clustering/3responses_both'
path_data='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/02v2_CD8/090_cytokine_bimatrix/23CD8allall_29CD8allall_02CD8v2_bimatrix.rds'
path_clustering_observables='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/02v2_CD8/090_cytokine_bimatrix/23CD8allall_29CD8allall_02CD8v2_clustering_observables.xls'
path_clustering='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/02v2_CD8/090_cytokine_bimatrix_frequencies_clustering/23CD8allall_29CD8allall_02CD8v2_cl49_clustering.xls'
path_clustering_labels='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/02v2_CD8/090_cytokine_bimatrix_frequencies_clustering/23CD8allall_29CD8allall_02CD8v2_cl49_clustering_labels.xls'
aggregate_fun='mean'
path_pvs='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/02v2_CD8/090_cytokine_bimatrix_frequencies_clustering/3responses_both/23CD8allall_29CD8allall_02CD8v2_cl49_frequencies_pvs_glmer_binomial_interglht_top10.xls'
FDR_cutoff='10'

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


if(!any(grepl("aggregate_fun=", args))){
  aggregate_fun='median'
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
labels

# clustering observables
clustering_observables <- read.table(path_clustering_observables, header = TRUE, sep = "\t", as.is = TRUE)
rownames(clustering_observables) <- clustering_observables$mass
clustering_observables

clust_observ <- clustering_observables[clustering_observables$clustering_observable, "mass"]
clust_observ

# ------------------------------------------------------------
# Load pvalues
# ------------------------------------------------------------

pvs <- read.table(path_pvs, header = TRUE, sep = "\t", as.is = TRUE)

comparisons <- colnames(pvs)[grep("adjp_", colnames(pvs))]
comparisons


# ------------------------------------------------------------
# Marker information
# ------------------------------------------------------------

# Get the isotope and antigen for fcs markers

m <- match(fcs_colnames, clustering_observables$mass)

fcs_panel <- data.frame(fcs_colname = fcs_colnames, Isotope = clustering_observables$mass[m], Antigen = clustering_observables$marker[m], stringsAsFactors = FALSE)


# Indeces of observables used for clustering 
scols <- which(fcs_colnames %in% clust_observ)
scols
# Indeces of other observables
xcols <- which(!fcs_colnames %in% clust_observ)
xcols

# Ordered by decreasing pca score
if("avg_score" %in% colnames(clustering_observables)){
  scols <- scols[order(clustering_observables[fcs_colnames[scols], "avg_score"], decreasing = TRUE)]
  xcols <- xcols[order(clustering_observables[fcs_colnames[xcols], "avg_score"], decreasing = TRUE)]
}


# ------------------------------------------------------------
# Get the median expression
# ------------------------------------------------------------

samp <- expr[, "sample_id"]
clust <- clustering[, "cluster"]
e <- expr[, fcs_colnames]


colnames(e) <- fcs_panel$Antigen

a <- aggregate(e, by = list(clust), FUN = aggregate_fun)

# get cluster frequencies
freq_clust <- table(clust)


# ------------------------------------------------------------
# Heatmaps of raw median expression
# ------------------------------------------------------------

### Use all markers for plotting
expr <- as.matrix(a[, fcs_panel$Antigen[c(scols, xcols)]])
rownames(expr) <- labels$label

# labels_row <- paste0(rownames(expr), " (", round(as.numeric(freq_clust)/sum(freq_clust)*100, 2), "%)")
labels_row <- as.character(labels$label)
names(labels_row) <- as.character(labels$label)
labels_col <- colnames(expr)


color <- colorRampPalette(brewer.pal(n = 8, name = 'Greys'))(100)



for(i in 1:length(comparisons)){
  # i = 1
  
  comparison <- comparisons[i]
  print(comparison)
  comparison_prefix <- paste0(gsub("adjp_", "", comparison), "_")
  
  pvs_sign <- pvs[pvs[, comparison] < FDR_cutoff & !is.na(pvs[, comparison]), , drop = FALSE]
  print(pvs_sign)
  
  prefix_sign <- paste0(prefix, comparison_prefix)
  
  if(nrow(pvs_sign) == 0){
    
    pdf(file.path(outdir, paste0(prefix_sign, "pheatmap.pdf")))
    plot(1, type="n", axes=F, xlab="", ylab="")
    dev.off()
    
  }else{
    
    pvs_sign <- pvs_sign[order(pvs_sign[, gsub("adjp_", "pval_", comparison)]), , drop = FALSE]
    
    rows_order <- as.character(pvs_sign$label)
    
    ## No row clustering
    pheatmap(expr[rows_order, , drop = FALSE], color = color, cellwidth = 24, cellheight = 24, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row[rows_order], display_numbers = FALSE, number_color = "black", fontsize_number = 8, gaps_col = length(scols), fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix_sign, "pheatmap.pdf")))
    
    
  }
  
}















sessionInfo()


