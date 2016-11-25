##############################################################################
## <<06_cytokines_bimatrix_top_selection_merged.R>>

# BioC 3.3
# Created 13 Nov 2016
# Updated 13 Nov 2016

# Make clusters out of the top frequent combinations of positive markers 
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
library(UpSetR)

##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-merged_23_29/02_CD4'
clust_prefix='23_29_CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_cl30_'
clust_outdir='10_cytokines_merged_top_combinations/01_clustering'

path_data=c('/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_02_CD4_merging2/060_cytokines_bimatrix/01_clustering/23CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_bimatrix.txt','/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-29_02_CD4_merging/060_cytokines_bimatrix/01_clustering/29CD4_02CD4_pca1_merging3_Tmem_cytCM_raw2_bimatrix.txt')
path_clustering_observables=c('/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_02_CD4_merging2/060_cytokines_bimatrix/01_clustering/23CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_clustering_observables.xls','/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-29_02_CD4_merging/060_cytokines_bimatrix/01_clustering/29CD4_02CD4_pca1_merging3_Tmem_cytCM_raw2_clustering_observables.xls')

data_name=c('data23','data29')
top_combinations=30


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

prefix <- clust_prefix
outdir <- clust_outdir

if(!file.exists(outdir)) 
  dir.create(outdir, recursive = TRUE)


# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------


expr <- lapply(1:length(data_name), function(i){
  
  if(grepl(".txt", path_data[i])){
    expr <- read.table(path_data[i], header = TRUE, sep = "\t", as.is = TRUE)
  }
  if(grepl(".rds", path_data[i])){
    expr <- readRDS(path_data[i])
  }
  
  expr$data_id <- data_name[i]
  return(expr)
})

expr <- rbind.fill(expr)

### cell_id is a paste of cell_id and sample_id because cell ids are the same for data 23 and 29 
cell_id <- paste0(expr[, "cell_id"])
samp <- expr[, "sample_id"]
fcs_colnames <- colnames(expr)[!grepl("cell_id|sample_id|data_id", colnames(expr))]
e <- expr[, fcs_colnames]

data_id <- paste0(expr[, "data_id"])


# ------------------------------------------------------------
# Load clustering_observables
# ------------------------------------------------------------
# Use an intersection of observables used in each dataset

clustering_observables <- lapply(1:length(data_name), function(i){
  
  clustering_observables <- read.table(path_clustering_observables[i], header = TRUE, sep = "\t", as.is = TRUE)
  return(clustering_observables)
  
})

clustering_observables <- Reduce(function(...) merge(..., by = c("mass", "marker"), all=TRUE, sort = FALSE), clustering_observables)

clustering_observables$clustering_observable <- apply(clustering_observables[, grep("clustering_observable", colnames(clustering_observables))], 1, all)

clust_observ <- clustering_observables[clustering_observables$clustering_observable, "mass"]




# selected columns for clustering 

scols <- which(fcs_colnames %in% clust_observ)

ef <- as.matrix(e[, scols])


# ------------------------------------------------------------
# Upsetr plots
# ------------------------------------------------------------

bidf <- data.frame(ef[, clustering_observables[clustering_observables$clustering_observable, "mass"]], row.names = 1:nrow(ef), check.names = FALSE)
colnames(bidf) <- clustering_observables[clustering_observables$clustering_observable, "marker"]

pdf(file.path(outdir, paste0(prefix, "upsetr.pdf")), w = 16, h = 6)
upset(bidf, sets = colnames(bidf), nintersects = 50, order.by = "freq")
dev.off()

# --------------------------------------------------------------------------
# Make clusters out of the top frequent combinations of positive markers 
# --------------------------------------------------------------------------


bivec <- apply(bidf, 1, paste0, collapse = "")

table_bivec <- table(bivec)

sort_freqs <- sort(table_bivec, decreasing = TRUE)
sort_comb <- names(sort_freqs)


# keep top_combinations and to the rest assign "00000000000" 
clust_comb <- factor(bivec, levels = sort_comb)
levels(clust_comb) <- c(sort_comb[1:top_combinations], rep(paste0(rep(0, ncol(bidf)), collapse = ""), length(sort_comb) - top_combinations))


clust <- as.numeric(clust_comb)


### Save clustering results

clust_out <- data.frame(cluster = clust, cell_id = cell_id, sample_id = samp, stringsAsFactors = FALSE)

clust_out <- split(clust_out, factor(data_id, levels = data_name))

for(i in 1:length(clust_out))
  write.table(clust_out[[i]], file = file.path(outdir, paste0(prefix, data_name[i], "_", "clustering.xls")), row.names=FALSE, quote=FALSE, sep="\t")




### Save data frames with labels

# get cluster frequencies
freq_clust <- table(clust)
prop_clust <- round(freq_clust/sum(freq_clust)*100, 2)

# get the names of cytokines that are positive in a given cluster

a <- aggregate(ef, by = list(clust), FUN = "mean")

data_new_labels <- a[, clustering_observables[clustering_observables$clustering_observable, "mass"]]
colnames(data_new_labels) <- clustering_observables[clustering_observables$clustering_observable, "marker"]
rownames(data_new_labels) <- a[, 1]

cluster_binames <- apply(data_new_labels, 1, function(x){
  
  paste0(rev(colnames(data_new_labels)[x > 0.75]), collapse = "/")
  
})

names(cluster_binames) <- a[, 1]

# the last cluster is a drop
cluster_binames[top_combinations + 1] <- "drop"



labels <- list()

for(i in 1:length(clust_out)){
  # i = 2
  
  clust_tmp <- clust_out[[i]]$cluster
  
  cluster <- as.character(sort(unique(clust_tmp)))
  
  label <- paste0(cluster, " - ", cluster_binames[cluster], " (", prop_clust[cluster], ")")
  
  # drop cluster must have only "drop" in the name
  label[grep("drop", label)] <- "drop"
  
  labels[[i]] <- data.frame(cluster = cluster, label = label, stringsAsFactors = FALSE)
  
  write.table(labels[[i]], file = file.path(outdir, paste0(prefix, data_name[i], "_", "clustering_labels.xls")), row.names=FALSE, quote=FALSE, sep="\t")
  
}


### Plot heatmap of generated clusters


color <- rev(colorRampPalette(brewer.pal(n = 8, name = "Spectral"))(100))
labels_col <- colnames(data_new_labels)

cluster_rows <- FALSE

cluster <- as.character(sort(unique(clust)))

labels_row <- paste0(cluster, " - ", cluster_binames[cluster], " (", prop_clust[cluster], ")")

annotation_row <- NA
annotation_colors <- NA

pheatmap(data_new_labels, color = color, cellwidth = 20, cellheight = 20, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, display_numbers = TRUE, number_color = "black", fontsize_number = 6, fontsize_row = 10, fontsize_col = 10, fontsize = 10, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "clusters_pheatmap.pdf")))



### Save cluster frequencies and the median expression

clusters_out <- data.frame(cluster = cluster, label = labels_row, counts = as.numeric(freq_clust), frequencies = as.numeric(freq_clust)/sum(freq_clust), data_new_labels)

write.table(clusters_out, file.path(outdir, paste0(prefix, "clusters.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)





# ---------------------------------------
# Compare different distance measures


# data_bi <- data_new_labels[1:top_combinations, ]
# 
# data_bi <- unique(bidf)
# 
# color <- rev(colorRampPalette(brewer.pal(n = 8, name = "Spectral"))(100))
# labels_col <- colnames(data_new_labels)
# 
# labels_row <- NULL
# 
# annotation_row <- NA
# annotation_colors <- NA
# 
# de <- as.matrix(dist(data_bi, method = "euclidean"))
# 
# cluster_rows <- hclust(dist(data_bi, method = "euclidean"), method = "average")
# 
# pheatmap(data_bi, color = color, cellwidth = 20, cellheight = 20, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, display_numbers = TRUE, number_color = "black", fontsize_number = 6, fontsize_row = 10, fontsize_col = 10, fontsize = 5, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "clusters_pheatmap_euclidean_average.pdf")))
# 
# 
# db <- as.matrix(dist(data_bi, method = "binary"))
# 
# cluster_rows <- hclust(dist(data_bi, method = "binary"), method = "average")
# 
# pheatmap(data_bi, color = color, cellwidth = 20, cellheight = 20, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, display_numbers = TRUE, number_color = "black", fontsize_number = 6, fontsize_row = 10, fontsize_col = 10, fontsize = 5, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "clusters_pheatmap_binary_average.pdf")))
# 
# 
# cluster_rows <- hclust(dist(data_bi, method = "binary"), method = "ward.D2")
# 
# pheatmap(data_bi, color = color, cellwidth = 20, cellheight = 20, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, display_numbers = TRUE, number_color = "black", fontsize_number = 6, fontsize_row = 10, fontsize_col = 10, fontsize = 5, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "clusters_pheatmap_binary_ward.pdf")))



# --------------------------------------------------------------------------
# Additionally, save the freq_out as counts.xls files which is necessary for the 08_frequencies_merged.R
# --------------------------------------------------------------------------


# ---------------------------------------
# Calculate the cluster frequencies per sample
# ---------------------------------------


for(i in 1:length(data_name)){
  
  sub_index <- data_id == data_name[i]
  
  # calculate frequencies
  freq <- table(clust[sub_index], samp[sub_index])
  
  prop <- t(t(freq) / colSums(freq)) * 100 # proportion of clusters in samples
  
  # use labels as names of clusters
  mlab <- match(rownames(freq), labels[[i]]$cluster)
  
  
  ### Save the frequencies and proportions
  prop_out <- data.frame(labels[[i]][mlab, c("cluster", "label")], as.data.frame.matrix(prop))
  
  freq_out <- data.frame(labels[[i]][mlab, c("cluster", "label")], as.data.frame.matrix(freq))
  
  write.table(prop_out, file=file.path(outdir, paste0(prefix, data_name[i], "_", "frequencies.xls")), row.names=FALSE, quote=FALSE, sep="\t")
  write.table(freq_out, file=file.path(outdir, paste0(prefix, data_name[i], "_", "counts.xls")), row.names=FALSE, quote=FALSE, sep="\t")
  
  
}






















sessionInfo()














################################################################
### 06_cytokines_bimatrix_top_selection_merged.R done!
################################################################