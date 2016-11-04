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
clust_prefix='23CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_cl16_'
clust_outdir='060_cytokines_bimatrix/01_clustering'
path_data='060_cytokines_bimatrix/01_clustering/23CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_bimatrix.txt'
path_clustering_observables='060_cytokines_bimatrix/01_clustering/23CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_clustering_observables.xls'
som_dim=4

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

clust <- fsom$mapping[,1]


### Save clustering results

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

data_order <- match(colnames(data), clustering_observables$mass)

color <- rev(colorRampPalette(brewer.pal(n = 8, name = "Spectral"))(100))
labels_col <- clustering_observables$marker[1:length(data_order)]
annotation_colors <- NA

cluster_rows <- hclust(dist(data), method = "ward.D2")

labels_row <- paste0(rownames(data), "  ( ", format(code_sizes_full, big.mark=",", scientific=FALSE), " )")
annotation_row <- NA


pheatmap(data[, data_order], color = color, cellwidth = 24, cellheight = 12, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, display_numbers = TRUE, number_color = "black", fontsize_number = 6, fontsize_row = 7, fontsize_col = 14, fontsize = 12, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "codes_pheatmap.pdf")))


### Save codes


codes <- data.frame(code_id = rownames(data), cluster = 1:nrow(data), size = code_sizes_full, data)

write.table(codes, file = file.path(outdir, paste0(prefix, "codes.xls")), row.names=FALSE, quote=FALSE, sep="\t")




























sessionInfo()














################################
### 06_cytokines_bimatrix_clustering done!
################################