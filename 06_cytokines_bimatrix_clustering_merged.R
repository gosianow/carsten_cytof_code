##############################################################################
## <<06_cytokines_bimatrix_clustering_merged.R>>

# BioC 3.3
# Created 6 Nov 2016
# Updated 6 Nov 2016

# It clusters the 2 bimatrices using SOM
# It saves also the freq_out as counts.xls files which is necessary for the 08_frequencies_merged.R
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

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-merged_23_29/02_CD4'
clust_prefix='23_29_CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_cl16_'
clust_outdir='10_cytokines_merged/01_clustering'

path_data=c('/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_02_CD4_merging2/060_cytokines_bimatrix/01_clustering/23CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_bimatrix.txt','/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-29_02_CD4_merging/060_cytokines_bimatrix/01_clustering/29CD4_02CD4_pca1_merging3_Tmem_cytCM_raw2_bimatrix.txt')
path_clustering_observables=c('/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_02_CD4_merging2/060_cytokines_bimatrix/01_clustering/23CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_clustering_observables.xls','/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-29_02_CD4_merging/060_cytokines_bimatrix/01_clustering/29CD4_02CD4_pca1_merging3_Tmem_cytCM_raw2_clustering_observables.xls')

data_name=c('data23','data29')
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



# --------------------------------------------------------------------------
# run FlowSOM 
# --------------------------------------------------------------------------


# -------------------------------------
# SOM
# -------------------------------------

set.seed(rand_seed)
fsom <- FlowSOM::SOM(ef, xdim = som_dim, ydim = som_dim, rlen = 10, mst = 1)


data <- fsom$codes

clust <- fsom$mapping[,1]


### Save clustering results

clust_out <- data.frame(cluster = clust, cell_id = cell_id, sample_id = samp, stringsAsFactors = FALSE)

clust_out <- split(clust_out, factor(data_id, levels = data_name))

for(i in 1:length(clust_out))
  write.table(clust_out[[i]], file = file.path(outdir, paste0(prefix, data_name[i], "_", "clustering.xls")), row.names=FALSE, quote=FALSE, sep="\t")



# make data frame with labels
labels <- data.frame(cluster = sort(unique(clust)), label = sort(unique(clust)))

for(i in 1:length(data_name))
  write.table(labels, file = file.path(outdir, paste0(prefix, data_name[i], "_", "clustering_labels.xls")), row.names=FALSE, quote=FALSE, sep="\t")



### Code sizes

code_sizes <- table(fsom$mapping[, 1]) 

# Sometimes not all the codes have mapped cells so they will have size 0
code_sizes_full <- rep(0, nrow(data))
names(code_sizes_full) <- 1:nrow(data)

code_sizes_full[names(code_sizes)] <- as.numeric(code_sizes)


### Plot codes as a heatmap
rownames(data) <- 1:nrow(data)

color <- rev(colorRampPalette(brewer.pal(n = 8, name = "Spectral"))(100))
labels_col <- clustering_observables[clustering_observables$clustering_observable, "marker"]
annotation_colors <- NA

cluster_rows <- hclust(dist(data), method = "ward.D2")

labels_row <- paste0(rownames(data), "  ( ", format(code_sizes_full, big.mark=",", scientific=FALSE), " )")
annotation_row <- NA


pheatmap(data[, clustering_observables[clustering_observables$clustering_observable, "mass"]], color = color, cellwidth = 24, cellheight = 12, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, display_numbers = TRUE, number_color = "black", fontsize_number = 6, fontsize_row = 7, fontsize_col = 14, fontsize = 12, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "codes_pheatmap.pdf")))


### Save codes


codes <- data.frame(code_id = rownames(data), cluster = 1:nrow(data), size = code_sizes_full, data)

write.table(codes, file = file.path(outdir, paste0(prefix, "codes.xls")), row.names=FALSE, quote=FALSE, sep="\t")






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
  mlab <- match(rownames(freq), labels$cluster)
  
  
  ### Save the frequencies and proportions
  prop_out <- data.frame(labels[mlab, c("cluster", "label")], as.data.frame.matrix(prop))
  
  freq_out <- data.frame(labels[mlab, c("cluster", "label")], as.data.frame.matrix(freq))
  
  write.table(prop_out, file=file.path(outdir, paste0(prefix, data_name[i], "_", "frequencies.xls")), row.names=FALSE, quote=FALSE, sep="\t")
  write.table(freq_out, file=file.path(outdir, paste0(prefix, data_name[i], "_", "counts.xls")), row.names=FALSE, quote=FALSE, sep="\t")
  
  
}






















sessionInfo()














################################################################
### 06_cytokines_bimatrix_clustering_merged.R done!
################################################################