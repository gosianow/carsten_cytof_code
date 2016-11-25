##############################################################################
## <<06_cytokines_bimatrix_clustering_merged.R>>

# BioC 3.3
# Created 6 Nov 2016
# Updated 12 Nov 2016

# It clusters the 2 bimatrices using SOM or SOM + cluster consensus
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
library(UpSetR)
library(limma)


##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-merged_23_29/02_CD4'
clust_prefix='23_29_CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_cl40_'
clust_outdir='10_cytokines_merged/01_clustering'

path_data=c('/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_02_CD4_merging2/060_cytokines_bimatrix/01_clustering/23CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_bimatrix.txt','/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-29_02_CD4_merging/060_cytokines_bimatrix/01_clustering/29CD4_02CD4_pca1_merging3_Tmem_cytCM_raw2_bimatrix.txt')
path_clustering_observables=c('/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_02_CD4_merging2/060_cytokines_bimatrix/01_clustering/23CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_clustering_observables.xls','/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-29_02_CD4_merging/060_cytokines_bimatrix/01_clustering/29CD4_02CD4_pca1_merging3_Tmem_cytCM_raw2_clustering_observables.xls')

data_name=c('data23','data29')
som_dim=10
nmetaclusts=40

linkage='average'
distance='euclidean'


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



if(all(!grepl("linkage=", args))){
  linkage <- "average"
}

if(all(!grepl("distance=", args))){
  distance <- "euclidean"
}


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
color_ramp <- c(colors_muted, gg_color_hue(max(1, nmetaclusts - length(colors_muted))))

colors_clusters <- color_ramp[1:nmetaclusts]
names(colors_clusters) <- 1:nmetaclusts


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
# run FlowSOM 
# --------------------------------------------------------------------------


# -------------------------------------
# SOM
# -------------------------------------

## Define the start codes as most frequent combinations
bivec <- apply(bidf, 1, paste0, collapse = "")

table_bivec <- table(bivec)

sort_freqs <- sort(table_bivec, decreasing = TRUE)
sort_comb <- names(sort_freqs)

start_codes <- strsplit2(sort_comb[1:(som_dim*som_dim)], "")
start_codes <- apply(start_codes, 2, as.numeric)

colnames(start_codes) <- clustering_observables[clustering_observables$clustering_observable, "mass"]
start_codes <- start_codes[, scols]

## Run SOM
set.seed(rand_seed)
fsom <- FlowSOM::SOM(ef, xdim = som_dim, ydim = som_dim, rlen = 10, mst = 1, distf = 2, codes = start_codes)


data <- fsom$codes
rownames(data) <- 1:nrow(data)

data_clust <- data


### Keep non zero codes and convert them to 01 to be able to use the binary distance
# data <- data[rownames(data) %in% names(table(fsom$mapping[,1])), ]
# data_01 <- data
# data_01[data_01 < 0.5] <- 0
# data_01[data_01 > 0.5] <- 1
# data_clust <- data_01



if(som_dim^2 > nmetaclusts){
  
  # -------------------------------------
  # metaClustering_consensus 
  # -------------------------------------
  # consensus clustering that is reproducible with seed
  
  ### Sometimes not all the codes are used in mapping
  k <- nmetaclusts
  
  dist_definition <- function(x)
    dist(x, method = distance)
  
  
  pdf(file.path(outdir, paste0(prefix, "ConsensusClusterPlus.pdf")), width = 7, height = 7)
  
  results <- ConsensusClusterPlus::ConsensusClusterPlus(t(data_clust),
    maxK = k, reps = 100, pItem = 0.9, pFeature = 1, title = tempdir(),
    plot = NULL, verbose = FALSE, clusterAlg = "hc", innerLinkage = linkage, finalLinkage = "average", distance = "dist_definition", seed = rand_seed)
  
  dev.off()
  
  # get cluster ids
  fsom_mc <- results[[k]]$consensusClass
  fsom_mc_tree <- results[[k]]$consensusTree
  
  clust <- fsom_mc[as.character(fsom$mapping[,1])]
  
  # elements for the heatmap
  annotation_row <- data.frame(cluster = factor(fsom_mc))
  rownames(annotation_row) <- rownames(data)
  annotation_colors <- list(cluster = colors_clusters)
  
  ### Save clustering results
  save(fsom, file = file.path(outdir, paste0(prefix, "fsom.rda")))
  save(fsom_mc, file = file.path(outdir, paste0(prefix, "fsom_mc.rda")))
  
}else{
  
  clust <- fsom$mapping[,1]
  
  fsom_mc <- 1:nrow(data)
  fsom_mc_tree <- hclust(dist(data_clust, method = distance), method = linkage)
  
  # elements for the heatmap
  annotation_row <- NA
  annotation_colors = NA
  
}




### Save clustering results

clust_out <- data.frame(cluster = clust, cell_id = cell_id, sample_id = samp, stringsAsFactors = FALSE)

clust_out <- split(clust_out, factor(data_id, levels = data_name))

for(i in 1:length(clust_out))
  write.table(clust_out[[i]], file = file.path(outdir, paste0(prefix, data_name[i], "_", "clustering.xls")), row.names=FALSE, quote=FALSE, sep="\t")



### Code sizes

code_sizes <- table(fsom$mapping[, 1]) 

# Sometimes not all the codes have mapped cells so they will have size 0
code_sizes_full <- rep(0, nrow(data))
names(code_sizes_full) <- rownames(data)

code_sizes_full[names(code_sizes)] <- as.numeric(code_sizes)


### Plot codes as a heatmap


color <- rev(colorRampPalette(brewer.pal(n = 8, name = "Spectral"))(100))
labels_col <- clustering_observables[clustering_observables$clustering_observable, "marker"]

cluster_rows <- fsom_mc_tree

labels_row <- paste0(rownames(data), "  ( ", format(code_sizes_full, big.mark=",", scientific=FALSE), " )")


pheatmap(data[, clustering_observables[clustering_observables$clustering_observable, "mass"]], color = color, cellwidth = 24, cellheight = 12, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, display_numbers = TRUE, number_color = "black", fontsize_number = 6, fontsize_row = 7, fontsize_col = 14, fontsize = 10, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "codes_pheatmap.pdf")))


### Save codes


codes <- data.frame(code_id = rownames(data), cluster = fsom_mc, size = code_sizes_full, data)

write.table(codes, file = file.path(outdir, paste0(prefix, "codes.xls")), row.names=FALSE, quote=FALSE, sep="\t")



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


labels <- list()

for(i in 1:length(clust_out)){
  # i = 2
  
  clust_tmp <- clust_out[[i]]$cluster
  
  cluster <- as.character(sort(unique(clust_tmp)))
  
  label <- paste0(cluster, " - ", cluster_binames[cluster], " (", prop_clust[cluster], ")")
  
  labels[[i]] <- data.frame(cluster = cluster, label = label, stringsAsFactors = FALSE)
  
  write.table(labels[[i]], file = file.path(outdir, paste0(prefix, data_name[i], "_", "clustering_labels.xls")), row.names=FALSE, quote=FALSE, sep="\t")
  
}




### Plot heatmap of generated clusters


color <- rev(colorRampPalette(brewer.pal(n = 8, name = "Spectral"))(100))
labels_col <- colnames(data_new_labels)

cluster_rows <- hclust(dist(data_new_labels), method = linkage)

cluster <- as.character(sort(unique(clust)))

labels_row <- paste0(cluster, " - ", cluster_binames[cluster], " (", prop_clust[cluster], ")")

annotation_row <- data.frame(cluster = factor(as.numeric(cluster)))
rownames(annotation_row) <- rownames(data_new_labels)
annotation_colors <- list(cluster = colors_clusters)

pheatmap(data_new_labels, color = color, cellwidth = 20, cellheight = 20, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, display_numbers = TRUE, number_color = "black", fontsize_number = 6, fontsize_row = 10, fontsize_col = 10, fontsize = 10, annotation_row = annotation_row, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "clusters_pheatmap.pdf")))



### Save cluster frequencies and the median expression

clusters_out <- data.frame(cluster = cluster, label = labels_row, counts = as.numeric(freq_clust), frequencies = as.numeric(freq_clust)/sum(freq_clust), data_new_labels)

write.table(clusters_out, file.path(outdir, paste0(prefix, "clusters.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



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
### 06_cytokines_bimatrix_clustering_merged.R done!
################################################################