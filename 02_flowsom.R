##############################################################################
## <<02_flowsom.R>>

# BioC 3.3
# Created 27 July 2016

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
library(ConsensusClusterPlus)


##############################################################################
# Test arguments
##############################################################################

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
# path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.xlsx'
# flowsom_prefix='23_01_pca1_cl20_'
# path_clustering_observables='23_01_pca1_clustering_observables.xls'
# rand_seed_consensus=123
# nmetaclusts=20


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
rand_seed <- 1234

prefix <- flowsom_prefix


# ### Colors for 20 clusters 
# # ggplot palette
# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length=n+1)
#   hcl(h=hues, l=60 , c=100)[1:n]
# }
# 
# color_ramp <- c(colorRampPalette(brewer.pal(12,"Paired"))(12)[-c(11)],  gg_color_hue(max(1, nmetaclusts-11)) )


# ------------------------------------------------------------
# define directories
# ------------------------------------------------------------

fcsDir <- "010_cleanfcs"; if( !file.exists(fcsDir) ) dir.create(fcsDir)
pcaDir <- "020_pcascores"; if( !file.exists(pcaDir) ) dir.create(pcaDir)
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


# ------------------------------------------------------------
# Load more data
# ------------------------------------------------------------

if(!grepl("/", path_clustering_observables)){
  clust_observ <- read.table(file.path(hmDir, path_clustering_observables), header = TRUE, sep = "\t", as.is = TRUE)
}else{
  clust_observ <- read.table(path_clustering_observables, header = TRUE, sep = "\t", as.is = TRUE)
}

clust_observ <- clust_observ[, 1]

# -------------------------------------

# selected columns for clustering 

scols <- which(colnames(fcs[[1]]) %in% clust_observ)

# arc-sin-h the columns specific 
fcsT <- lapply(fcs, function(u) {
  e <- exprs(u)
  e[ , scols] <- asinh( e[ , scols] / 5 )
  exprs(u) <- e
  u
})

# -------------------------------------
# run FlowSOM and make heatmaps
# -------------------------------------



# Number of clusters


fs <- as(fcsT,"flowSet")

fsom <- FlowSOM::ReadInput(fs, transform = FALSE, scale = FALSE)

# SOM
set.seed(rand_seed)
fsom <- FlowSOM::BuildSOM(fsom, colsToUse = scols)


# consensus clustering
# fsom_mc <- FlowSOM::metaClustering_consensus(fsom$map$codes, k=nmetaclusts)

# consensus clustering that is reproducible with seed
data <- fsom$map$codes
k <- nmetaclusts

pdf(file.path(hmDir, paste0(prefix, "ConsensusClusterPlus.pdf")), width = 7, height = 7)

results <- ConsensusClusterPlus::ConsensusClusterPlus(t(data),
  maxK = k, reps = 100, pItem = 0.9, pFeature = 1, title = tempdir(),
  plot = NULL, verbose = FALSE, clusterAlg = "hc", distance = "euclidean", seed = rand_seed_consensus)

dev.off()



# ## For plotting tracking plot with my colors
# pdf(file.path(hmDir, paste0(prefix, "ConsensusClusterPlus_bar.pdf")), width = 14, height = 4)
# barplot(rep(1, length(results[[k]]$clrs[[1]])), col = results[[k]]$clrs[[1]][results[[k]]$consensusTree$order])
# dev.off()
# 
# results[[k]]$consensusClass[results[[k]]$consensusTree$order]



fsom_mc <- results[[k]]$consensusClass

# get cluster ids
clust <- fsom_mc[fsom$map$mapping[,1]]

# get cluster frequencies
freq_clust <- table(clust)




### Save clustering results

save(fsom, file = file.path(hmDir, paste0(prefix, "fsom.rda")))
save(fsom_mc, file = file.path(hmDir, paste0(prefix, "fsom_mc.rda")))


clust_out <- data.frame(cluster = clust, stringsAsFactors = FALSE)
write.table(clust_out, file = file.path(hmDir, paste0(prefix, "clustering.xls")), row.names=FALSE, quote=FALSE, sep="\t")

# make data frame with labels
labels <- data.frame(cluster = 1:nmetaclusts, label = sprintf("%02d", 1:nmetaclusts))
write.table(labels, file = file.path(hmDir, paste0(prefix, "clustering_labels.xls")), row.names=FALSE, quote=FALSE, sep="\t")


sessionInfo()













################################
### 02_flowsom done!
################################