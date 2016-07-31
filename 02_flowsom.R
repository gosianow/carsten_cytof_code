##############################################################################
## <<02_flowsom.R>>

# BioC 3.3
# Created 27 July 2016
# Updated 

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

##############################################################################
# Test arguments
##############################################################################

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
# flowsom_prefix=''
# path_clustering_observables='pca1_clustering_observables.xls'

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


### Set seed to be reproducible with FlowSOM clustering
options(myseed = rand_seed)
set.seed(rand_seed)

getOption("myseed")

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
md <- read.xls("metadata.xlsx",stringsAsFactors=FALSE)

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
nmetaclusts <- 20

set.seed(rand_seed)

fs <- as(fcsT,"flowSet")

system.time( fsom <- FlowSOM::ReadInput(fs, transform = FALSE, scale = FALSE) )
system.time( fsom <- FlowSOM::BuildSOM(fsom, colsToUse = scols) )
system.time( fsom_mc <- FlowSOM::metaClustering_consensus(fsom$map$codes, k=nmetaclusts) )

# get cluster ids
clust <- fsom_mc[fsom$map$mapping[,1]]

# get cluster frequencies
freq_clust <- table(clust)




### Save clustering results

save(fsom, file = file.path(hmDir, paste0(flowsom_prefix, "fsom.rda")))
save(fsom_mc, file = file.path(hmDir, paste0(flowsom_prefix, "fsom_mc.rda")))


freq_clust_out <- data.frame(cluster = names(freq_clust), freq = as.numeric(freq_clust))
write.table(freq_clust_out, file=file.path(hmDir, paste0(flowsom_prefix, "cluster_counts.xls")), row.names=FALSE, quote=FALSE, sep="\t")


clust_out <- data.frame(cluster = clust, stringsAsFactors = FALSE)
write.table(clust_out, file = file.path(hmDir, paste0(flowsom_prefix, "clustering.xls")), row.names=FALSE, quote=FALSE, sep="\t")

# make data frame with labels
labels <- data.frame(cluster = 1:nmetaclusts, label = sprintf("%02d", 1:nmetaclusts))
write.table(labels, file = file.path(hmDir, paste0(flowsom_prefix, "clustering_labels.xls")), row.names=FALSE, quote=FALSE, sep="\t")


sessionInfo()













################################
### 02_flowsom done!
################################