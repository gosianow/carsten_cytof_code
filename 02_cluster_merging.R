##############################################################################
## <<02_cluster_merging.R>>

# BioC 3.3
# Created 28 July 2016
# Updated 

##############################################################################
Sys.time()
##############################################################################

library(gdata)

##############################################################################
# Test arguments
##############################################################################

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
# merging_prefix='pca1_mergingNEW_'
# path_cluster_merging='cluster_mergingNEW.xlsx'
# path_clustering='pca1_cl20_clustering.xls'

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


# ------------------------------------------------------------
# define directories
# ------------------------------------------------------------

fcsDir <- "010_cleanfcs"; if( !file.exists(fcsDir) ) dir.create(fcsDir)
pcaDir <- "020_pcascores"; if( !file.exists(pcaDir) ) dir.create(pcaDir)
hmDir <- "030_heatmaps"; if( !file.exists(hmDir) ) dir.create(hmDir)


# ------------------------------------------------------------
# Load clustering data
# ------------------------------------------------------------

# read cluster merging file
cm <- read.xls(path_cluster_merging)

cm

# remove spaces in labels bcs they are problematic...
cm$label <- factor(cm$label, labels = gsub(" ", "_", levels(cm$label))) 


# read original clustering
clust <- read.table(file.path(hmDir, path_clustering), header = TRUE, sep = "\t", as.is = TRUE)
clust <- clust[, 1]



# ------------------------------------------------------------
# Get new merged clustering
# ------------------------------------------------------------

# new clustering
clustm <- factor(clust, levels = cm$old_cluster)
levels(clustm) <- cm$new_cluster
clustm <- as.numeric(as.character(clustm))


# get cluster frequencies
freq_clust <- table(clustm)


# new labels
labels <- unique(cm[,c("new_cluster","label")])
colnames(labels) <- c("cluster", "label")

labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))
labels


if(sum(duplicated(labels$label)) > 0 | sum(duplicated(labels$cluster)) > 0)
  stop("Wrong merging file")



### Save cluster merging results


freq_clust_out <- data.frame(cluster = names(freq_clust), freq = as.numeric(freq_clust))
write.table(freq_clust_out, file=file.path(hmDir, paste0(merging_prefix, "cluster_counts.xls")), row.names=FALSE, quote=FALSE, sep="\t")

clust_out <- data.frame(cluster = clustm, stringsAsFactors = FALSE)
write.table(clust_out, file = file.path(hmDir, paste0(merging_prefix, "clustering.xls")), row.names=FALSE, quote=FALSE, sep="\t")


write.table(labels, file = file.path(hmDir, paste0(merging_prefix, "clustering_labels.xls")), row.names=FALSE, quote=FALSE, sep="\t")









sessionInfo()













################################
### 02_cluster_merging done!
################################
