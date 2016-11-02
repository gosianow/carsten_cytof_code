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
# merging_outdir='030_heatmaps'
# path_cluster_merging='cluster_mergingNEW.xlsx'
# path_clustering='030_heatmaps/pca1_cl20_clustering.xls'

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

prefix <- merging_prefix
outdir <- merging_outdir

if( !file.exists(outdir) ) 
  dir.create(outdir)


# ------------------------------------------------------------
# Load clustering data
# ------------------------------------------------------------

# read cluster merging file
cm <- read.xls(path_cluster_merging)

cm

if(!all(c("old_cluster", "label", "new_cluster") %in% colnames(cm)))
  stop("Merging file must contain 'old_cluster', 'label' and 'new_cluster' columns!")


# remove spaces in labels bcs they are problematic...
cm$label <- factor(cm$label, labels = gsub(" ", "_", levels(cm$label))) 


# read original clustering
clustering <- read.table(path_clustering, header = TRUE, sep = "\t", as.is = TRUE)
clust <- clustering[, "cluster"]



# ------------------------------------------------------------
# Get new merged clustering
# ------------------------------------------------------------

# new labels
labels <- unique(cm[,c("new_cluster","label")])
colnames(labels) <- c("cluster", "label")

labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))

if(sum(duplicated(labels$label)) > 0 | sum(duplicated(labels$cluster)) > 0)
  stop("Wrong merging file")


# new clustering
clustm <- factor(clust, levels = cm$old_cluster)
levels(clustm) <- cm$new_cluster
clustm <- as.numeric(as.character(clustm))


# get cluster frequencies
freq_clust <- table(clustm)
labels$counts <- as.numeric(freq_clust)
labels$proportions <- round(labels$counts/sum(labels$counts) * 100, 2)

labels




### Save cluster merging results

clust_out <- data.frame(cluster = clustm, cell_id = clustering$cell_id, sample_id = clustering$sample_id, stringsAsFactors = FALSE)
write.table(clust_out, file = file.path(outdir, paste0(merging_prefix, "clustering.xls")), row.names=FALSE, quote=FALSE, sep="\t")


write.table(labels, file = file.path(outdir, paste0(merging_prefix, "clustering_labels.xls")), row.names=FALSE, quote=FALSE, sep="\t")









sessionInfo()













################################
### 02_cluster_merging done!
################################
