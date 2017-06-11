
Sys.time()


library(flowCore)
library(gdata)


##############################################################################
# Test arguments
##############################################################################

outdir='../carsten_cytof/PD1_project/CK_2016-06-23_01_CD4_mergingNEW2/010_cleanfcs'
dir_fcs='../carsten_cytof/PD1_project/CK_2016-06-23_01/010_cleanfcs'
path_metadata='../carsten_cytof/PD1_project/CK_metadata/metadata_23_01.xlsx'
path_clustering='../carsten_cytof/PD1_project/CK_2016-06-23_01/030_heatmaps/23_01_pca1_mergingNEW2_clustering.xls'
path_clustering_labels='../carsten_cytof/PD1_project/CK_2016-06-23_01/030_heatmaps/23_01_pca1_mergingNEW2_clustering_labels.xls'
extract_cluster='CD4'


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


# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

# read in metadata
md <- gdata::read.xls(path_metadata, stringsAsFactors=FALSE)


# define FCS file names
f <- file.path(dir_fcs, md$filename)
names(f) <- md$shortname


# read in raw FCS files 
fcs <- lapply(f, read.FCS)



# ------------------------------------------------------------
# Load clustering data
# ------------------------------------------------------------

clustering <- read.table(path_clustering, header = TRUE, sep = "\t", as.is = TRUE)
clust <- clustering[, "cluster"]

labels <- read.table(path_clustering_labels, header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))


# ------------------------------------------------------------
# Save fcs files for a cluster to export
# ------------------------------------------------------------

## use labels as cluster names
mlab <- match(clust, labels$cluster)
clust <- labels$label[mlab]
  
clustList <- split(clust, factor(clustering[, "sample_id"], levels = names(fcs)))


writeOutCluster <- function(u, v, z, keep, outdir) {
  # u - cluster; v - flowFrame; z - original filename
  
  fn <- file.path(outdir, basename(z))
  
  # out <- v[u == keep, ]
  # out <- v[grep(keep, u), ]
  out <- v[u %in% keep, ]
  
  write.FCS(out, fn)
  
  return(NULL)
  
}



dummy <- mapply(writeOutCluster, u = clustList, v = fcs, z = f, keep = extract_cluster, outdir = outdir)







sessionInfo()




