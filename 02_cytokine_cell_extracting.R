
Sys.time()


library(flowCore)
library(gdata)
library(limma)

##############################################################################
# Test arguments
##############################################################################

outdir='../carsten_cytof/PD1_project/CK_2016-06-23_02_CD4_merging2_Tmem_merging2_CD69/010_cleanfcs'
dir_fcs='../carsten_cytof/PD1_project/CK_2016-06-23_02_CD4_merging2/010_cleanfcs'
path_metadata='../carsten_cytof/PD1_project/CK_metadata/metadata_23_02.xlsx'
path_clustering='../carsten_cytof/PD1_project/CK_2016-06-23_02_CD4_merging2/030_heatmaps/23CD4_02CD4_pca1_merging2_clustering.xls'
path_clustering_labels='../carsten_cytof/PD1_project/CK_2016-06-23_02_CD4_merging2/030_heatmaps/23CD4_02CD4_pca1_merging2_clustering_labels.xls'
path_extract_cluster='../carsten_cytof/PD1_project/CK_2016-06-23_02_CD4_merging2/010_helpfiles/23CD4_02CD4_pca1_merging2_extract_cluster_Tmem.txt'
path_extract_marker='../carsten_cytof/PD1_project/CK_panels/panel2CD4_23_cytokines_CD69.xlsx'


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
# Load metadata
# ------------------------------------------------------------

md <- gdata::read.xls(path_metadata, stringsAsFactors=FALSE)

# ------------------------------------------------------------
# Load clustering data
# ------------------------------------------------------------

clustering <- read.table(path_clustering, header = TRUE, sep = "\t", as.is = TRUE)
clust <- clustering[, "cluster"]

labels <- read.table(path_clustering_labels, header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))
labels

# ------------------------------------------------------------
# Load clusters to extract
# ------------------------------------------------------------

extract_cluster <- read.table(path_extract_cluster, stringsAsFactors = FALSE)[, 1]
extract_cluster

stopifnot(all(extract_cluster %in% labels$label))

# ------------------------------------------------------------
# Load marker expression to extract
# ------------------------------------------------------------



extract_marker <- gdata::read.xls(path_extract_marker, stringsAsFactors = FALSE)
extract_marker <- extract_marker[complete.cases(extract_marker), , drop = FALSE]
extract_marker

if(nrow(extract_marker) > 0){
  extract_marker <- extract_marker[, c("fcs_colname", "positive_cutoff_raw_base", "positive_cutoff_raw_tx"), drop = FALSE]
  colnames(extract_marker) <- c("fcs_colname", "base", "tx")
  extract_marker
}





# ------------------------------------------------------------
# Load FCS files
# ------------------------------------------------------------

f <- file.path(dir_fcs, md$filename)
names(f) <- md$shortname
names(f)

fcs <- lapply(f, read.FCS)
sapply(fcs, nrow)

# ------------------------------------------------------------
# Save fcs files for a cluster to export
# ------------------------------------------------------------

## Use labels as cluster names
mm <- match(clust, labels$cluster)
clust <- labels$label[mm]

clust_split <- split(clust, factor(clustering[, "sample_id"], levels = names(fcs)))
sapply(clust_split, length)

## Some checks
stopifnot(all(sapply(fcs, nrow) == sapply(clust_split, length)))

if(!all(names(fcs) == names(clust_split) & names(fcs) == md$shortname))
  stop("Wrong order of samples!")


cells2keep <- lapply(names(fcs), function(i){
  # i = "base_HD1"
  
  ## Filtering based on marker expression
  if(nrow(extract_marker) > 0){
    cols <- extract_marker$fcs_colname
    e <- asinh(exprs(fcs[[i]])[, cols, drop = FALSE] / 5 )
    
    day <- strsplit2(i, "_")[,1]
    
    cells2keep1 <- colSums(t(e) > extract_marker[, day]) == nrow(extract_marker)
    table(cells2keep1)
  }else{
    cells2keep1 <- rep(TRUE, nrow(fcs[[i]]))
  }
  
  ## Filtering based on clusters
  cells2keep2 <- clust_split[[i]] %in% extract_cluster
  table(cells2keep2)
  
  cells2keep <- cells2keep1 & cells2keep2
  
  out <- fcs[[i]][cells2keep, ]
  
  fn <- file.path(outdir, basename(f[i]))
  write.FCS(out, fn)
  
  return(cells2keep)
  
})


cells2keep_out <- data.frame(cell_id = clustering$cell_id, sample_id = clustering$sample_id, cells2keep = as.numeric(unlist(cells2keep))) 

write.table(cells2keep_out, file.path(outdir, paste0("cells2keep.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)









sessionInfo()




