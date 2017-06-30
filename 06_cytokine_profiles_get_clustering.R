

Sys.time()


##############################################################################
# Test arguments
##############################################################################


prefix='23CD4TmemCD69_29CD4TmemCD69_02CD4_cl49_'

outdir='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/02_CD4/090_cytokine_bimatrix_frequencies_clustering/cytokine_profiles'

path_bimatrix_clustering='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/02_CD4/090_cytokine_bimatrix_frequencies_clustering/23CD4TmemCD69_29CD4TmemCD69_02CD4_cl49_clustering.xls'

path_clustering=c('../carsten_cytof/PD1_project/CK_2016-06-23_02_CD4_merging2/030_heatmaps/23CD4_02CD4_pca1_merging2_clustering.xls','../carsten_cytof/PD1_project/CK_2016-06-29_02_CD4_merging/030_heatmaps/29CD4_02CD4_pca1_merging3_clustering.xls')

path_cells2keep=c('../carsten_cytof/PD1_project/CK_2016-06-23_02_CD4_merging2_Tmem_merging2_CD69/010_cleanfcs/cells2keep.txt','../carsten_cytof/PD1_project/CK_2016-06-29_02_CD4_merging_Tmem_merging3_CD69/010_cleanfcs/cells2keep.txt')

data=c('23CD4','29CD4')


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

if( !file.exists(outdir) ) 
  dir.create(outdir, recursive = TRUE)


# ------------------------------------------------------------
# Load bimatrix clustering
# ------------------------------------------------------------

bimatrix_clustering <- read.table(path_bimatrix_clustering, header = TRUE, sep = "\t", as.is = TRUE)


# ------------------------------------------------------------
# Load clustering from which the Tmem cells were extracted
# ------------------------------------------------------------


clustering_list <- lapply(1:length(path_clustering), function(i){
  clustering <- read.table(path_clustering[i], header = TRUE, sep = "\t", as.is = TRUE)
  clustering
})

clustering <- plyr::rbind.fill(clustering_list)
head(clustering)

lapply(clustering_list, dim)

# ------------------------------------------------------------
# Load info about which cells were extracted as Tmem
# ------------------------------------------------------------


cells2keep_list <- lapply(1:length(path_cells2keep), function(i){
  cells2keep <- read.table(path_cells2keep[i], header = TRUE, sep = "\t", as.is = TRUE)
  cells2keep$data <- paste0("data", data[i])
  cells2keep$cells2keep <- as.logical(cells2keep$cells2keep)
  cells2keep
})

lapply(cells2keep_list, dim)

cells2keep <- plyr::rbind.fill(cells2keep_list)
head(cells2keep)

cells2keep$data <- factor(cells2keep$data, levels = paste0("data", data))
table(cells2keep$data)

### Some checks
stopifnot(all(clustering$sample_id == cells2keep$sample_id))


# ------------------------------------------------------------
# Split the bimatrix clustering among the datasets that were merged
# ------------------------------------------------------------

if(length(clustering_list) > 1){
  
  bimatrix_clustering_list <- split(bimatrix_clustering, cells2keep$data[cells2keep$cells2keep])
  
}else{
  bimatrix_clustering_list <- list(bimatrix_clustering)
}

lapply(bimatrix_clustering_list, dim)

# ------------------------------------------------------------
# Save the clustering of Tmem cells for each of the bimatrix clusters
# ------------------------------------------------------------


for(i in 1:length(path_clustering)){
  # i = 1
  
  bimatrix_clusters <- unique(bimatrix_clustering_list[[i]]$cluster)
  table(bimatrix_clustering_list[[i]]$cluster)
  
  clustering_sub <- clustering_list[[i]][cells2keep_list[[i]]$cells2keep, , drop = FALSE]
  
  ### Some checks
  stopifnot(nrow(bimatrix_clustering_list[[i]]) == nrow(clustering_sub))
  
  for(j in 1:length(bimatrix_clusters)){
    # j = 1
    
    clustering_out <- clustering_sub[bimatrix_clustering_list[[i]]$cluster == bimatrix_clusters[j], ,drop = FALSE]
    
    write.table(clustering_out, file.path(outdir, paste0(prefix, "clustering_data", data[i], "_cl", bimatrix_clusters[j] ,".txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    
  }
  
}

































Sys.time()

sessionInfo()


