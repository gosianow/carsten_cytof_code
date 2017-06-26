

Sys.time()

# Load packages
library(FlowSOM)
library(ConsensusClusterPlus)
library(plyr) # for rbind.fill
library(limma)

##############################################################################
# Test arguments
##############################################################################


prefix='23CD4TmemCD69_02CD4_'
outdir='../carsten_cytof/PD1_project/CK_2016-06-23_02_CD4_merging2_Tmem_merging2_CD69/090_cytokine_bimatrix_frequencies_clustering'
path_data='../carsten_cytof/PD1_project/CK_2016-06-23_02_CD4_merging2_Tmem_merging2_CD69/090_cytokine_bimatrix/23CD4TmemCD69_02CD4_bimatrix.rds'
path_clustering_observables='../carsten_cytof/PD1_project/CK_2016-06-23_02_CD4_merging2_Tmem_merging2_CD69/090_cytokine_bimatrix/23CD4TmemCD69_02CD4_clustering_observables.xls'
rand_seed_consensus=1234
som_dim=7
nmetaclusts=49



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


rand_seed <- 1234
linkage <- "average"


if(!file.exists(outdir)) 
  dir.create(outdir, recursive = TRUE)


# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------


expr <- readRDS(path_data)


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

## Define the start codes as most frequent combinations
bivec <- apply(ef[, clust_observ, drop = FALSE], 1, paste0, collapse = "")

table_bivec <- table(bivec)

sort_freqs <- sort(table_bivec, decreasing = TRUE)
sort_comb <- names(sort_freqs)

start_codes <- strsplit2(sort_comb[1:(som_dim*som_dim)], "")
start_codes <- apply(start_codes, 2, as.numeric)

colnames(start_codes) <- clust_observ
start_codes <- start_codes[, scols]


# -------------------------------------
# SOM
# -------------------------------------

set.seed(rand_seed)
fsom <- FlowSOM::SOM(ef, xdim = som_dim, ydim = som_dim, rlen = 10, mst = 1, distf = 2, codes = start_codes)


if(som_dim^2 > nmetaclusts){
  
  # -------------------------------------
  # metaClustering_consensus 
  # -------------------------------------
  # consensus clustering that is reproducible with seed
  
  ### Sometimes not all the codes are used in mapping
  codes <- fsom$codes
  k <- nmetaclusts
  
  pdf(file.path(outdir, paste0(prefix, "ConsensusClusterPlus.pdf")), width = 7, height = 7)
  
  fccp <- ConsensusClusterPlus::ConsensusClusterPlus(t(codes),
    maxK = k, reps = 100, pItem = 0.9, pFeature = 1, title = tempdir(),
    plot = NULL, verbose = FALSE, clusterAlg = "hc", innerLinkage = linkage, finalLinkage = "average", distance = "euclidean", seed = rand_seed)
  
  dev.off()
  
  
  ### Get cluster ids
  fsom_mc <- fccp[[k]]$consensusClass
  
  clust <- fsom_mc[fsom$mapping[,1]]
  
  # Save clustering results
  saveRDS(fccp, file = file.path(outdir, paste0(prefix, "fccp.rds")))
  
}else{
  
  clust <- fsom$mapping[,1]
  
}



# Save clustering results

saveRDS(fsom, file = file.path(outdir, paste0(prefix, "fsom.rds")))

clust_out <- data.frame(cluster = clust, cell_id = cell_id, sample_id = samp, stringsAsFactors = FALSE)
write.table(clust_out, file = file.path(outdir, paste0(prefix, "clustering.xls")), row.names=FALSE, quote=FALSE, sep="\t")


# -------------------------------------
# Save clustering labels
# -------------------------------------

# get cluster frequencies
freq_clust <- table(clust)

# make data frame with labels
labels <- data.frame(cluster = sort(unique(clust)), label = sort(unique(clust)), counts = as.numeric(freq_clust))
labels$proportions <- round(labels$counts/sum(labels$counts) * 100, 2)


write.table(labels, file = file.path(outdir, paste0(prefix, "clustering_labels.xls")), row.names=FALSE, quote=FALSE, sep="\t")




























sessionInfo()



