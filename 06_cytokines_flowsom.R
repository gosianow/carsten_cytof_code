##############################################################################
## <<06_cytokines.R>>

# BioC 3.3
# Created 24 Aug 2016
# Updated 24 Aug 2016

##############################################################################
Sys.time()
##############################################################################

# Load packages
library(gdata)
library(Repitools)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2)
library(limma)
library(FlowSOM)
library(ConsensusClusterPlus)


##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_02_CD4_merging2'
cytokines_prefix='23CD4_02CD4_pca1_merging_Tmem_cytCM_raw_cl40_'
path_data=
nmetaclusts=40


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

prefix <- cytokines_prefix
suffix <- ""

prefix_clust <- paste0("cl", nmetaclusts, "_")

if(!data2analyse %in% c("raw", "norm"))
  stop("data2analyse must be 'raw' or 'norm'!")


# ------------------------------------------------------------
# define directories
# ------------------------------------------------------------

fcsDir <- "010_cleanfcs"; if( !file.exists(fcsDir) ) dir.create(fcsDir)
hmDir <- "030_heatmaps"; if( !file.exists(hmDir) ) dir.create(hmDir)
sneDir <- "040_tsnemaps"; if( !file.exists(sneDir) ) dir.create(sneDir)
cyDir <- "060_cytokines"; if( !file.exists(cyDir) ) dir.create(cyDir)

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

fcs_colnames <- colnames(fcs[[1]])

### Create sample info
samp <- rep(names(fcs), sapply(fcs, nrow))

# ------------------------------------------------------------
# Load more data
# ------------------------------------------------------------

## positive cutoffs for cytokines
cytokines_cutoffs <- read.xls(path_cytokines_cutoffs, stringsAsFactors=FALSE)

if(!cutoff_colname %in% colnames(cytokines_cutoffs))
  stop("There are no such column with cutoffs!")


# -------------------------------------
# get the isotope and antigen for fcs markers

m <- match(fcs_colnames, cytokines_cutoffs$fcs_colname)

fcs_panel <- data.frame(colnames = fcs_colnames, Isotope = cytokines_cutoffs$Isotope[m], Antigen = cytokines_cutoffs$Antigen[m], stringsAsFactors = FALSE)


# -------------------------------------

### Indeces of observables used for positive-negative analysis

pncols <- which(fcs_panel$colnames %in% cytokines_cutoffs[!is.na(cytokines_cutoffs[, cutoff_colname]), "fcs_colname"])

## Order them in the way as in panel
mm <- match(cytokines_cutoffs[!is.na(cytokines_cutoffs[, cutoff_colname]), "Antigen"], fcs_panel$Antigen[pncols])
pncols <- pncols[mm]


# ------------------------------------------------------------

## clustering results
clust <- read.table(file.path(hmDir, path_clustering), header = TRUE, sep = "\t", as.is = TRUE)
clust <- clust[, "cluster"]

## cluster labels
labels <- read.table(file.path(hmDir, path_clustering_labels), header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))


if(!all(clsubset %in% labels$label))
  stop("Cluster labels are wrong!")


# load marker selection for plotting on the heatmaps
marker_selection <- NULL

if(file.exists(file.path(path_marker_selection))){
  
  marker_selection <- read.table(file.path(path_marker_selection), header = TRUE, sep = "\t", as.is = TRUE)
  marker_selection <- marker_selection[, 1]
  
  if(!all(marker_selection %in% fcs_panel$Antigen[pncols]))
    stop("Marker selection is wrong")
  
}



# ------------------------------------------------------------
# Cell clustering with FlowSOM
# ------------------------------------------------------------

# Number of clusters
rand_seed <- 1234

# SOM
set.seed(rand_seed)
fsom <- FlowSOM::SOM(bm)


# consensus clustering that is reproducible with seed
data <- fsom$codes
k <- nmetaclusts

pdf(file.path(cyDir, paste0(prefix, prefix_clust, "ConsensusClusterPlus", suffix, ".pdf")), width = 7, height = 7)

results <- ConsensusClusterPlus::ConsensusClusterPlus(t(data),
  maxK = k, reps = 100, pItem = 0.9, pFeature = 1, title = tempdir(),
  plot = NULL, verbose = FALSE, clusterAlg = "hc", distance = "euclidean", seed = rand_seed)

dev.off()


# get cluster ids
fsom_mc <- results[[k]]$consensusClass
clust <- fsom_mc[fsom$mapping[,1]]


which_cells_kept <- seq(length(cells2keep))[cells2keep][cells2keep2]

### Save clustering results
clust_out <- data.frame(cluster = clust, cell_id = which_cells_kept, bm, stringsAsFactors = FALSE)

write.table(clust_out, file = file.path(cyDir, paste0(prefix, prefix_clust, "clustering", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")



### Save cluster frequencies and the mean expression
a <- aggregate(bm, by=list(clust), FUN = mean)

# get cluster frequencies
freq_clust <- table(clust)

clusters_out <- data.frame(cluster = names(freq_clust), label = names(freq_clust), counts = as.numeric(freq_clust), frequencies = as.numeric(freq_clust)/sum(freq_clust) * 100, a[, -1])

write.table(clusters_out, file.path(cyDir, paste0(prefix, prefix_clust, "clusters", suffix, ".xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)





























################################
### 06_cytokines done!
################################