##############################################################################
## <<06_cytokines_expression.R>>

# BioC 3.3
# Created 27 Sep 2016
# Updated 27 Sep 2016

# Prepare the expression of cytokines from the Tmem clusters

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
library(RColorBrewer)
library(UpSetR)
library(limma)


##############################################################################
# Test arguments
##############################################################################

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_02_CD4_merging2'
# cytokines_prefix='23CD4_02CD4_pca1_merging_Tmem_cytCM_'
# cytokines_outdir='060_cytokines_expression'
# path_data='010_data/23CD4_02CD4_expr_raw.rds'
# path_cytokines_cutoffs='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_panels/panel2CD4_cytokines_CM.xlsx'
# path_clustering='030_heatmaps/23CD4_02CD4_pca1_merging2_clustering.xls'
# path_clustering_labels='030_heatmaps/23CD4_02CD4_pca1_merging2_clustering_labels.xls'
# clsubset=c('CM','EM','TM','TE')
# cutoff_colname='positive_cutoff_raw_base'

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
outdir <- cytokines_outdir

if(!file.exists(outdir)) 
  dir.create(outdir)

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

# read expression data
if(grepl(".txt", path_data)){
  expr <- read.table(path_data, header = TRUE, sep = "\t", as.is = TRUE)
}
if(grepl(".rds", path_data)){
  expr <- readRDS(path_data)
}

cell_id <- expr[, "cell_id"]
samp <- expr[, "sample_id"]
fcs_colnames <- colnames(expr)[!grepl("cell_id|sample_id", colnames(expr))]
e <- expr[, fcs_colnames]


# ------------------------------------------------------------
# Load clustering results
# ------------------------------------------------------------


## clustering results
clustering <- read.table(path_clustering, header = TRUE, sep = "\t", as.is = TRUE)
clust <- clustering[, "cluster"]

## cluster labels
labels <- read.table(path_clustering_labels, header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))


if(!all(clsubset %in% labels$label))
  stop("Cluster labels are wrong!")


# ------------------------------------------------------------
# Load cytokines_cutoffs
# ------------------------------------------------------------

## positive cutoffs for cytokines
cytokines_cutoffs <- read.xls(path_cytokines_cutoffs, stringsAsFactors=FALSE)

if(!all(cutoff_colname %in% colnames(cytokines_cutoffs)))
  stop("There are no such column with cutoffs!")


# ------------------------------------------------------------
# Do the subsetting of cells and markers used for the analysis
# ------------------------------------------------------------


# -------------------------------------
# Get the isotope and antigen for fcs markers

m <- match(fcs_colnames, cytokines_cutoffs$fcs_colname)

fcs_panel <- data.frame(fcs_colname = fcs_colnames, Isotope = cytokines_cutoffs$Isotope[m], Antigen = cytokines_cutoffs$Antigen[m], stringsAsFactors = FALSE)


# -------------------------------------
# Indeces of observables used for positive-negative analysis

pn_markers <- complete.cases(cytokines_cutoffs[, cutoff_colname, drop = FALSE]) 

pncols <- which(fcs_colnames %in% cytokines_cutoffs[pn_markers, "fcs_colname"])

## Order them in the way as in panel
mm <- match(cytokines_cutoffs[pn_markers, "Antigen"], fcs_panel$Antigen[pncols])
pncols <- pncols[mm]


# ------------------------------------------------------------
# Keep only the data from the specified clusters

cells2keep_clust <- clust %in% labels[labels$label %in% clsubset, "cluster"]


# ------------------------------------------------------------
# Save the expression of cytokines for Tmem cells
# ------------------------------------------------------------

exprb <- expr[cells2keep_clust, c("cell_id", "sample_id", fcs_colnames[pncols])]

saveRDS(exprb, file.path(outdir, paste0(prefix, "expr_raw.rds")))


# ------------------------------------------------------------
# Create a table with clustering observables, clustering and clustering_labels - needed for 04_expression.R to work 
# ------------------------------------------------------------


clustering_observables <- data.frame(mass = fcs_panel$fcs_colname[pncols], marker = fcs_panel$Antigen[pncols], clustering_observable = TRUE, stringsAsFactors = FALSE)

write.table(clustering_observables, file.path(outdir, paste0(prefix, "clustering_observables.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



clust_out <- data.frame(cluster = 1, cell_id = exprb$cell_id, sample_id = exprb$sample_id, stringsAsFactors = FALSE)
write.table(clust_out, file = file.path(outdir, paste0(prefix, "clustering.xls")), row.names=FALSE, quote=FALSE, sep="\t")


# get cluster frequencies
clust <- clust_out$cluster
freq_clust <- table(clust)

# make data frame with labels
labels <- data.frame(cluster = 1, label = "Tmem", counts = as.numeric(freq_clust))
write.table(labels, file = file.path(outdir, paste0(prefix, "clustering_labels.xls")), row.names=FALSE, quote=FALSE, sep="\t")


































################################
### 06_cytokines_expression.R done!
################################