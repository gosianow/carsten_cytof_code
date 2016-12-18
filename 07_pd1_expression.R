##############################################################################
## <<07_pd1_expression.R>>

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

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_02_CD4_merging2'
pd1_prefix='23CD4_02CD4_pca1_merging_Tmem_cytCM_'
pd1_outdir='070_pd1_expression'
path_data='010_data/23CD4_02CD4_expr_raw.rds'
path_cytokines_cutoffs='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_panels/panel2CD4_cytokines_CM.xlsx'
path_clustering='030_heatmaps/23CD4_02CD4_pca1_merging2_clustering.xls'
path_clustering_labels='030_heatmaps/23CD4_02CD4_pca1_merging2_clustering_labels.xls'
clsubset=c('CM','EM','TM','TE')
cutoff_colname=c('positive_cutoff_raw_base','positive_cutoff_raw_tx')

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

setwd(rwd)

prefix <- pd1_prefix
suffix <- ""
outdir <- pd1_outdir

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
# Indeces of observables used for PD1 and positive-negative analysis

pd1_marker <- cytokines_cutoffs$Antigen == "PD-1"

if(any(is.na(cytokines_cutoffs[pd1_marker, cutoff_colname])))
  stop("NAs in PD-1 cutoffs")

pn_markers <- complete.cases(cytokines_cutoffs[, cutoff_colname, drop = FALSE]) & !pd1_marker


pd1col <- which(fcs_colnames %in% cytokines_cutoffs[pd1_marker, "fcs_colname"])

pncols <- which(fcs_colnames %in% cytokines_cutoffs[pn_markers, "fcs_colname"])

## Order them in the way as in panel
mm <- match(cytokines_cutoffs[pn_markers, "Antigen"], fcs_panel$Antigen[pncols])
pncols <- pncols[mm]


# ------------------------------------------------------------
# Keep only the data from the specified clusters

cells2keep_clust <- clust %in% labels[labels$label %in% clsubset, "cluster"]


# ----------------------------------------------------------------------------------------
# Create the bimatrix - for PD1
# ----------------------------------------------------------------------------------------

epd1 <- e[cells2keep_clust, pd1col, drop = FALSE]
sampb <- samp[cells2keep_clust]

## get the corresponding cutoffs
mm <- match(colnames(epd1), cytokines_cutoffs$fcs_colname)

cytcut <- cytokines_cutoffs[mm, cutoff_colname, drop = FALSE]
rownames(cytcut) <- colnames(epd1)
print(cytcut)


### create the bimatrix

## use one cutoff
if(length(cutoff_colname) == 1){
  
  bimatrix_pd1 <- epd1 > cytcut[, cutoff_colname]
  
}

## use base and tx cutoffs
if(length(cutoff_colname) == 2){
  
  cutoff_colname_base <- cutoff_colname[grep("base", cutoff_colname)]
  cutoff_colname_tx <- cutoff_colname[grep("tx", cutoff_colname)]
  
  cytcut_samp <- cytcut[, ifelse(grepl("base", sampb), cutoff_colname_base, cutoff_colname_tx)]
  
  bimatrix_pd1 <- epd1 > as.numeric(cytcut_samp)
  
}


bimatrix_pd1 <- apply(bimatrix_pd1, 2, as.numeric)


# ------------------------------------------------------------
# Create a table with clustering and clustering_labels - needed for 04_frequencies.R to work 
# ------------------------------------------------------------

clustering <- data.frame(cluster = bimatrix_pd1[, 1] + 1, cell_id = cell_id[cells2keep_clust], sample_id = samp[cells2keep_clust], stringsAsFactors = FALSE)

write.table(clustering, file.path(outdir, paste0(prefix, "clustering.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

freq_clust <- table(clustering$cluster)

clustering_labels <- data.frame(cluster = c(1, 2), label = c("PD1-", "PD1+"), counts = as.numeric(freq_clust), stringsAsFactors = FALSE)

write.table(clustering_labels, file.path(outdir, paste0(prefix, "clustering_labels.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)




# ------------------------------------------------------------
# Save the expression of cytokines for Tmem cells
# ------------------------------------------------------------

exprb <- expr[cells2keep_clust, c("cell_id", "sample_id", fcs_colnames[pncols])]

saveRDS(exprb, file.path(outdir, paste0(prefix, "expr_raw.rds")))


# ------------------------------------------------------------
# Create a table with clustering observables - needed for 04_expression.R to work 
# ------------------------------------------------------------


clustering_observables <- data.frame(mass = fcs_panel$fcs_colname[pncols], marker = fcs_panel$Antigen[pncols], clustering_observable = TRUE, stringsAsFactors = FALSE)

write.table(clustering_observables, file.path(outdir, paste0(prefix, "clustering_observables.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)































################################
### 07_pd1_expression.R done!
################################