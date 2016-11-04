##############################################################################
## <<06_cytokines_bimatrix.R>>

# BioC 3.3
# Created 24 Aug 2016
# Updated 25 Aug 2016

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
# cytokines_prefix='23CD4_02CD4_pca1_merging_Tmem_cytCM_raw2_'
# cytokines_outdir='060_cytokines_bimatrix'
# path_data='010_data/23CD4_02CD4_expr_raw.rds'
# path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_02.xlsx'
# path_cytokines_cutoffs='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_panels/panel2CD4_cytokines_CM.xlsx'
# path_clustering='030_heatmaps/23CD4_02CD4_pca1_merging_clustering.xls'
# path_clustering_labels='030_heatmaps/23CD4_02CD4_pca1_merging_clustering_labels.xls'
# clsubset=c('CM','EM')
# cutoff_colname=c('positive_cutoff_raw_base','positive_cutoff_raw_tx')

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
  dir.create(outdir, recursive = TRUE)

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
# Load metadata
# ------------------------------------------------------------

md <- read.xls(path_metadata, stringsAsFactors=FALSE)

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

eb <- e[cells2keep_clust, pncols]
sampb <- samp[cells2keep_clust]

# ------------------------------------------------------------
# Create the bimatrix - TRUE when cells are positive expressed for a given marker
# ------------------------------------------------------------

## get the corresponding cutoffs
mm <- match(colnames(eb), cytokines_cutoffs$fcs_colname)

cytcut <- cytokines_cutoffs[mm, cutoff_colname, drop = FALSE]
rownames(cytcut) <- colnames(eb)
print(cytcut)


### create the bimatrix

## use one cutoff
if(length(cutoff_colname) == 1){
  
  bimatrix <- t(t(eb) > cytcut[, cutoff_colname])
  
}

## use base and tx cutoffs
if(length(cutoff_colname) == 2){
  
  cutoff_colname_base <- cutoff_colname[grep("base", cutoff_colname)]
  cutoff_colname_tx <- cutoff_colname[grep("tx", cutoff_colname)]
  
  cytcut_samp <- cytcut[, ifelse(grepl("base", sampb), cutoff_colname_base, cutoff_colname_tx)]
  
  bimatrix <- t(t(eb) > cytcut_samp)
  
}

bimatrix <- apply(bimatrix, 2, as.numeric)
bm <- bimatrix

## Keep only cells that are positive for at least one marker
cells2keep_pos <- rowSums(bimatrix) > 0
table(cells2keep_pos)


bimatrix_pos <- bimatrix[cells2keep_pos, ]



### Save the bmatrix (it has to contain cell_id and sample_id)

bimatrix_out <- data.frame(cell_id = cell_id[cells2keep_clust][cells2keep_pos], sample_id = samp[cells2keep_clust][cells2keep_pos], bimatrix_pos, check.names = FALSE)


write.table(bimatrix_out, file.path(outdir, paste0(prefix, "bimatrix", suffix, ".txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



# ------------------------------------------------------------
# Upsetr plots
# ------------------------------------------------------------

bidf <- data.frame(bimatrix, row.names = 1:nrow(bimatrix), check.names = FALSE)
colnames(bidf) <- fcs_panel$Antigen[pncols]

pdf(file.path(outdir, paste0(prefix, "upsetr", suffix, ".pdf")), w = 16, h = 6)
upset(bidf, sets = colnames(bidf), nintersects = 50, order.by = "freq")
dev.off()


# ------------------------------------------------------------
# Create a table with clustering observables - needed for 02_flowsom.R to work 
# ------------------------------------------------------------

clustering_observables <- data.frame(mass = fcs_panel$fcs_colname[pncols], marker = fcs_panel$Antigen[pncols], clustering_observable = TRUE, stringsAsFactors = FALSE)

write.table(clustering_observables, file.path(outdir, paste0(prefix, "clustering_observables.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# ------------------------------------------------------------
# Save the subsets of bimatrix for diff. conditions separately (needed for plotting the heatmaps per condition) 
# ------------------------------------------------------------

# mm <- match(bimatrix_out$sample_id, md$shortname)
# 
# split_condition <- factor(md$condition[mm])
# bimatrix_split <- split(bimatrix_out, split_condition)
# 
# split_levels <- levels(split_condition)
# 
# dummy <- lapply(1:length(split_levels), function(i){
#   # i = 1
#   
#   bimatrix_split_out <- bimatrix_split[[split_levels[i]]]
#   
#   write.table(bimatrix_split_out, file.path(outdir, paste0(prefix, "bimatrix_", split_levels[i], suffix, ".txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
#   
#   return(NULL)
#   
# })


































################################
### 06_cytokines_bimatrix.R done!
################################