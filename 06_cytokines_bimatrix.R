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
# cytokines_prefix='23CD4_02CD4_pca1_merging_Tmem_cytCM_raw_'
# cytokines_outdir='060_cytokines'
# path_data='010_data/23CD4_02CD4_expr_raw.rds'
# path_cytokines_cutoffs='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_panels/panel2CD4_cytokines_CM.xlsx'
# path_clustering='030_heatmaps/23CD4_02CD4_pca1_merging_clustering.xls'
# path_clustering_labels='030_heatmaps/23CD4_02CD4_pca1_merging_clustering_labels.xls'
# clsubset=c('CM','EM','TE')
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
# Load cytokines_cutoffs
# ------------------------------------------------------------

## positive cutoffs for cytokines
cytokines_cutoffs <- read.xls(path_cytokines_cutoffs, stringsAsFactors=FALSE)

if(!cutoff_colname %in% colnames(cytokines_cutoffs))
  stop("There are no such column with cutoffs!")


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
# Do the subsetting of cells and markers used for the analysis
# ------------------------------------------------------------


# -------------------------------------
# Get the isotope and antigen for fcs markers

m <- match(fcs_colnames, cytokines_cutoffs$fcs_colname)

fcs_panel <- data.frame(fcs_colname = fcs_colnames, Isotope = cytokines_cutoffs$Isotope[m], Antigen = cytokines_cutoffs$Antigen[m], stringsAsFactors = FALSE)


# -------------------------------------
# Indeces of observables used for positive-negative analysis

pncols <- which(fcs_colnames %in% cytokines_cutoffs[!is.na(cytokines_cutoffs[, cutoff_colname]), "fcs_colname"])

## Order them in the way as in panel
mm <- match(cytokines_cutoffs[!is.na(cytokines_cutoffs[, cutoff_colname]), "Antigen"], fcs_panel$Antigen[pncols])
pncols <- pncols[mm]


# ------------------------------------------------------------
# Keep only the data from the specified clusters

cells2keep_clust <- clust %in% labels[labels$label %in% clsubset, "cluster"]

eb <- e[cells2keep_clust, pncols]

# ------------------------------------------------------------
# Create the bimatrix - TRUE when cells are positive expressed for a given marker
# ------------------------------------------------------------

## get the corresponding cutoffs
mm <- match(colnames(eb), cytokines_cutoffs$fcs_colname)

ecut <- cytokines_cutoffs[mm, cutoff_colname]
names(ecut) <- colnames(eb)

ecut


## create the bimatrix
bimatrix <- t(t(eb) > ecut)
bimatrix <- apply(bimatrix, 2, as.numeric)

## Keep only cells that are positive for at least one marker
cells2keep_pos <- rowSums(bimatrix) > 0

bimatrix_pos <- bimatrix[cells2keep_pos, ]



### Save the bmatrix (it has to contain cell_id and sample_id)

bimatrix_out <- data.frame(cell_id = cell_id[cells2keep_clust][cells2keep_pos], sample_id = samp[cells2keep_clust][cells2keep_pos], bimatrix_pos, check.names = FALSE)


write.table(bimatrix_out, file.path(outdir, paste0(prefix, "bimatrix", suffix, ".txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



# ------------------------------------------------------------
# Upsetr plots
# ------------------------------------------------------------

bidf <- data.frame(bimatrix, row.names = 1:nrow(bimatrix), check.names = FALSE)

pdf(file.path(outdir, paste0(prefix, "upsetr", suffix, ".pdf")), w = 16, h = 6)
upset(bidf, sets = colnames(bidf), nintersects = 50, order.by = "freq")
dev.off()




# ------------------------------------------------------------
# Create a table with clustering observables - needed for 02_flowsom.R to work 
# ------------------------------------------------------------

clustering_observables <- data.frame(mass = fcs_panel$fcs_colname[pncols], marker = fcs_panel$Antigen[pncols], clustering_observable = TRUE, stringsAsFactors = FALSE)

write.table(clustering_observables, file.path(outdir, paste0(prefix, "clustering_observables.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


















################################
### 06_cytokines_bimatrix.R done!
################################