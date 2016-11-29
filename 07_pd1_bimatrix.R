##############################################################################
## <<07_pd1_bimatrix.R>>

# BioC 3.3
# Created 27 Aug 2016
# Updated 31 Aug 2016

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
pd1_prefix='23CD4_02CD4_pca1_merging_Tmem_cytCM_raw2_pd1_'
pd1_outdir='070_pd1'
path_data='010_data/23CD4_02CD4_expr_raw.rds'
path_cytokines_cutoffs='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_panels/panel2CD4_cytokines_CM.xlsx'
path_clustering='030_heatmaps/23CD4_02CD4_pca1_merging_clustering.xls'
path_clustering_labels='030_heatmaps/23CD4_02CD4_pca1_merging_clustering_labels.xls'
clsubset=c('CM','EM')
cutoff_colname=c('positive_cutoff_raw_base','positive_cutoff_raw_tx')
marker='PD-1'


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

if(length(marker) != 1 || !any(marker %in% cytokines_cutoffs$Antigen))
  stop("Marker name is wrong!")

pd1_marker <- cytokines_cutoffs$Antigen == marker

if(any(is.na(cytokines_cutoffs[pd1_marker, cutoff_colname])))
  stop(paste0("NAs in ", marker," cutoffs"))

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

### Save the bmatrix (it has to contain cell_id and sample_id)

bimatrix_out <- data.frame(cell_id = cell_id[cells2keep_clust], sample_id = samp[cells2keep_clust], bimatrix_pd1)

write.table(bimatrix_out, file.path(outdir, paste0(prefix, "bimatrix", suffix, ".txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# ------------------------------------------------------------
# Create a table with clustering and clustering_labels - needed for 04_frequencies.R to work 
# ------------------------------------------------------------

clustering <- data.frame(cluster = bimatrix_pd1[, 1] + 1, cell_id = cell_id[cells2keep_clust], sample_id = samp[cells2keep_clust], stringsAsFactors = FALSE)

write.table(clustering, file.path(outdir, paste0(prefix, "clustering.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

freq_clust <- table(clustering$cluster)

clustering_labels <- data.frame(cluster = c(1, 2), label = paste0(marker, c("-", "+")), counts = as.numeric(freq_clust), stringsAsFactors = FALSE)

write.table(clustering_labels, file.path(outdir, paste0(prefix, "clustering_labels.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# ----------------------------------------------------------------------------------------
# Create a bimatrix for other cytokines for positive and negative PD1
# ----------------------------------------------------------------------------------------

eb <- e[cells2keep_clust, pncols, drop = FALSE]
sampb <- samp[cells2keep_clust]


## get the corresponding cutoffs
mm <- match(colnames(eb), cytokines_cutoffs$fcs_colname)

cytcut <- cytokines_cutoffs[mm, cutoff_colname, drop = FALSE]
rownames(cytcut) <- colnames(eb)
print(cytcut)



### Treat PD-1 positive or PD-1 negative cells as 100%
cells2keep_pd1 <- list()

cells2keep_pd1[["positive"]] <- bimatrix_pd1 == 1
cells2keep_pd1[["negative"]] <- bimatrix_pd1 == 0

pd1_type <- names(cells2keep_pd1)
keep_pos_cells <- c(FALSE, TRUE)

for(i in 1:length(pd1_type)){
  # i = 1
  
  prefix_type <- paste0(pd1_type[i], "_")
  
  eb_tmp <- eb[cells2keep_pd1[[i]], ]
  sampb_tmp <- sampb[cells2keep_pd1[[i]]]
  
  ### create the bimatrix
  
  ## use one cutoff
  if(length(cutoff_colname) == 1){
    
    bimatrix <- t(t(eb_tmp) > cytcut[, cutoff_colname])
    
  }
  
  ## use base and tx cutoffs
  if(length(cutoff_colname) == 2){
    
    cutoff_colname_base <- cutoff_colname[grep("base", cutoff_colname)]
    cutoff_colname_tx <- cutoff_colname[grep("tx", cutoff_colname)]
    
    cytcut_samp <- cytcut[, ifelse(grepl("base", sampb_tmp), cutoff_colname_base, cutoff_colname_tx)]
    
    bimatrix <- t(t(eb_tmp) > cytcut_samp)
    
  }
  
  bimatrix <- apply(bimatrix, 2, as.numeric)
  
  ### Keep only cells that are positive for at least one marker
  
  # for PD1+, I keep all the cells
  cells2keep_pos <- rep(TRUE, nrow(bimatrix))
  bimatrix_pos <- bimatrix
  
  # for PD1-, I keep the cells that are positive for at least one marker
  if(keep_pos_cells[i]){
    cells2keep_pos <- rowSums(bimatrix) > 0
    table(cells2keep_pos)
    bimatrix_pos <- bimatrix[cells2keep_pos, ]
  }
  

  ### Save the bmatrix (it has to contain cell_id and sample_id)
  
  bimatrix_out <- data.frame(cell_id = cell_id[cells2keep_clust][cells2keep_pd1[[i]]][cells2keep_pos], sample_id = samp[cells2keep_clust][cells2keep_pd1[[i]]][cells2keep_pos], bimatrix_pos)
  
  write.table(bimatrix_out, file.path(outdir, paste0(prefix, prefix_type, "bimatrix", suffix, ".txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  

  # ------------------------------------------------------------
  # Upsetr plots
  # ------------------------------------------------------------
  
  bidf <- data.frame(bimatrix, row.names = 1:nrow(bimatrix), check.names = FALSE)
  colnames(bidf) <- fcs_panel$Antigen[pncols]
  
  pdf(file.path(outdir, paste0(prefix, prefix_type, "upsetr", suffix, ".pdf")), w = 16, h = 6)
  upset(bidf, sets = colnames(bidf), nintersects = 50, order.by = "freq")
  dev.off()
  
}



# ------------------------------------------------------------
# Create a table with clustering observables - needed for 02_flowsom.R to work 
# ------------------------------------------------------------

clustering_observables <- data.frame(mass = fcs_panel$fcs_colname[pncols], marker = fcs_panel$Antigen[pncols], clustering_observable = TRUE, stringsAsFactors = FALSE)

write.table(clustering_observables, file.path(outdir, paste0(prefix, "clustering_observables.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

for(i in 1:length(pd1_type))
  write.table(clustering_observables, file.path(outdir, paste0(prefix, pd1_type[i], "_", "clustering_observables.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)














################################
### 07_pd1_bimatrix.R done!
################################