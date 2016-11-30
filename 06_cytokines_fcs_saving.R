##############################################################################
## <<06_cytokines_fcs_saving.R>>

# BioC 3.3
# Created 30 Nov 2016

##############################################################################
Sys.time()
##############################################################################

library(flowCore)
library(gdata)
library(Repitools)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2)


##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-merged_23_29/02_CD4'
save_prefix='23_29_CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_cl40_data23_'
save_outdir='10_cytokines_merged_top_combinations/06_dumpfcs'

path_panel='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_panels/panel2CD4_23_cytokines_CM.xlsx'
path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_02.xlsx'

path_fcs='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_02_CD4_merging2/010_cleanfcs'
path_clustering='10_cytokines_merged_top_combinations/01_clustering/23_29_CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_cl40_data23_clustering.xls'
path_clustering_labels='10_cytokines_merged_top_combinations/01_clustering/23_29_CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_cl40_data23_clustering_labels.xls'


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

fcsDir <- path_fcs

outdir <- save_outdir
if(!file.exists(outdir)) 
  dir.create(outdir, recursive = TRUE)

prefix <- save_prefix

# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

# read metadata
md <- read.xls(path_metadata, stringsAsFactors=FALSE)

md$condition <- factor(md$condition)



# ------------------------------------------------------------
# Load cluster data
# ------------------------------------------------------------

## clustering labels
labels <- read.table(path_clustering_labels, header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))

## clustering
clustering <- read.table(path_clustering, header = TRUE, sep = "\t", as.is = TRUE)

## drop the "drop" cluster
if("drop" %in% labels$label){
  
  clust2drop <- labels$cluster[labels$label == "drop"]
  cells2drop <- clustering$cluster != clust2drop
  clustering <- clustering[cells2drop, , drop = FALSE]
  labels <- labels[labels$label != "drop", ,drop = FALSE]
  labels$label <- factor(labels$label)
  
}



# ------------------------------------------------------------
# Load fcs files
# ------------------------------------------------------------

# define FCS file names
f <- file.path(fcsDir, md$filename)
names(f) <- md$shortname

# read raw FCS files in
fcs <- lapply(f, read.FCS)

fcs_colnames <- colnames(fcs[[1]])
fcs_colnames

## Create sample info
samp <- rep(names(fcs), sapply(fcs, nrow))

cells <- data.frame(sample_id = samp, cell_id = 1:length(samp), stringsAsFactors = FALSE)


# ------------------------------------------------------------
# Load panel
# ------------------------------------------------------------

# read panel, pick which columns to use
panel <- read.xls(path_panel, stringsAsFactors=FALSE)

if(!all(c("fcs_colname", "Isotope", "Antigen", "transform") %in% colnames(panel)))
  stop("Wrong columns in panel!!!")


# cols - get fcs columns that are in the panel with transform = 1
cols <- which(fcs_colnames %in% panel$fcs_colname[panel$transform==1])

# get the isotope and antigen for fcs markers
m <- match(fcs_colnames, panel$fcs_colname)

fcs_panel <- data.frame(fcs_colname = fcs_colnames, Isotope = panel$Isotope[m], Antigen = panel$Antigen[m], stringsAsFactors = FALSE)


# --------------------------------------------------------------------------
# Subset the fcs files 
# --------------------------------------------------------------------------

# arc-sin-h the columns specific 
fcsT <- lapply(fcs, function(u) {
  e <- exprs(u)
  e[,cols] <- asinh( e[,cols] / 5 )
  exprs(u) <- e
  u
})


cells_split <- split(cells, factor(cells$sample_id, levels = names(fcs)))
clustering_split <- split(clustering, factor(clustering$sample_id, levels = names(fcs)))

conditions <- levels(md$condition)


for(i in 1:nrow(labels)){
  # i = 1
  
  extract_cluster <- labels[i, "cluster"]
  
  for(j in 1:nlevels(md$condition)){
    # j = 1
    
    samp_sub <- md$shortname[md$condition == conditions[j]]
    
    fcs_sub <- lapply(samp_sub, function(k){
      # k = "base_HD1"
      
      fcs_tmp <- fcsT[[k]][cells_split[[k]][, "cell_id"] %in% clustering_split[[k]][clustering_split[[k]][, "cluster"] == extract_cluster, "cell_id"], ]
      
      fcs_tmp
      
    })

    fcs_sub <- as(as(fcs_sub, "flowSet"), "flowFrame")
    
    ### Save the new fcs 
    write.FCS(fcs_sub, file.path(outdir, paste0(prefix, "arcsineh_cl", extract_cluster, "_", conditions[j], ".fcs")))
    
  }
  
}








sessionInfo()






################################
### 06_cytokines_fcs_saving.R done!
################################