##############################################################################
## <<02_select_observables.R>>

# BioC 3.3
# Created 27 July 2016


##############################################################################
Sys.time()
##############################################################################

# Load packages

##############################################################################
# Test arguments
##############################################################################

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
# observ_prefix='pca1_'
# pca_score_cutoff=3
# pca_skip_top=0
# path_pca_score


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

prefix <- observ_prefix

# ------------------------------------------------------------
# define directories
# ------------------------------------------------------------

fcsDir <- "010_cleanfcs"; if( !file.exists(fcsDir) ) dir.create(fcsDir)
pcaDir <- "020_pcascores"; if( !file.exists(pcaDir) ) dir.create(pcaDir)
hmDir <- "030_heatmaps"; if( !file.exists(hmDir) ) dir.create(hmDir)


# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------


# ------------------------------------------------------------
# Load more data
# ------------------------------------------------------------


if(!grepl("/", path_pca_score)){
  prs <- read.table(file.path(pcaDir, path_pca_score), header = TRUE, sep = "\t", as.is = TRUE)
}else{
  prs <- read.table(file.path(path_pca_score), header = TRUE, sep = "\t", as.is = TRUE)
}

prs <- prs[order(prs$avg_score, decreasing = TRUE), ]


# -------------------------------------
# selected observables with PCA scores higher than a cutoff but skipp top X scores 
# -------------------------------------


clustering_observables <- prs[, c("mass", "marker", "avg_score")]

clustering_observables$clustering_observable <- prs$avg_score > pca_score_cutoff


if(pca_skip_top > 0)
  clustering_observables$clustering_observable[1:pca_skip_top] <- FALSE



### Save the observables

write.table(clustering_observables, file = file.path(hmDir, paste0(prefix, "clustering_observables.xls")), row.names=FALSE, quote=FALSE, sep="\t")




sessionInfo()


























################################
### 02_select_observables done!
################################