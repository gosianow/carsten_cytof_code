##############################################################################
## <<02_select_observables.R>>

# BioC 3.3
# Created 27 July 2016
# Updated 

##############################################################################
Sys.time()
##############################################################################

# Load packages
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

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
# pca_score_cutoff=3
# pca_skip_top=0
# observ_prefix='pca1_'

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

# read metadata
md <- read.xls("metadata.xlsx",stringsAsFactors=FALSE)

# read panel, pick which columns to use
panel <- read.xls("panel.xlsx",stringsAsFactors=FALSE)


# define FCS file names
f <- file.path(fcsDir, md$filename)
names(f) <- md$shortname


# read raw FCS files in
fcs <- lapply(f, read.FCS)


# get isotope mass of columns in fcs files.. to match against the panel
panel_mass <- as.numeric(gsub("[[:alpha:]]", "", colnames(fcs[[1]])))


# cols - get fcs columns that are in the panel with transform = 1
cols <- which(panel_mass %in% panel$Isotope[panel$transform==1])

# Antigen - get the antigen name
m <- match(panel_mass, panel$Isotope)

fcs_panel <- data.frame(colnames = colnames(fcs[[1]]), Isotope = panel_mass, cols = panel_mass %in% panel$Isotope[panel$transform==1], Antigen = panel$Antigen[m], stringsAsFactors = FALSE)

fcs_panel$Antigen[is.na(fcs_panel$Antigen)] <- ""


# arc-sin-h the columns specific 
fcsT <- lapply(fcs, function(u) {
  e <- exprs(u)
  e[,cols] <- asinh( e[,cols] / 5 )
  exprs(u) <- e
  u
})



# ------------------------------------------------------------
# Load more data
# ------------------------------------------------------------


prs <- read.table(file.path(pcaDir,"princompscore_by_sample.xls"), header = TRUE, sep = "\t", as.is = TRUE)
prs <- prs[order(prs$avg_score, decreasing = TRUE), ]


# -------------------------------------
# selected observables with PCA scores higher than a cutoff but skipp top X scores 
# -------------------------------------


if(pca_skip_top > 0)
  prs <- prs[-c(1:pca_skip_top), ]

scols <- which(colnames(fcsT[[1]]) %in% prs[prs$avg_score > pca_score_cutoff, "mass"])

clustering_observables <- colnames(fcsT[[1]])[scols]


### Save the observables

scols_out <- data.frame(clustering_observables = clustering_observables, stringsAsFactors = FALSE)
write.table(scols_out, file = file.path(hmDir, paste0(prefix, "clustering_observables.xls")), row.names=FALSE, quote=FALSE, sep="\t")




sessionInfo()













################################
### Done!
################################