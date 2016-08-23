##############################################################################
## <<03_plottsne_expr.R>>

# BioC 3.3
# Created 22 Aug 2016
# Updated 

##############################################################################
Sys.time()
##############################################################################

# Load packages
library(flowCore)
library(gdata)
library(FlowSOM)
library(Repitools)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2)
library(coop)
library(RColorBrewer)
library(Rtsne)
library(coop) # cosine

##############################################################################
# Test arguments
##############################################################################

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_03'
# path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_03.xlsx'
# tsnep_prefix='23_03_pca1_mark_'
# path_rtsne_out='23_03_pca1_rtsne_out_norm.rda'
# path_rtsne_data='23_03_pca1_rtsne_data_norm.xls'
# path_clustering_observables='23_03_pca1_clustering_observables.xls'
# tsne_cmin=1000
# pdf_width=15
# pdf_height=10
# tsnep_suffix='_norm'
# marker_selection=c('ICAM-1','CD274_PDL1','CD86')

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

prefix <- tsnep_prefix
suffix <- tsnep_suffix


# ------------------------------------------------------------
# define directories
# ------------------------------------------------------------

fcsDir <- "010_cleanfcs"; if( !file.exists(fcsDir) ) dir.create(fcsDir)
pcaDir <- "020_pcascores"; if( !file.exists(pcaDir) ) dir.create(pcaDir)
hmDir <- "030_heatmaps"; if( !file.exists(hmDir) ) dir.create(hmDir)
sneDir <- "040_tsnemaps"; if( !file.exists(sneDir) ) dir.create(sneDir)


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

# ------------------------------------------------------------
# Load cluster data
# ------------------------------------------------------------


# clustering_observables
if(!grepl("/", path_clustering_observables)){
  clustering_observables <- read.table(file.path(hmDir, path_clustering_observables), header = TRUE, sep = "\t", as.is = TRUE)
}else{
  clustering_observables <- read.table(path_clustering_observables, header = TRUE, sep = "\t", as.is = TRUE)
}

rownames(clustering_observables) <- clustering_observables$mass


# marker_selection

if(!all(marker_selection %in% clustering_observables$marker))
  stop("Marker selection is invalid!")


# ------------------------------------------------------------
# Load tSNE data
# ------------------------------------------------------------

load(file.path(sneDir, path_rtsne_out))

rtsne_data <- read.table(file.path(sneDir, path_rtsne_data), header = TRUE, sep = "\t", as.is = TRUE)


# -------------------------------------
# get the isotope and antigen for fcs markers

m <- match(fcs_colnames, clustering_observables$mass)

fcs_panel <- data.frame(colnames = fcs_colnames, Isotope = clustering_observables$mass[m], Antigen = clustering_observables$marker[m], stringsAsFactors = FALSE)


# -------------------------------------

cols <- which(fcs_colnames %in% clustering_observables[clustering_observables$marker %in% marker_selection, "mass"])

# arc-sin-h the columns specific 
fcsT <- lapply(fcs, function(u) {
  e <- exprs(u)
  e[ , cols] <- asinh( e[ , cols] / 5 )
  exprs(u) <- e
  u
})



# ------------------------------------------------------------
# Get expression
# ------------------------------------------------------------

### Get marker expression

es <- lapply(fcsT, function(u) {
  exprs(u)[, cols]
})

e <- do.call("rbind",es)
colnames(e) <- fcs_panel$Antigen[cols]


# normalize to 0-1
el <- e
rng <- apply(el,2,quantile,p=c(.01,.99))
for(i in 1:ncol(el)) {
  el[,i] <- (el[,i]-rng[1,i])/(rng[2,i]-rng[1,i])
}
el[el<0] <- 0
el[el>1] <- 1



# ------------------------------------------------------------
# tSNE plots
# ------------------------------------------------------------

## get expression for cells used in tSNE
expr <- el[rtsne_data$cell_index, ]


for(i in 1:length(marker_selection)){
  # i = 1
  
  ggdf <- data.frame(tSNE1 = rtsne_out$Y[,1], tSNE2 = rtsne_out$Y[,2], sample = rtsne_data$sample_name, marker_selection = expr[, marker_selection[i]])
  
  # add group info
  mm <- match(ggdf$sample, md$shortname)
  ggdf$group <- md$condition[mm]
  
  
  ### Plot of tsne
  
  ## facet per group
  ggp <- ggplot(ggdf,  aes(x = tSNE1, y = tSNE2, color = marker_selection)) +
    geom_point(size=1) +
    facet_wrap(~ group) +
    labs(x = "tSNE 1", y="tSNE 2")+ 
    theme_bw() +
    theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold")) +
    scale_colour_gradientn(marker_selection[i], colours = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(100))
  
  pdf(file.path(sneDir, paste0(prefix, marker_selection[i], "_norm_", "tSNE", suffix, ".pdf")), width = pdf_width, height = pdf_height)                 
  print(ggp)
  dev.off()
  
  
  ## one plot 
  ggp <- ggplot(ggdf,  aes(x = tSNE1, y = tSNE2, color = marker_selection)) +
    geom_point(size = 1) +
    labs(x = "tSNE 1", y="tSNE 2")+ 
    theme_bw() +
    theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold")) +
    scale_colour_gradientn(marker_selection[i], colours = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(100))
  
  pdf(file.path(sneDir, paste0(prefix, marker_selection[i], "_norm_", "tSNEone", suffix, ".pdf")), width = 9, height = 7)                 
  print(ggp)
  dev.off()
  
  
}








sessionInfo()















################################
### 03_plottsne_expr done!
################################