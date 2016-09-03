##############################################################################
## <<03_plottsne_expr.R>>

# BioC 3.3
# Created 22 Aug 2016
# Updated 30 Aug 2016

##############################################################################
Sys.time()
##############################################################################

# Load packages
library(gdata)
library(Repitools)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2) # melt
library(RColorBrewer)
library(Rtsne)
library(coop) # cosine

##############################################################################
# Test arguments
##############################################################################

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
# tsnep_prefix='23_01_pca1_raw_mark_raw_'
# tsnep_outdir='040_tsnemaps_expr'
# path_data='010_data/23_01_expr_raw.rds'
# path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.xlsx'
# path_rtsne_out='040_tsnemaps/23_01_pca1_raw_rtsne_out.rda'
# path_rtsne_data='040_tsnemaps/23_01_pca1_raw_rtsne_data.xls'
# path_clustering_observables='030_heatmaps/23_01_pca1_clustering_observables.xls'
# pdf_width=15
# pdf_height=10

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
suffix <- ""
outdir <- tsnep_outdir

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
# Load metadata
# ------------------------------------------------------------

md <- read.xls(path_metadata, stringsAsFactors=FALSE)


### Colors 
colors <- unique(md[, c("condition", "color")])
colors$condition <- factor(colors$condition)
## replace _ with \n
levels(colors$condition) <- gsub("_", "\n", levels(colors$condition ))

color_groups <- colors$color
names(color_groups) <- colors$condition

color_samples <- md$color
names(color_samples) <- md$shortname


# ------------------------------------------------------------
# Load clustering data
# ------------------------------------------------------------

# clustering observables
clustering_observables <- read.table(path_clustering_observables, header = TRUE, sep = "\t", as.is = TRUE)
rownames(clustering_observables) <- clustering_observables$mass

clust_observ <- clustering_observables[clustering_observables$clustering_observable, "mass"]


# ------------------------------------------------------------
# get the isotope and antigen for fcs markers

m <- match(fcs_colnames, clustering_observables$mass)

fcs_panel <- data.frame(fcs_colname = fcs_colnames, Isotope = clustering_observables$mass[m], Antigen = clustering_observables$marker[m], stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Load tSNE data
# ------------------------------------------------------------

load(path_rtsne_out)

rtsne_data <- read.table(path_rtsne_data, header = TRUE, sep = "\t", as.is = TRUE)


# ------------------------------------------------------------
# tSNE plots with marker expression as a heat color
# ------------------------------------------------------------

## get expression for cells used in tSNE
rownames(e) <- expr$cell_id
expr_tsne <- e[rtsne_data$cell_index, fcs_panel$fcs_colname]

colnames(expr_tsne) <- fcs_panel$Antigen

marker_selection <- colnames(expr_tsne)


for(i in 1:length(marker_selection)){
  # i = 1
  
  ggdf <- data.frame(tSNE1 = rtsne_out$Y[,1], tSNE2 = rtsne_out$Y[,2], sample = rtsne_data$sample_name, marker_selection = expr_tsne[, marker_selection[i]])
  
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
  
  pdf(file.path(outdir, paste0(prefix, gsub("[[:punct:]]", "", marker_selection[i]), "_tSNEgroup", suffix, ".pdf")), width = pdf_width, height = pdf_height)                 
  print(ggp)
  dev.off()
  
  
  ## one plot 
  ggp <- ggplot(ggdf,  aes(x = tSNE1, y = tSNE2, color = marker_selection)) +
    geom_point(size = 1) +
    labs(x = "tSNE 1", y="tSNE 2")+ 
    theme_bw() +
    theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold")) +
    scale_colour_gradientn(marker_selection[i], colours = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(100))
  
  pdf(file.path(outdir, paste0(prefix, gsub("[[:punct:]]", "", marker_selection[i]), "_tSNEone", suffix, ".pdf")), width = 9, height = 7)      
  print(ggp)
  dev.off()
  
  
}


# ------------------------------------------------------------
# tSNE plots with samples as a heat color
# ------------------------------------------------------------


ggdf <- data.frame(tSNE1 = rtsne_out$Y[,1], tSNE2 = rtsne_out$Y[,2], sample = rtsne_data$sample_name)

# add group info
mm <- match(ggdf$sample, md$shortname)
ggdf$group <- md$condition[mm]


### Plot of tsne

## facet per group
ggp <- ggplot(ggdf,  aes(x = tSNE1, y = tSNE2, color = sample)) +
  geom_point(size=1) +
  facet_wrap(~ group) +
  labs(x = "tSNE 1", y="tSNE 2")+ 
  theme_bw() +
  theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold")) +
  scale_color_manual(values = color_samples)

pdf(file.path(outdir, paste0(prefix, "sample_tSNEgroup", suffix, ".pdf")), width = pdf_width, height = pdf_height)    
print(ggp)
dev.off()


## one plot 
ggp <- ggplot(ggdf,  aes(x = tSNE1, y = tSNE2, color = sample)) +
  geom_point(size = 1) +
  labs(x = "tSNE 1", y="tSNE 2")+ 
  theme_bw() +
  theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold")) +
  scale_color_manual(values = color_samples)

pdf(file.path(outdir, paste0(prefix, "sample_tSNEone", suffix, ".pdf")), width = 9, height = 7)      
print(ggp)
dev.off()





sessionInfo()















################################
### 03_plottsne_expr done!
################################