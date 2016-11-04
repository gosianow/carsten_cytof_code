##############################################################################
## <<03_plottsne.R>>

# BioC 3.3
# Created 28 July 2016
# Updated 

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
# tsnep_prefix='23_01_pca1_cl20_raw_'
# tsnep_outdir='040_tsnemaps'
# path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.xlsx'
# path_rtsne_out='040_tsnemaps/23_01_pca1_raw_rtsne_out.rda'
# path_rtsne_data='040_tsnemaps/23_01_pca1_raw_rtsne_data.xls'
# path_clustering='030_heatmaps/23_01_pca1_cl20_clustering.xls'
# path_clustering_labels='030_heatmaps/23_01_pca1_cl20_clustering_labels.xls'
# tsne_cmin=1000
# tsne_distse=1
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
  dir.create(outdir, recursive = TRUE)


# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

md <- read.xls(path_metadata, stringsAsFactors=FALSE)

# ------------------------------------------------------------
# Load cluster data
# ------------------------------------------------------------

## clustering
clustering <- read.table(path_clustering, header = TRUE, sep = "\t", as.is = TRUE)

clust <- clustering[, "cluster"]

## clustering labels
labels <- read.table(path_clustering_labels, header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))


# ------------------------------------------------------------
### Colors for tSNE maps
# ------------------------------------------------------------

# ------------------------------ 
# palette 1

# ggplot palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=60 , c=100)[1:n]
}

# ------------------------------ 
# color blind palette

colors_muted <- c("#DC050C", "#E8601C", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6", "#F1932D", "#F6C141", "#F7EE55", "#4EB265", "#90C987", "#CAEDAB")
color_ramp <- c(colors_muted, gg_color_hue(max(1, nlevels(labels$label) - length(colors_muted))))

colors_tsne <- color_ramp[1:nlevels(labels$label)]
names(colors_tsne) <- levels(labels$label)


# ------------------------------------------------------------
# Load tSNE data
# ------------------------------------------------------------

load(path_rtsne_out)

rtsne_data <- read.table(path_rtsne_data, header = TRUE, sep = "\t", as.is = TRUE)

## keep the tSNE cells that were used in clustering
cells2keep_rtsne <- rtsne_data$cell_index %in% clustering[, "cell_id"]


# ------------------------------------------------------------
# tSNE plots
# ------------------------------------------------------------

# get clustering for cells that were used in tSNE
names(clust) <- clustering[, "cell_id"]
clust_tsne <- clust[as.character(rtsne_data$cell_index[cells2keep_rtsne])]

ggdf <- data.frame(tSNE1 = rtsne_out$Y[cells2keep_rtsne,1], tSNE2 = rtsne_out$Y[cells2keep_rtsne,2], cluster = clust_tsne, sample = rtsne_data$sample_name[cells2keep_rtsne])

# add group info
mm <- match(ggdf$sample, md$shortname)
ggdf$group <- md$condition[mm]


# use cluster labels instead of numbers
mm <- match(ggdf$cluster, labels$cluster)
ggdf$cluster <- labels$label[mm]
ggdf$cluster <- factor(ggdf$cluster, levels = levels(labels$label))


# skipp the "drop" cluster

ggdf <- ggdf[ggdf$cluster != "drop", ]
ggdf$cluster <- factor(ggdf$cluster)



# -----------------------------------
### Plot of tsne - all cells, all clusters

## facet per group
# ggp <- ggplot(ggdf,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
#   geom_point(size=1) +
#   facet_wrap(~ group) +
#   labs(x = "tSNE 1", y="tSNE 2")+ 
#   theme_bw() +
#   theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold"), legend.key = element_blank()) +
#   scale_color_manual(values = colors_tsne[levels(ggdf$cluster)]) +
#   guides(colour = guide_legend(override.aes = list(size = 5)))
# 
# pdf(file.path(outdir, paste0(prefix, "tSNEgroup", suffix, ".pdf")), width = pdf_width, height = pdf_height)  
# print(ggp)
# dev.off()


## one plot 
ggp <- ggplot(ggdf,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
  geom_point(size = 1) +
  labs(x = "tSNE 1", y="tSNE 2")+ 
  theme_bw() +
  theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold"), legend.key = element_blank()) +
  scale_color_manual(values = colors_tsne[levels(ggdf$cluster)]) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

pdf(file.path(outdir, paste0(prefix, "tSNEone", suffix, ".pdf")), width = 9, height = 7)                 
print(ggp)
dev.off()


## facet per sample
# ggp <- ggplot(ggdf,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
#   geom_point(size=1) +
#   facet_wrap(~ sample) +
#   labs(x = "tSNE 1", y="tSNE 2")+ 
#   theme_bw() +
#   theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold"), legend.key = element_blank()) +
#   scale_color_manual(values = colors_tsne[levels(ggdf$cluster)]) +
#   guides(colour = guide_legend(override.aes = list(size = 5)))
# 
# pdf(file.path(outdir, paste0(prefix, "tSNEsample", suffix, ".pdf")), width = pdf_width, height = pdf_height)          
# print(ggp)
# dev.off()



# -----------------------------------
### Plot of tsne - large clusters
if(any(grepl("tsne_cmin", args))){
  
  tc <- table(ggdf$cluster)
  tc
  
  pdf(file.path(outdir, paste0(prefix, "cmin", suffix, ".pdf")), width = 10, height = 5)        
  barplot(tc)
  abline(h = tsne_cmin, col = "blue")
  dev.off()
  
  kk <- names(tc[tc > tsne_cmin])
  
  ggdf_sub <- ggdf[ggdf$cluster %in% kk,]
  ggdf_sub$cluster <- factor(ggdf_sub$cluster)
  
  
  ## facet per group
  # ggp <- ggplot(ggdf_sub,  aes(x = tSNE1, y = tSNE2, color = cluster )) +
  #   geom_point(size=1) + 
  #   facet_wrap( ~ group)+ 
  #   labs(x = "tSNE 1", y="tSNE 2")+ 
  #   theme_bw() +
  #   theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold"), legend.key = element_blank()) +
  #   scale_color_manual(values = colors_tsne[levels(ggdf_sub$cluster)]) + 
  #   guides(colour = guide_legend(override.aes = list(size = 5)))
  # 
  # pdf(file.path(outdir, paste0(prefix, "tSNEgroup_subset", suffix, ".pdf")), width = pdf_width, height = pdf_height)             
  # print(ggp)
  # dev.off()
  
}




# -----------------------------------
### Plot of tsne - cells that are close to the centers of clusters (no outliers), all clusters
if(any(grepl("tsne_distse", args))){
  
  el <- rtsne_data[cells2keep_rtsne, -grep("cell_index|sample_name", colnames(rtsne_data))]
  a <- aggregate(el, by = list(clust_tsne), FUN=median)
  
  dists <- distse <- rep(NA, nrow(el))
  
  clust_lev <- a$Group.1
  
  for(i in 1:length(clust_lev)) {
    # i = 1
    cent <- a[i,-1,drop=FALSE]
    data <- el[clust_tsne == clust_lev[i], , drop=FALSE]
    d <- 1 - cosine( t(rbind(cent,data)) )[1,-1]
    dists[clust_tsne == clust_lev[i]] <- d
    cent <- matrix(rep( as.numeric(cent), nrow(data)), byrow=TRUE, ncol=ncol(cent))
    #distse[clust_tsne[s]==i] <- sqrt(rowSums((cent-data)^2))
    distse[clust_tsne == clust_lev[i]] <- rowMeans(abs(cent-data))
  }
  
  
  pdf(file.path(outdir, paste0(prefix, "distse", suffix, ".pdf")), width = 7, height = 7)        
  hist(distse, breaks = 100)
  abline(v = tsne_distse, col = "blue")
  dev.off()
  
  
  ggdf_sub <- ggdf[distse < tsne_distse, ]
  ggdf_sub$cluster <- factor(ggdf_sub$cluster)
  
  
  ## facet per group
  # ggp <- ggplot(ggdf_sub,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
  #   geom_point(size=1) + 
  #   facet_wrap( ~ group) +
  #   labs(x = "tSNE 1", y="tSNE 2")+ 
  #   theme_bw() +
  #   theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold"), legend.key = element_blank())+
  #   scale_color_manual(values = colors_tsne[levels(ggdf_sub$cluster)]) +
  #   guides(colour = guide_legend(override.aes = list(size = 5)))
  # 
  # 
  # pdf(file.path(outdir, paste0(prefix, "tSNEgroup_filtered", suffix, ".pdf")), width = pdf_width, height = pdf_height)
  # print(ggp)
  # dev.off()
  
  
  
  ## one plot
  ggp <- ggplot(ggdf_sub,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
    geom_point(size=1) + 
    labs(x = "tSNE 1", y="tSNE 2")+ 
    theme_bw() +
    theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold"), legend.key = element_blank())+
    scale_color_manual(values = colors_tsne[levels(ggdf_sub$cluster)]) +
    guides(colour = guide_legend(override.aes = list(size = 5)))
  
  
  pdf(file.path(outdir, paste0(prefix, "tSNEone_filtered", suffix, ".pdf")), width = 9, height = 7)                 
  print(ggp)
  dev.off()
  
}








sessionInfo()















################################
### 03_plottsne done!
################################