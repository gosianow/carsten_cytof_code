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
# tsnep_prefix='pca1_mergingNEW2_'
# path_rtsne_out='pca1_cl20_rtsne_out.rda'
# path_rtsne_data='pca1_cl20_rtsne_data.xls'
# path_clustering='pca1_mergingNEW2_clustering.xls'
# path_clustering_labels='pca1_mergingNEW2_clustering_labels.xls'
# tsne_cmin=1000
# pdf_width=15
# pdf_height=10
# path_metadata

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


# ------------------------------------------------------------
# Load cluster data
# ------------------------------------------------------------

clust <- read.table(file.path(hmDir, path_clustering), header = TRUE, sep = "\t", as.is = TRUE)
clust <- clust[, 1]

labels <- read.table(file.path(hmDir, path_clustering_labels), header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))



### Colors for 20 clusters 
# ggplot palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=60 , c=100)[1:n]
}

color_ramp <- c(colorRampPalette(brewer.pal(12,"Paired"))(12)[-c(11)],  gg_color_hue(max(1, nlevels(labels$label)-11)) )

tsne_colors <- color_ramp[1:nlevels(labels$label)]
names(tsne_colors) <- levels(labels$label)

# ------------------------------------------------------------
# Load tSNE data
# ------------------------------------------------------------

load(file.path(sneDir, path_rtsne_out))

rtsne_data <- read.table(file.path(sneDir, path_rtsne_data), header = TRUE, sep = "\t", as.is = TRUE)



# ------------------------------------------------------------
# tSNE plots
# ------------------------------------------------------------

# get clustering for cells that were used in tSNE
clust <- clust[rtsne_data$cell_index]

ggdf <- data.frame(tSNE1 = rtsne_out$Y[,1], tSNE2 = rtsne_out$Y[,2], cluster = clust, sample = rtsne_data$sample_name)

# add group info
mm <- match(ggdf$sample, md$shortname)
ggdf$group <- md$condition[mm]

# use cluster labels instead of numbers
mm <- match(ggdf$cluster, labels$cluster)
ggdf$cluster <- labels$label[mm]
ggdf$cluster <- factor(ggdf$cluster, levels = levels(labels$label))


# skipp drop cluster

ggdf <- ggdf[ggdf$cluster != "drop", ]
ggdf$cluster <- factor(ggdf$cluster)




### Plot of tsne - all cells, all clusters

## facet per group
ggp <- ggplot(ggdf,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
  geom_point(size=1) +
  facet_wrap(~ group) +
  labs(x = "tSNE 1", y="tSNE 2")+ 
  theme_bw() +
  theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold")) +
  scale_color_manual(values = tsne_colors[levels(ggdf$cluster)]) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

pdf(file.path(sneDir, paste0(prefix, "tSNE.pdf")), width = pdf_width, height = pdf_height)                 
print(ggp)
dev.off()


## one plot 
ggp <- ggplot(ggdf,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
  geom_point(size = 1) +
  labs(x = "tSNE 1", y="tSNE 2")+ 
  theme_bw() +
  theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold")) +
  scale_color_manual(values = tsne_colors[levels(ggdf$cluster)]) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

pdf(file.path(sneDir, paste0(prefix, "tSNEone.pdf")), width = 9, height = 7)                 
print(ggp)
dev.off()



### Plot of tsne - large clusters
tc <- table(ggdf$cluster)
kk <- names(tc[tc > tsne_cmin])

ggdf_sub <- ggdf[ggdf$cluster %in% kk,]
ggdf_sub$cluster <- factor(ggdf_sub$cluster)

## facet per group
ggp <- ggplot(ggdf_sub,  aes(x = tSNE1, y = tSNE2, color = cluster )) +
  geom_point(size=1) + 
  facet_wrap( ~ group)+ 
  labs(x = "tSNE 1", y="tSNE 2")+ 
  theme_bw() +
  theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold")) +
  scale_color_manual(values = tsne_colors[levels(ggdf_sub$cluster)]) + 
  guides(colour = guide_legend(override.aes = list(size = 5)))

pdf(file.path(sneDir, paste0(prefix, "tSNE_subset.pdf")), width = pdf_width, height = pdf_height)                 
print(ggp)
dev.off()





### Plot of tsne - cells that are close to the centers of clusters (no outliers), all clusters
el <- rtsne_data[, -c(1:2)]
a <- aggregate(el, by = list(clust), FUN=median)

dists <- distse <- rep(NA, nrow(el))

clust_lev <- a$Group.1

for(i in 1:length(clust_lev)) {
  # i = 1
  cent <- a[i,-1,drop=FALSE]
  data <- el[clust == clust_lev[i], , drop=FALSE]
  d <- 1 - cosine( t(rbind(cent,data)) )[1,-1]
  dists[clust == clust_lev[i]] <- d
  cent <- matrix(rep( as.numeric(cent), nrow(data)), byrow=TRUE, ncol=ncol(cent))
  #distse[clust[s]==i] <- sqrt(rowSums((cent-data)^2))
  distse[clust == clust_lev[i]] <- rowMeans(abs(cent-data))
}


ggdf_sub <- ggdf[distse < 0.25, ]
ggdf_sub$cluster <- factor(ggdf_sub$cluster)

## facet per group
ggp <- ggplot(ggdf_sub,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
  geom_point(size=1) + 
  facet_wrap( ~ group) +
  labs(x = "tSNE 1", y="tSNE 2")+ 
  theme_bw() +
  theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold"))+
  scale_color_manual(values = tsne_colors[levels(ggdf_sub$cluster)]) +
  guides(colour = guide_legend(override.aes = list(size = 5)))


pdf(file.path(sneDir, paste0(prefix, "tSNE_filtered.pdf")), width = pdf_width, height = pdf_height)                 
print(ggp)
dev.off()



## one plot
ggp <- ggplot(ggdf_sub,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
  geom_point(size=1) + 
  labs(x = "tSNE 1", y="tSNE 2")+ 
  theme_bw() +
  theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold"))+
  scale_color_manual(values = tsne_colors[levels(ggdf_sub$cluster)]) +
  guides(colour = guide_legend(override.aes = list(size = 5)))


pdf(file.path(sneDir, paste0(prefix, "tSNEone_filtered.pdf")), width = 9, height = 7)                 
print(ggp)
dev.off()








sessionInfo()















################################
### 03_plottsne done!
################################