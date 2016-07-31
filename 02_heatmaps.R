##############################################################################
## <<02_heatmaps.R>>

# BioC 3.3
# Created 27 July 2016
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
library(pheatmap)
library(RColorBrewer)
library(ggdendro)

##############################################################################
# Test arguments
##############################################################################

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
# heatmap_prefix='pca1_cl20_'
# path_panel='panel.xlsx'
# path_clustering_observables='pca1_clustering_observables.xls'
# path_clustering='pca1_cl20_clustering.xls'
# path_clustering_labels='pca1_cl20_clustering_labels.xls'
# path_pcascore='princompscore_by_sample.xls'
# xspace=2

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

prefix <- heatmap_prefix


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

# define FCS file names
f <- file.path(fcsDir, md$filename)
names(f) <- md$shortname

# read raw FCS files in
fcs <- lapply(f, read.FCS)


# read panel, pick which columns to use
panel <- read.xls(path_panel,stringsAsFactors=FALSE)

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


if(!grepl("/", path_pcascore)){
  prs <- read.table(file.path(pcaDir, path_pcascore), header = TRUE, sep = "\t", as.is = TRUE)
}else{
  prs <- read.table(path_pcascore, header = TRUE, sep = "\t", as.is = TRUE)
}

rownames(prs) <- prs$mass


if(!grepl("/", path_clustering_observables)){
  clust_observ <- read.table(file.path(hmDir, path_clustering_observables), header = TRUE, sep = "\t", as.is = TRUE)
}else{
  clust_observ <- read.table(path_clustering_observables, header = TRUE, sep = "\t", as.is = TRUE)
}

clust_observ <- clust_observ[, 1]


clust <- read.table(file.path(hmDir, path_clustering), header = TRUE, sep = "\t", as.is = TRUE)
clust <- clust[, 1]

labels <- read.table(file.path(hmDir, path_clustering_labels), header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))

# ------------------------------------------------------------

### Indeces of observables used for clustering 

scols <- which(fcs_panel$colnames %in% clust_observ)

# ordered by decreasing pca score
scols <- scols[order(prs[fcs_panel$colnames[scols], "avg_score"], decreasing = TRUE)]


### Indeces of other observables

xcols <- setdiff(cols, scols)

# ordered by decreasing pca score
xcols <- xcols[order(prs[fcs_panel$colnames[xcols], "avg_score"], decreasing = TRUE)]




# ------------------------------------------------------------
# Prepare data
# ------------------------------------------------------------



### Calculate cluster frequencies

samp <- rep( names(fcsT), sapply(fcsT,nrow) )
freq <- table( cluster = clust, samp )
# Use labels as new cluster names
mlab <- match(rownames(freq), labels$cluster)
rownames(freq) <- labels$label[mlab]

prop <- t(t(freq) / colSums(freq))
# normalize 
rp <- sweep(prop, 1, STATS=rowMeans(prop), FUN="-")




### Get marker expression

# re-extract the data as list of tables
es <- lapply(fcsT, function(u) {
  exprs(u)[,scols]
})

esX <- lapply(fcsT, function(u) {
  exprs(u)[,xcols]
})

e <- do.call("rbind",es)
eX <- do.call("rbind",esX)

colnames(e) <- fcs_panel$Antigen[scols]
colnames(eX) <- fcs_panel$Antigen[xcols]


a <- aggregate( e, by=list(clust), FUN=median)
aX <- aggregate( eX, by=list(clust), FUN=median)

mlab <- match(a$Group.1, labels$cluster)
rownames(a) <- rownames(aX) <- labels$label[mlab]



### Normalize the data to 01

# normalize to 0-1
el <- e
rng <- apply(el,2,quantile,p=c(.01,.99))
for(i in 1:ncol(el)) {
  el[,i] <- (el[,i]-rng[1,i])/(rng[2,i]-rng[1,i])
}
el[el<0] <- 0
el[el>1] <- 1


# normalize to 0-1
elX <- eX
rng <- apply(elX,2,quantile,p=c(.01,.99))
for(i in 1:ncol(elX)) {
  elX[,i] <- (elX[,i]-rng[1,i])/(rng[2,i]-rng[1,i])
}
elX[elX<0] <- 0
elX[elX>1] <- 1


al <- aggregate( el, by=list(clust), FUN=median)
alX <- aggregate( elX, by=list(clust), FUN=median)

mlab <- match(al$Group.1, labels$cluster)
rownames(al) <- rownames(alX) <- labels$label[mlab]

colnames(al)[1] <- "cluster"
colnames(alX)[1] <- "cluster"


### Save the median expression

write.table(al, file.path(hmDir, paste0(prefix, "cluster_medians_in.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(alX, file.path(hmDir, paste0(prefix, "cluster_medians_ex.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# ------------------------------------------------------------
# Multiheatmaps
# ------------------------------------------------------------


### Prepare colors and breaks for the heatmaps
cL <- NULL
ncol <- 20

palette <- colorRampPalette(c("green", "black", "red"),space = "Lab")
br <- seq(-20,20,length=ncol)
col <- palette(ncol-1)
cL[[1]] <- list(breaks=br,colors=col)

palette <- colorRampPalette(c("blue", "orange", "red"),space = "Lab")
col <- palette(ncol-1)
br <- seq(0,40,length=ncol)
cL[[2]] <- list(breaks=br,colors=col)

palette <- colorRampPalette(c("darkblue", "yellow", "red"),space = "Lab")
col <- palette(ncol-1)
br <- seq(0,1,length=ncol)
cL[[4]] <- cL[[3]] <- list(breaks=br,colors=col)


testD <- list(rp*100, prop*100, al[,-1], alX[,-1])

testD[[1]][testD[[1]]>max(cL[[1]]$breaks)] <- max(cL[[1]]$breaks)
testD[[1]][testD[[1]]<min(cL[[1]]$breaks)] <- min(cL[[1]]$breaks)
testD[[2]][testD[[2]]>max(cL[[2]]$breaks)] <- max(cL[[2]]$breaks)
testD[[2]][testD[[2]]<min(cL[[2]]$breaks)] <- min(cL[[2]]$breaks)


### Heatmaps will all the clusters
pdf(file.path(hmDir, paste0(prefix, "multiheatmap.pdf")), width = 10, height = 6)
multiHeatmap(testD,cL, xspace = xspace, ystarts=c(.25,.9,.925,.95,.98), clabelcex=.5)
dev.off()


### Heatmaps with larger clusters
w <- which(rowSums(prop>=.05)>=5)
testDs <- lapply(testD, function(u) u[w,,drop=FALSE])

pdf(file.path(hmDir, paste0(prefix, "multiheatmap_subset.pdf")), width = 10, height = 6)
multiHeatmap(testDs,cL, xspace = xspace, ystarts=c(.25,.9,.925,.95,.98), clabelcex=.5)
dev.off()


# ------------------------------------------------------------
# Plot expression of markers in each cluster
# ------------------------------------------------------------

df <- data.frame(e, clust = factor(clust), samp)
dfm <- melt(df, id.vars = c("clust","samp"))

dfm$clust <- factor(dfm$clust, labels = labels$label)


ggp <- ggplot(dfm, aes(x=value)) + 
  geom_density(size=1,adjust=5) + 
  facet_grid(clust~variable, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf(file.path(hmDir, paste0(prefix, "distros.pdf")), w = ncol(el), h = nrow(labels))
print(ggp)
dev.off()



# ------------------------------------------------------------
# pheatmaps
# ------------------------------------------------------------


expr <- cbind(al[,-1], alX[,-1])

labels_row <- labels$label
labels_col <- colnames(expr)


cluster_rows <- hclust(dist(al[,-1]), method = "average")


pheatmap(mat = expr, color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(100), cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, breaks = seq(from = 0, to = 1, length.out = 101), legend_breaks = seq(from = 0, to = 1, by = 0.2), gaps_col = length(scols), fontsize_row = 10, fontsize_col = 10, fontsize = 6, filename = file.path(hmDir, paste0(prefix, "pheatmap_row_clust.pdf")), width = 10, height = 7)


pheatmap(mat = expr, color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(100), cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, breaks = seq(from = 0, to = 1, length.out = 101), legend_breaks = seq(from = 0, to = 1, by = 0.2), gaps_col = length(scols), fontsize_row = 10, fontsize_col = 10, fontsize = 6, filename = file.path(hmDir, paste0(prefix, "pheatmap.pdf")), width = 10, height = 7)




# ------------------------------------------------------------
# ggplot tile
# ------------------------------------------------------------

ggdf <- data.frame(cluster = labels$label, cbind(al[,-1], alX[,-1]))

ggdf$cluster <- factor(ggdf$cluster, levels = rev(levels(labels$label)))

ggdfm <- melt(ggdf, id.var = "cluster", variable.name = "observable")
ggdfm$observable <- factor(ggdfm$observable, levels = colnames(ggdf)[-1])


ggp <- ggplot(ggdfm, aes(x = observable, y = cluster, fill = value)) + 
  geom_tile() + 
  geom_text(aes(label = round(value, 2)), color = "black", size = 2.5) + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  xlab("") + 
  ylab("") + 
  theme_bw() +
  scale_fill_gradientn("", colours = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(100), limits=c(0, 1)) +
  theme(panel.background = element_rect(fill = NA, colour = NA), axis.ticks = element_blank(), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size = 14)) 


pdf(file.path(hmDir, paste0(prefix, "ggtile.pdf")), width = 10, height = 7)
print(ggp)
dev.off()



### Row order from clustering

cluster_rows <- hclust(dist(al[,-1]), method = "average")

ggdfm$cluster <- factor(ggdfm$cluster, levels = rev(levels(labels$label)[cluster_rows$order]))

ggp <- ggplot(ggdfm, aes(x = observable, y = cluster, fill = value)) + 
  geom_tile() + 
  geom_text(aes(label = round(value, 2)), color = "black", size = 2.5) + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  xlab("") + 
  ylab("") + 
  theme_bw() +
  scale_fill_gradientn("", colours = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(100), limits=c(0, 1)) +
  theme(panel.background = element_rect(fill = NA, colour = NA), axis.ticks = element_blank(), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size = 14)) 


pdf(file.path(hmDir, paste0(prefix, "ggtile_row_clust.pdf")), width = 10, height = 7)
print(ggp)
dev.off()



# ------------------------------------------------------------
# ggplot tile + ggdendogram - not finished yet
# ------------------------------------------------------------


# ggp <- ggdendrogram(cluster_rows, rotate = FALSE, size = 2) 
# 
# 
# pdf(file.path(hmDir, paste0(prefix, "ggdendro_row_clust.pdf")), width = 7, height = 4)
# print(ggp)
# dev.off()




# dhc <- as.dendrogram(cluster_rows)
# ddata <- dendro_data(dhc, type = "rectangle")
# 
# 
# ggp <- ggplot(segment(ddata)) + 
#   geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
#   geom_text(data = label(ddata), aes(x=x, y=y, label=label, hjust=0), size=3) +
#   coord_flip() + 
#   scale_y_reverse(expand=c(0.2, 0)) +
#   theme_dendro()
# 
# 
# 
# pdf(file.path(hmDir, paste0(prefix, "ggdendro_row_clust.pdf")), width = 2, height = 7)
# print(ggp)
# dev.off()
# 
# 
# 
# geom_segment(data=segment(hcdata), aes(x=x, y=y, xend=xend, yend=yend)) +
#   geom_text(data=label(hcdata), aes(x=x, y=y, label=label, hjust=0), size=3) +
#   coord_flip() + 
#   scale_y_reverse(expand=c(0.2, 0))
# 











sessionInfo()













################################
### 02_heatmaps done!
################################