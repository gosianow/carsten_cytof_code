##############################################################################
## <<04_frequencies.R>>

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
# freq_prefix='pca1_cl20_'
# path_clustering='pca1_cl20_clustering.xls'
# path_clustering_labels='pca1_cl20_clustering_labels.xls'


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

prefix <- freq_prefix


# ------------------------------------------------------------
# define directories
# ------------------------------------------------------------

fcsDir <- "010_cleanfcs"; if( !file.exists(fcsDir) ) dir.create(fcsDir)
pcaDir <- "020_pcascores"; if( !file.exists(pcaDir) ) dir.create(pcaDir)
hmDir <- "030_heatmaps"; if( !file.exists(hmDir) ) dir.create(hmDir)
frqDir <- "050_frequencies"; if( !file.exists(frqDir) ) dir.create(frqDir)

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

clust <- read.table(file.path(hmDir, path_clustering), header = TRUE, sep = "\t", as.is = TRUE)
clust <- clust[, 1]

labels <- read.table(file.path(hmDir, path_clustering_labels), header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))


# ---------------------------------------
# plots and analysis of cluster frequencies
# ---------------------------------------


# calculate frequencies
samp <- rep( names(fcsT), sapply(fcsT, nrow) )
freq <- table( cluster = clust, samp )

# use labels as names of clusters
mlab <- match(rownames(freq), labels$cluster)
rownames(freq) <- labels$label[mlab]

prop <- t(t(freq) / colSums(freq))



### Save the frequencies and proportions
prop_out <- data.frame(cluster = rownames(prop), as.data.frame.matrix(prop * 100))
freq_out <- data.frame(cluster = rownames(freq), as.data.frame.matrix(freq))

write.table(prop_out, file=file.path(frqDir, paste0(prefix, "frequencies.xls")), row.names=FALSE, quote=FALSE, sep="\t")
write.table(freq_out, file=file.path(frqDir,paste0(prefix, "counts.xls")), row.names=FALSE, quote=FALSE, sep="\t")





### Data frame preparation for ggplot

# ggdf <- melt(prop, value.name = "prop")
# 
# # use labels as clusters
# ggdf$cluster <- factor(ggdf$cluster, levels = labels$cluster, labels = labels$label)
# 
# # add group info
# mm <- match(ggdf$samp, md$shortname)
# ggdf$group <- md$condition[mm]
# 
# ggds <- ddply(ggdf, .(group, cluster), summarise, mean = mean(prop), sd = sd(prop))
# 
# ggp <- list()
# 
# clusters <- levels(ggdf$cluster)
# 
# for(i in 1:nlevels(ggdf$cluster)){
#   
#   df <- ggdf[ggdf$cluster == clusters[i], , drop = FALSE]
#   ds <- ggds[ggds$cluster == clusters[i], , drop = FALSE]
#   
#   ggp[[i]] <- ggplot(df, aes(x = group, y = prop)) + 
#     geom_point(size=3, shape = 1, stroke = 1) + 
#     geom_point(data=ds, aes(x=group, y=mean), colour='blue', size=2) +
#     geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean-sd, ymax=mean+sd), colour='blue', width=0.05) +
#     ggtitle(clusters[i]) + 
#     theme_bw() +
#     ylab("proportion") + 
#     xlab("") +
#     theme(axis.text.x  = element_text(size=13, face="bold"),
#       axis.title.y  = element_text(size=13, face="bold"))
#   
# }
# 
# 
# pdf(file.path(frqDir, paste0(prefix, "frequencies.pdf")), w=7, h=5, onefile=TRUE)
# for (i in seq(length(ggp)))
#   print(ggp[[i]])
# dev.off()






sessionInfo()













################################
### Done!
################################