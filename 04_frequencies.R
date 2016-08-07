##############################################################################
## <<04_frequencies.R>>

# BioC 3.3
# Created 27 July 2016


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

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_02_CD4_merging2'
# freq_prefix='pnlCD4_pca1_merging_CD4_'
# path_clustering='pnlCD4_pca1_merging_CD4_clustering.xls'
# path_clustering_labels='pnlCD4_pca1_merging_CD4_clustering_labels.xls'


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

# define FCS file names
f <- file.path(fcsDir, md$filename)
names(f) <- md$shortname

# read raw FCS files in
fcs <- lapply(f, read.FCS)


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
samp <- rep( names(fcs), sapply(fcs, nrow) )
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





### Plot frequencies

ggdf <- melt(prop, value.name = "prop")

# use labels as clusters
ggdf$cluster <- factor(ggdf$cluster, levels = labels$label)

# add group info
mm <- match(ggdf$samp, md$shortname)
ggdf$group <- factor(md$condition[mm])

## merge base_HD and tx_HD into one level - HD
new_levels <- levels(ggdf$group)
new_levels[grep("HD", new_levels)] <- "HD"
levels(ggdf$group) <- new_levels

## replace _ with \n
levels(ggdf$group) <- gsub("_", "\n", levels(ggdf$group))

## calculate mean and sd for the error bars on the plot
ggds <- ddply(ggdf, .(group, cluster), summarise, mean = mean(prop), sd = sd(prop))



ggp <- list()

clusters <- levels(ggdf$cluster)

for(i in 1:nlevels(ggdf$cluster)){

  df <- ggdf[ggdf$cluster == clusters[i], , drop = FALSE]
  ds <- ggds[ggds$cluster == clusters[i], , drop = FALSE]

  ggp[[i]] <- ggplot(df, aes(x = group, y = prop)) +
    geom_jitter(size=2.5, shape = 17, width = 0.5) +
    geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean, ymax=mean), colour='black', width=0.4) +
    geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean-sd, ymax=mean+sd), colour='black', width=0.25) +
    ggtitle(clusters[i]) +
    theme_bw() +
    ylab("Frequency") +
    xlab("") +
    theme(axis.text.x = element_text(size=12, face="bold"), 
      axis.title.y = element_text(size=12, face="bold"), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(), 
      axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"), 
      axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))

}


pdf(file.path(frqDir, paste0(prefix, "frequencies.pdf")), w=5, h=4, onefile=TRUE)
for(i in seq(length(ggp)))
  print(ggp[[i]])
dev.off()


## plot using facet wrap

# ggp <- ggplot(ggdf, aes(x = group, y = prop)) +
#   geom_jitter(size=2.5, shape = 17, width = 0.5) +
#   # geom_point(data = ggds, aes(x=group, y=mean), colour='black', size=2) +
#   geom_errorbar(data=ggds, aes(x=group, y=mean, ymin=mean, ymax=mean), colour='black', width=0.4) +
#   geom_errorbar(data=ggds, aes(x=group, y=mean, ymin=mean-sd, ymax=mean+sd), colour='black', width=0.25) +
#   facet_wrap(~ cluster, ncol = 1, scales = "free") +
#   theme_bw() +
#   ylab("Frequency") +
#   xlab("") +
#   theme(axis.text.x = element_text(size=12, face="bold"), 
#     axis.title.y = element_text(size=12, face="bold"), 
#     panel.grid.major = element_blank(), 
#     panel.grid.minor = element_blank(), 
#     panel.border = element_blank(), 
#     axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"), 
#     axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"), 
#     strip.text = element_text(size = 14, face="bold"), 
#     strip.background = element_rect(colour = "white",  fill="white"),
#     panel.margin.y = unit(2, "lines"))
# 
# 
# pdf(file.path(frqDir, paste0(prefix, "frequencies_facet.pdf")), w=5, h = 4 * nlevels(ggdf$cluster), onefile=TRUE)
#   print(ggp)
# dev.off()












sessionInfo()













################################
### 04_frequencies done!
################################