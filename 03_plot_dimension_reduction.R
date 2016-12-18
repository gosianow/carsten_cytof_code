##############################################################################
## <<03_plot_dimension_reduction.R.R>>

# BioC 3.3
# Created 2 Nov 2016
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

##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
dr_prefix='23_01_pca1_cl20_raw_'
dr_outdir='040_dimension_reduction'
path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.xlsx'
path_dr_data='040_dimension_reduction/23_01_pca1_raw_pca_data.xls'
path_clustering='030_heatmaps/23_01_pca1_cl20_clustering.xls'
path_clustering_labels='030_heatmaps/23_01_pca1_cl20_clustering_labels.xls'
dr_suffix='_pca'
pdf_width=15
pdf_height=10


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

prefix <- dr_prefix
suffix <- dr_suffix
outdir <- dr_outdir

if(!file.exists(outdir)) 
  dir.create(outdir)


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
# Load dimension reduction data
# ------------------------------------------------------------


dr_data <- read.table(path_dr_data, header = TRUE, sep = "\t", as.is = TRUE)

## keep cells that were used in clustering
cells2keep <- dr_data$cell_index %in% clustering[, "cell_id"]


# ------------------------------------------------------------
# Dimension reduction plots
# ------------------------------------------------------------

# get clustering for cells that were used in tSNE
names(clust) <- clustering[, "cell_id"]
clust_sub <- clust[as.character(dr_data$cell_index[cells2keep])]

ggdf <- data.frame(dim1 = dr_data[cells2keep, "dim1"], dim2 = dr_data[cells2keep, "dim2"], cluster = clust_sub, sample = dr_data$sample_name[cells2keep])

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
### Plot - all cells, all clusters


## one plot 
ggp <- ggplot(ggdf,  aes(x = dim1, y = dim2, color = cluster)) +
  geom_point(size = 1) +
  theme_bw() +
  theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold"), legend.key = element_blank()) +
  scale_color_manual(values = colors_tsne[levels(ggdf$cluster)]) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

pdf(file.path(outdir, paste0(prefix, "dim_red", suffix, ".pdf")), width = 9, height = 7)                 
print(ggp)
dev.off()









sessionInfo()















################################
### 03_plot_dimension_reduction.R done!
################################