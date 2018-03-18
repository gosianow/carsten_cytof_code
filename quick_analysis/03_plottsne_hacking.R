

Sys.time()

# Load packages
library(gdata)
library(Rtsne)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(coop) # cosine
library(limma) 


##############################################################################
# Test arguments
##############################################################################

args <- NULL

# prefix='23all_03v2_pca1_merging4_'
# outdir='../carsten_cytof/PD1_project/CK_2016-06-23_03all/040_tsnemaps_hacking'
# path_metadata='../carsten_cytof/PD1_project/CK_metadata/metadata_23_03all.xlsx'
# path_rtsne_out='../carsten_cytof/PD1_project/CK_2016-06-23_03all/040_tsnemaps/23all_03v2_pca1_rtsne_out.rds'
# path_rtsne_data='../carsten_cytof/PD1_project/CK_2016-06-23_03all/040_tsnemaps/23all_03v2_pca1_rtsne_data.xls'
# path_clustering='../carsten_cytof/PD1_project/CK_2016-06-23_03all/030_heatmaps/23all_03v2_pca1_merging5_clustering.xls'
# path_clustering_labels='../carsten_cytof/PD1_project/CK_2016-06-23_03all/030_heatmaps/23all_03v2_pca1_merging5_clustering_labels.xls'
# path_old_clustering='../carsten_cytof/PD1_project_back_2018_03_08/CK_2016-06-23_03all/030_heatmaps/23all_03v2_pca1_merging4_clustering.xls'
# path_old_clustering_labels='../carsten_cytof/PD1_project_back_2018_03_08/CK_2016-06-23_03all/030_heatmaps/23all_03v2_pca1_merging4_clustering_labels.xls'


prefix='23_03_pca1_merging4_'
outdir='../carsten_cytof/PD1_project/CK_2016-06-23_03/040_tsnemaps_hacking'
path_metadata='../carsten_cytof/PD1_project/CK_metadata/metadata_23_03.xlsx'
path_rtsne_out='../carsten_cytof/PD1_project/CK_2016-06-23_03/040_tsnemaps/23_03_pca1_rtsne_out.rds'
path_rtsne_data='../carsten_cytof/PD1_project/CK_2016-06-23_03/040_tsnemaps/23_03_pca1_rtsne_data.xls'
path_clustering='../carsten_cytof/PD1_project/CK_2016-06-23_03/030_heatmaps/23_03_pca1_merging5_clustering.xls'
path_clustering_labels='../carsten_cytof/PD1_project/CK_2016-06-23_03/030_heatmaps/23_03_pca1_merging5_clustering_labels.xls'
path_old_clustering='../carsten_cytof/PD1_project_back_2018_03_08/CK_2016-06-23_03/030_heatmaps/23_03_pca1_merging4_clustering.xls'
path_old_clustering_labels='../carsten_cytof/PD1_project_back_2018_03_08/CK_2016-06-23_03/030_heatmaps/23_03_pca1_merging4_clustering_labels.xls'


##############################################################################
# Read in the arguments
##############################################################################

# rm(list = ls())
# 
# args <- (commandArgs(trailingOnly = TRUE))
# for (i in 1:length(args)) {
#   eval(parse(text = args[[i]]))
# }
# 
# cat(paste0(args, collapse = "\n"), fill = TRUE)

##############################################################################

if(!file.exists(outdir)) 
  dir.create(outdir, recursive = TRUE)


suffix <- ""
pdf_width=22
pdf_height=7


# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

md <- read.xls(path_metadata, stringsAsFactors=FALSE)

# add more info about samples
cond_split <- strsplit2(md$condition, "_")
colnames(cond_split) <- c("day", "response")

md[, c("day", "response")] <- cond_split
md$response <- factor(md$response, levels = c("NR", "R", "HD"))
md$response <- factor(md$response)
md$day <- factor(md$day, levels = c("base", "tx"))
md$day <- factor(md$day)
md$patient_id <- factor(md$patient_id)


# ------------------------------------------------------------
# Load cluster data
# ------------------------------------------------------------

## clustering
clustering <- read.table(path_clustering, header = TRUE, sep = "\t", as.is = TRUE)

## clustering labels
labels <- read.table(path_clustering_labels, header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))
labels

mm <- match(clustering$cluster, labels$cluster)
clustering$label <- labels$label[mm]

clustering <- clustering[order(clustering$cell_id), , drop = FALSE]


# ------------------------------------------------------------
# Load old cluster data
# ------------------------------------------------------------

## clustering
clustering_old <- read.table(path_old_clustering, header = TRUE, sep = "\t", as.is = TRUE)

## clustering labels
labels_old <- read.table(path_old_clustering_labels, header = TRUE, sep = "\t", as.is = TRUE)
labels_old <- labels_old[order(labels_old$cluster, decreasing = FALSE), ]
labels_old$label <- factor(labels_old$label, levels = unique(labels_old$label))
labels_old

mm <- match(clustering_old$cluster, labels_old$cluster)
clustering_old$label <- labels_old$label[mm]

clustering_old <- clustering_old[order(clustering_old$cell_id), , drop = FALSE]


table(clustering$label[clustering$sample_id == "base_HD1"], clustering_old$label[clustering_old$sample_id == "base_HD1"])


# ------------------------------------------------------------
# Colors for tSNE maps
# ------------------------------------------------------------


# ggplot palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=60 , c=100)[1:n]
}

# color blind palette

# colors_muted <- c("#DC050C", "#E8601C", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6", "#F1932D", "#F6C141", "#F7EE55", "#4EB265", "#90C987", "#CAEDAB")
colors_muted <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6", "#4EB265", "#CAEDAB", "#E7298A", "#E78AC3", "#666666", "#999999", "#FF7F00", "#FDB462")
color_ramp <- c(colors_muted, gg_color_hue(max(1, nlevels(labels$label) - length(colors_muted))))

colors_tsne <- color_ramp[1:nlevels(labels$label)]
names(colors_tsne) <- levels(labels$label)


# ------------------------------------------------------------
# Load tSNE data
# ------------------------------------------------------------

rtsne_out <- readRDS(path_rtsne_out)

rtsne_data <- read.table(path_rtsne_data, header = TRUE, sep = "\t", as.is = TRUE)

table(rtsne_data$sample_name)

# ------------------------------------------------------------
# Match cells between new t-SNE map and old clustering 
# In old clustering HD2 fcs is used instead of HD3 fcs, so we should skipp clusters for base_HD3 and tx_HD3
# ------------------------------------------------------------

# --------------
### Clustering 
# --------------

### Add new cell IDs that for each sample start at 1

clustering_cell_id_split <- split(clustering$cell_id, clustering$sample_id)
sapply(clustering_cell_id_split, range)

clustering_cell_id_new <- unlist(lapply(clustering_cell_id_split, function(x){ 1:length(x) }))

clustering$cell_id_new <- clustering_cell_id_new
clustering$cell_id_new2 <- paste0(clustering$sample_id, "_", clustering_cell_id_new)


## Check
clustering_cell_id_new_split <- split(clustering$cell_id_new, clustering$sample_id)
sapply(clustering_cell_id_new_split, range)


table(clustering$cell_id - clustering$cell_id_new)


# --------------
### Old clustering 
# --------------

### Add new cell IDs that for each sample start at 1

clustering_old_cell_id_split <- split(clustering_old$cell_id, clustering_old$sample_id)
sapply(clustering_old_cell_id_split, range)

clustering_old_cell_id_new <- unlist(lapply(clustering_old_cell_id_split, function(x){ 1:length(x) }))

clustering_old$cell_id_new <- clustering_old_cell_id_new
clustering_old$cell_id_new2 <- paste0(clustering_old$sample_id, "_", clustering_old_cell_id_new)


## Check
clustering_old_cell_id_new_split <- split(clustering_old$cell_id_new, clustering_old$sample_id)
sapply(clustering_old_cell_id_new_split, range)


### Skipp base_HD3 and tx_HD3 samples

clustering_old_sub <- clustering_old[! clustering_old$sample_id %in% c("base_HD3", "tx_HD3"), ]


# ------------------------------------------------------------
# Prepare t-SNE data
# ------------------------------------------------------------


ggdf <- data.frame(tSNE1 = rtsne_out$Y[, 1], tSNE2 = rtsne_out$Y[, 2], cell_index = rtsne_data$cell_index, sample = rtsne_data$sample_name, stringsAsFactors = FALSE)

table(ggdf$sample)


### Match cell ids from tsne with cell ids in clustering and add new cell ids (cell_id_new2)
mm <- match(ggdf$cell_index, clustering$cell_id)
ggdf$cell_id_new2 <- clustering$cell_id_new2[mm]
ggdf$cell_id_new <- clustering$cell_id_new[mm]

table(is.na(mm))

table(ggdf$cell_index - ggdf$cell_id_new)


### Add info about OLD clustering using cell_id_new2
mm <- match(ggdf$cell_id_new2, clustering_old_sub$cell_id_new2)
ggdf$cluster <- as.character(clustering_old_sub$label)[mm]

table(is.na(mm))


### Use grey color for cells from HD3 samples
table(ggdf$cluster, useNA = "always")

ggdf$cluster[is.na(ggdf$cluster)] <- "HD3_cells"

ggdf$cluster <- factor(ggdf$cluster, levels = c(levels(labels_old$label), "HD3_cells"))


### Update colors
colors_tsne_HD3_cells <- c(HD3_cells = "grey")
colors_tsne <- c(colors_tsne, colors_tsne_HD3_cells)


# skipp the "drop" cluster

skipp_drop <- ggdf$cluster != "drop"
ggdf <- ggdf[skipp_drop, ]
ggdf$cluster <- factor(ggdf$cluster)



# ------------------------------------------------------------
# tSNE plots
# ------------------------------------------------------------


### Plot of tsne - all cells, all clusters

## one plot 
ggp <- ggplot(ggdf,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
  geom_point(size = 1) +
  labs(x = "t-SNE1", y="t-SNE2")+ 
  theme_bw() +
  theme(axis.text = element_text(size = 12), 
    axis.title  = element_text(size = 15),
    strip.text = element_text(size = 15, hjust = 0),
    strip.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(size = 0.3, linetype = "solid", color = "black")) +
  scale_color_manual(values = colors_tsne[levels(ggdf$cluster)]) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

pdf(file.path(outdir, paste0(prefix, "tSNEone", suffix, ".pdf")), width = 9, height = 7)                 
print(ggp)
dev.off()




























