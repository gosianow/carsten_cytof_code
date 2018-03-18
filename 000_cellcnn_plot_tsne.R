

Sys.time()

# Load packages
library(flowCore)
library(gdata)
library(Rtsne)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(coop) # cosine
library(limma) 
library(tools)

##############################################################################
# Test arguments
##############################################################################

args <- NULL

prefix='23_01_pca1_merging6_data23_'
outdir='../PD1_CellCnn_Lukas_June2017/cellcnn_tsne'
path_metadata='../carsten_cytof/PD1_project/CK_metadata/metadata_23_01.xlsx'
path_rtsne_out='../carsten_cytof/PD1_project/CK_2016-06-23_01/040_tsnemaps/23_01_pca1_rtsne_out.rds'
path_rtsne_data='../carsten_cytof/PD1_project/CK_2016-06-23_01/040_tsnemaps/23_01_pca1_rtsne_data.xls'
path_clustering='../carsten_cytof/PD1_project/CK_2016-06-23_01/030_heatmaps/23_01_pca1_merging6_clustering.xls'
path_clustering_labels='../carsten_cytof/PD1_project/CK_2016-06-23_01/030_heatmaps/23_01_pca1_merging6_clustering_labels.xls'
dir_cellcnn_files='../PD1_CellCnn_Lukas_June2017/panel1_base_data23/selected_cells'
day='base'



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


if(!file.exists(outdir)) 
  dir.create(outdir, recursive = TRUE)


suffix <- ""
pdf_width=15
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
head(clustering)


## clustering labels
labels <- read.table(path_clustering_labels, header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))
labels

mm <- match(clustering$cluster, labels$cluster)
clustering$label <- labels$label[mm]

clustering <- clustering[order(clustering$cell_id), , drop = FALSE]


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
colors_tsne


# ------------------------------------------------------------
# Load cell cnn results
# ------------------------------------------------------------

### CellCNN was run only on the R and NR samples 
samps2keep <- md$day == day & md$response %in% c("R", "NR")
samps2keep

fcnn <- file.path(dir_cellcnn_files, paste0(tools::file_path_sans_ext(md$filename[samps2keep]), "_transf_selected_cells.csv"))
names(fcnn) <- md$shortname[samps2keep]
fcnn


cnn_filters <- lapply(fcnn, function(x){
  
  cnn_filt <- read.csv(x, header = TRUE, sep = ",")
  
  out <- cnn_filt[, grep("binary", colnames(cnn_filt)), drop = FALSE]
  
  colnames(out) <- gsub("_binary", "", colnames(out))
  
  return(out)
  
})

sapply(cnn_filters, nrow)




# ------------------------------------------------------------
# Load tSNE data
# ------------------------------------------------------------

rtsne_out <- readRDS(path_rtsne_out)

rtsne_data <- read.table(path_rtsne_data, header = TRUE, sep = "\t", as.is = TRUE)

head(rtsne_data)


# ------------------------------------------------------------
# Identify cells to plot and create data to plot
# ------------------------------------------------------------


cell_ids_split <- split(clustering$cell_id, factor(clustering$sample_id, levels = md$shortname))

sapply(cell_ids_split, length)



tsne_cells <- lapply(cell_ids_split, function(x){
  
  x %in% rtsne_data$cell_index
  
})

## Checks
sapply(tsne_cells, table)



### Create a data frame with t-SNE output

rtsne <- data.frame(tSNE1 = rtsne_out$Y[, 1], tSNE2 = rtsne_out$Y[, 2], sample = rtsne_data$sample_name, cell_index = rtsne_data$cell_index, stringsAsFactors = FALSE)

rtsne <- rtsne[order(rtsne$cell_index), , drop = FALSE]


### Add clustering results

mm <- match(rtsne$cell_index, clustering$cell_id)
rtsne$label <- clustering$label[mm]


### Add cell cnn results

rtsne_split <- split(rtsne, factor(rtsne_data$sample_name, levels = md$shortname))

ggdf <- lapply(names(cnn_filters), function(i){
  
  cnn_filters_sub <- cnn_filters[[i]][tsne_cells[[i]], , drop = FALSE]
  
  out <- data.frame(rtsne_split[[i]], cnn_filters_sub)
  
  return(out)
})


ggdf <- plyr::rbind.fill(ggdf)
head(ggdf)


# add group info
mm <- match(ggdf$sample, md$shortname)
ggdf$group <- md$response[mm]

### Drop the "drop" cluster

ggdf <- ggdf[ggdf$label != "drop", , drop = FALSE]
ggdf$label <- factor(ggdf$label)

# ----------------------------------------------------------------------
# tSNE plots
# ----------------------------------------------------------------------


colors_filter <- c("#b2beb5", "#9932cc")
names(colors_filter) <- c("0", "1")


name_filters <- colnames(ggdf)[grep("filter", colnames(ggdf))]
name_filters

for(i in 1:length(name_filters))
  ggdf[, name_filters[i]] <- factor(ggdf[, name_filters[i]], levels = c("0", "1"))


for(i in 1:length(name_filters)){
  # i = 1
  
  suffix <- paste0("_", name_filters[i])
  
  ## one plot 
  ggp <- ggplot(ggdf, aes_string(x = "tSNE1", y = "tSNE2", color = name_filters[i])) +
    geom_point(size = 2) +
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
    scale_color_manual(values = colors_filter) +
    guides(colour = guide_legend(override.aes = list(size = 5)))
  
  pdf(file.path(outdir, paste0(prefix, "tSNEone", suffix, "_type1.pdf")), width = 9, height = 7)                 
  print(ggp)
  dev.off()
  

  ## one plot + colored by cluster
  ggp <- ggplot() +
    geom_point(data = ggdf[ggdf[, name_filters[i]] == "0", , drop = FALSE], aes(x = tSNE1, y = tSNE2, color = label), shape = 16, alpha = 1, size = 2) +
    geom_point(data = ggdf[ggdf[, name_filters[i]] == "1", , drop = FALSE], aes(x = tSNE1, y = tSNE2, fill = label), shape = 21, show.legend = FALSE, size = 2) +
    labs(x = "t-SNE1", y="t-SNE2") + 
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
    scale_fill_manual(values = colors_tsne) +
    scale_color_manual(values = colors_tsne) +
    guides(color = guide_legend(override.aes = list(size = 5)))
  
  pdf(file.path(outdir, paste0(prefix, "tSNEone", suffix, "_type2.pdf")), width = 9, height = 7)                 
  print(ggp)
  dev.off()
  
  
  
  ## facet per group
  ggp <- ggplot(ggdf, aes_string(x = "tSNE1", y = "tSNE2", color = name_filters[i])) +
    geom_point(size = 2) +
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
    scale_color_manual(values = colors_filter) +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
  facet_wrap(~ group) 
  
  pdf(file.path(outdir, paste0(prefix, "tSNEgroup", suffix, "_type1.pdf")), width = pdf_width, height = pdf_height)
  print(ggp)
  dev.off()
  
  
  ## facet per group + colored by cluster
  ggp <- ggplot() +
    geom_point(data = ggdf[ggdf[, name_filters[i]] == "0", , drop = FALSE], aes(x = tSNE1, y = tSNE2, color = label), shape = 16, alpha = 1, size = 2) +
    geom_point(data = ggdf[ggdf[, name_filters[i]] == "1", , drop = FALSE], aes(x = tSNE1, y = tSNE2, fill = label), shape = 21, show.legend = FALSE, size = 2, stroke = 1.5) +
    labs(x = "t-SNE1", y="t-SNE2") + 
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
    scale_fill_manual(values = colors_tsne) +
    scale_color_manual(values = colors_tsne) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    facet_wrap(~ group)
  
  pdf(file.path(outdir, paste0(prefix, "tSNEgroup", suffix, "_type2.pdf")), width = pdf_width, height = pdf_height)                 
  print(ggp)
  dev.off()
  
  
}




















sessionInfo()











