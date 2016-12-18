##############################################################################
## <<02_som_codes.R>>

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
library(cytofkit)
library(igraph)

##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
codes_prefix='23_01_pca1_cl20_'
codes_outdir='030_heatmaps'
path_codes='030_heatmaps/23_01_pca1_cl20_codes.xls'
pdf_width=15
pdf_height=10


rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
codes_prefix='23_01_pca1_mergingNEW2_'
codes_outdir='030_heatmaps'
path_codes='030_heatmaps/23_01_pca1_cl20_codes.xls'
pdf_width=15
pdf_height=10
path_cluster_merging='23_01_pca1_cl20_cluster_mergingNEW2.xlsx'


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

prefix <- codes_prefix
suffix <- ""
outdir <- codes_outdir

if(!file.exists(outdir)) 
  dir.create(outdir, recursive = TRUE)

rand_seed <- 1234

# ------------------------------------------------------------
# Load codes data
# ------------------------------------------------------------

codes <- read.table(path_codes, header = TRUE, sep = "\t", as.is = TRUE)


labels <- data.frame(cluster = sort(unique(codes$cluster)), label = factor(sort(unique(codes$cluster))))

data <- codes[, !grepl("code_id|cluster|size", colnames(codes)), drop = FALSE]

code_sizes <- codes$size

mm <- match(codes$cluster, labels$cluster)

code_clust <- labels$label[mm]

  
# ------------------------------------------------------------
# Load cluster merging data
# ------------------------------------------------------------

if(any(grepl("path_cluster_merging=", args))){
  
  # read cluster merging file
  cm <- read.xls(path_cluster_merging)
  
  cm
  
  if(!all(c("old_cluster", "label", "new_cluster") %in% colnames(cm)))
    stop("Merging file must contain 'old_cluster', 'label' and 'new_cluster' columns!")
  
  
  # remove spaces in labels bcs they are problematic...
  cm$label <- factor(cm$label, labels = gsub(" ", "_", levels(cm$label))) 
  
  # new labels
  labels <- unique(cm[,c("new_cluster","label")])
  colnames(labels) <- c("cluster", "label")
  
  labels <- labels[order(labels$cluster, decreasing = FALSE), ]
  labels$label <- factor(labels$label, levels = unique(labels$label))
  
  labels
  
  if(sum(duplicated(labels$label)) > 0 | sum(duplicated(labels$cluster)) > 0)
    stop("Wrong merging file")
  
  # new code clustering
  
  clustm <- factor(codes$cluster, levels = cm$old_cluster)
  levels(clustm) <- cm$new_cluster
  clustm <- as.numeric(as.character(clustm))
  
  mm <- match(clustm, labels$cluster)
  
  code_clust <- labels$label[mm]
  
}



# ------------------------------------------------------------
### Colors for the clusters
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

colors_clusters <- color_ramp[1:nlevels(labels$label)]
names(colors_clusters) <- levels(labels$label)

colors_clusters

# -------------------------------------
# BuildMST
# -------------------------------------


adjacency <- stats::dist(data, method = "euclidean")

fullGraph <- igraph::graph.adjacency(as.matrix(adjacency), mode = "undirected", weighted = TRUE)

MST_graph <- igraph::minimum.spanning.tree(fullGraph)

MST_l <- igraph::layout.kamada.kawai(MST_graph) 


layout <- MST_l
lty <- 1

vertex_sizes <- code_sizes/max(code_sizes) * 10
# vertex_sizes[vertex_sizes < 2] <- 2

### Plot MST

vertex_colors <- colors_clusters[code_clust]

pdf(file.path(outdir, paste0(prefix, "codes_mst.pdf")), width = 7, height = 7)

igraph::plot.igraph(MST_graph, layout = layout, vertex.size = vertex_sizes, vertex.label = NA, vertex.label.cex = 0.5, vertex.color = vertex_colors, edge.lty = lty)

dev.off()


# -------------------------------------
# Dimension reduction 
# -------------------------------------

dr <- list()

set.seed(rand_seed)
dr[["isomap"]] <- cytofkit::cytof_dimReduction(data, method="isomap")
set.seed(rand_seed)
dr[["tsne"]] <- cytofkit::cytof_dimReduction(data, method="tsne", perplexity = 20, tsneSeed = rand_seed)
set.seed(rand_seed)
dr[["diffusion"]] <- cytofkit::cytof_dimReduction(data, method="diffusionmap")
set.seed(rand_seed)
dr[["pca"]] <- cytofkit::cytof_dimReduction(data, method="pca")

dr[["graph"]] <- MST_l

dr_methods <- names(dr)



for(i in 1:length(dr)){
  # i = 1
  
  dr_data <- data.frame(dr[[i]])
  colnames(dr_data) <- c("dim1", "dim2")
  
  dr_data$cluster <- code_clust
  dr_data$size <- code_sizes
  
  ggdf <- dr_data
  
  ggp <- ggplot(ggdf,  aes(x = dim1, y = dim2, color = cluster, size = size)) +
    geom_point(alpha = 0.7) +
    labs(x = "dim1", y="dim2") + 
    theme_bw() +
    theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold"), legend.key = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_color_manual(values = colors_clusters[code_clust]) +
    guides(colour = guide_legend(override.aes = list(size = 5)))
  
  pdf(file.path(outdir, paste0(prefix, "codes_", dr_methods[[i]], ".pdf")), width = 9, height = 7)                 
  print(ggp)
  dev.off()
  
}









sessionInfo()















################################
### 02_som_codes.R done!
################################