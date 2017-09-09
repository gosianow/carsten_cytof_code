

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
library(dplyr)
library(cowplot)
library(scales)

##############################################################################
# Test arguments
##############################################################################

args <- NULL

prefix='23_01_pca1_merging6_'
outdir='../carsten_cytof/PD1_project/CK_2016-06-23_01/099_manual/tsne'

path_metadata='../carsten_cytof/PD1_project/CK_metadata/metadata_23_01.xlsx'
path_rtsne_out='../carsten_cytof/PD1_project/CK_2016-06-23_01/040_tsnemaps/23_01_pca1_rtsne_out.rds'
path_rtsne_data='../carsten_cytof/PD1_project/CK_2016-06-23_01/040_tsnemaps/23_01_pca1_rtsne_data.xls'
path_clustering='../carsten_cytof/PD1_project/CK_2016-06-23_01/030_heatmaps/23_01_pca1_merging6_clustering.xls'
path_clustering_labels='../carsten_cytof/PD1_project/CK_2016-06-23_01/030_heatmaps/23_01_pca1_merging6_clustering_labels.xls'

dir_manual_gating_files='../carsten_cytof/PD1_project/CK_2016-06-23_01/099_manual/populations'
dir_fcs_files='../carsten_cytof/PD1_project/CK_2016-06-23_01/099_manual/data'

day='base'
response='HD'


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

md$sample_number <- 1:nrow(md)

# ------------------------------------------------------------
# Load the fcs files, add cell and sample ids, save new fcs files
# ------------------------------------------------------------


fcs_files <- md[md$day == day & md$response %in% response, c("filename", "shortname", "sample_number")]
fcs_files

# for(i in 1:nrow(fcs_files)){
#   # i = 2
# 
#   fn <- paste0(dir_fcs_files, "/", fcs_files[i, "filename"])
# 
#   man_fcs <- read.FCS(fn)
# 
#   e <- exprs(man_fcs)
# 
#   # e[, "SampleID"] <- fcs_files[i, "sample_number"]
#   # e[, "Event_length"] <- 1:nrow(e)
#   # exprs(man_fcs) <- e
# 
#   cell_ids <- matrix(1:nrow(e), ncol = 1)
#   colnames(cell_ids) <- "cell_id"
#   e <- cbind(e, cell_ids)
# 
#   sample_ids <- matrix(fcs_files[i, "sample_number"], ncol = 1, nrow = nrow(e))
#   colnames(sample_ids) <- "sample_id"
#   e <- cbind(e, sample_ids)
# 
#   out_fcs <- flowFrame(exprs = e)
#   tail(exprs(out_fcs))
#   
#   write.FCS(out_fcs, paste0(file_path_sans_ext(fn), "_new2.fcs"))
# 
# }


fcs_files <- md[md$day == day & md$response %in% response, "filename"]
fcs_files

ncells <- sapply(1:length(fcs_files), function(i){
  fcs <- read.FCS(paste0(dir_fcs_files, "/",fcs_files[i]), transformation = FALSE)
  e <- exprs(fcs)
  return(nrow(e))
})

ncells


fcs_files <- md[md$day == day & md$response %in% response, "filename"]
fcs_files

ncells <- sapply(1:length(fcs_files), function(i){
  # i = 2
  fcs <- read.FCS(paste0(dir_fcs_files, "/",file_path_sans_ext(fcs_files[i]), "_new2.fcs"), transformation = FALSE, truncate_max_range = FALSE)
  e <- exprs(fcs)
  tail(e)
  
  return(nrow(e))
})

ncells



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


cell_ids_split <- split(clustering$cell_id, factor(clustering$sample_id, levels = md$shortname))

md$cell_counts <- sapply(cell_ids_split, length)
md$cell_counts_cumsum <- cumsum(md$cell_counts)

# ------------------------------------------------------------
# Colors for tSNE maps
# ------------------------------------------------------------


# ggplot palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=60 , c=100)[1:n]
}

# color blind palette

colors_muted <- c("#DC050C", "#E8601C", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6", "#F1932D", "#F6C141", "#F7EE55", "#4EB265", "#90C987", "#CAEDAB")
color_ramp <- c(colors_muted, gg_color_hue(max(1, nlevels(labels$label) - length(colors_muted))))

colors_tsne <- color_ramp[1:nlevels(labels$label)]
names(colors_tsne) <- levels(labels$label)
colors_tsne


# ------------------------------------------------------------
# Load tSNE data
# ------------------------------------------------------------

rtsne_out <- readRDS(path_rtsne_out)

rtsne_data <- read.table(path_rtsne_data, header = TRUE, sep = "\t", as.is = TRUE)

head(rtsne_data)


### Create a data frame with t-SNE output

rtsne <- data.frame(tSNE1 = rtsne_out$Y[, 1], tSNE2 = rtsne_out$Y[, 2], sample = rtsne_data$sample_name, cell_index = rtsne_data$cell_index, stringsAsFactors = FALSE)

rtsne <- rtsne[order(rtsne$cell_index), , drop = FALSE]


### Add clustering results

mm <- match(rtsne$cell_index, clustering$cell_id)
rtsne$label <- clustering$label[mm]


# ------------------------------------------------------------
# Load manual gating results and do tSNE plots
# ------------------------------------------------------------

samps2keep <- "base_HD2"
samps2keep

population_files <- list.files(dir_manual_gating_files, pattern = ".fcs", full.names = TRUE)
population_files

name_filters <- basename(file_path_sans_ext(population_files))
name_filters


ncells_manual <- sapply(1:length(population_files), function(i){
  man_fcs <- read.FCS(population_files[i], transformation = FALSE)
  e <- exprs(man_fcs)
  return(nrow(e))
})
ncells_manual
sum(ncells_manual)


ggp_list <- list()


for(i in 1:length(population_files)){
  # i = 1
  man_fcs <- read.FCS(population_files[i], transformation = FALSE, truncate_max_range = FALSE)
  
  e <- exprs(man_fcs)
  
  # ------------------------------------------------------------
  # Identify cells to plot and create data to plot
  # ------------------------------------------------------------
  
  mang <- data.frame(e[, c("sample_id", "cell_id")])
  colnames(mang) <- c("sample_number", "cell_id")
  
  ### Add sample name
  mm <- match(mang$sample_number, md$sample_number)
  mang$sample_id <- md$shortname[mm]
  
  mang_split <- split(mang, factor(mang$sample_id, levels = md$shortname))
  sapply(mang_split, nrow)
  
  ### Update cell ids in the man gating
  mang_split_update <- lapply(1:length(mang_split), function(i){
    
    if(nrow(mang_split[[i]]) > 0){
      if(i > 1){
        mang_split[[i]]$cell_id <- mang_split[[i]]$cell_id + md$cell_counts_cumsum[i-1]
      }
    }
    return(mang_split[[i]])
  })
  
  mang_updated <- plyr::rbind.fill(mang_split_update)
  
  
  ### Add manual gating results
  
  ggdf <- rtsne
  
  ggdf$manual_gating <- as.numeric(ggdf$cell_index %in% mang_updated$cell_id)
  
  
  # add group info
  mm <- match(ggdf$sample, md$shortname)
  ggdf$group <- md$response[mm]
  
  ### Drop the "drop" cluster
  
  ggdf <- ggdf[ggdf$label != "drop", , drop = FALSE]
  ggdf$label <- factor(ggdf$label)
  
  # ----------------------------------------------------------------------
  # tSNE plots
  # ----------------------------------------------------------------------
  
  
  colors_filter <- c("#b2beb5", "#3cb371")
  names(colors_filter) <- c("0", "1")
  
  ggdf$manual_gating <- factor(ggdf$manual_gating, levels = c("0", "1"))
  
  
  suffix <- paste0("_", name_filters[i])
  
  ## one plot 
  ggp <- ggplot(ggdf, aes_string(x = "tSNE1", y = "tSNE2", color = "manual_gating")) +
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
    scale_color_manual(values = colors_filter) +
    guides(colour = guide_legend(override.aes = list(size = 5)))
  
  pdf(file.path(outdir, paste0(prefix, "tSNEone", suffix, "_type1.pdf")), width = 9, height = 7)                 
  print(ggp)
  dev.off()
  
  
  ## one plot + colored by cluster
  ggp <- ggplot() +
    geom_point(data = ggdf[ggdf[, "manual_gating"] == "0", , drop = FALSE], aes(x = tSNE1, y = tSNE2, color = label), shape = 16, alpha = 0.8) +
    geom_point(data = ggdf[ggdf[, "manual_gating"] == "1", , drop = FALSE], aes(x = tSNE1, y = tSNE2, fill = label), shape = 21, show.legend = FALSE) +
    labs(x = "t-SNE1", y="t-SNE2") + 
    theme_bw() +
    ggtitle(paste0(name_filters[i], " by manual gating")) +
    theme(axis.text = element_text(size = 14, color = "black"), 
      axis.title  = element_text(size = 16),
      title = element_text(size = 16),
      strip.text = element_text(size = 16, hjust = 0),
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
  
  ### ------------------------
  ### Plot only the samples that were used for manual gating 
  ### ------------------------
  
  ## one plot 
  ggp <- ggplot(ggdf[ggdf$sample %in% samps2keep, , drop = FALSE], aes_string(x = "tSNE1", y = "tSNE2", color = "manual_gating")) +
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
  
  pdf(file.path(outdir, paste0(prefix, "tSNEone", suffix, "_type1_subset.pdf")), width = 9, height = 7)                 
  print(ggp)
  dev.off()
  
  
  ## one plot + colored by cluster
  ggp_list[[i]] <- ggp <- ggplot() +
    geom_point(data = ggdf[ggdf[, "manual_gating"] == "0" & ggdf$sample %in% samps2keep, , drop = FALSE], aes(x = tSNE1, y = tSNE2, color = label), shape = 16, alpha = 0.8, size = 2) +
    geom_point(data = ggdf[ggdf[, "manual_gating"] == "1" & ggdf$sample %in% samps2keep, , drop = FALSE], aes(x = tSNE1, y = tSNE2, fill = label), shape = 21, size = 2, stroke = 1.5, show.legend = FALSE) +
    labs(x = "t-SNE1", y="t-SNE2") + 
    theme_bw() +
    ggtitle(paste0(name_filters[i])) +
    theme(axis.text = element_text(size = 12, color = "black"), 
      axis.title  = element_text(size = 14),
      title = element_text(size = 14),
      strip.text = element_text(size = 16, hjust = 0),
      strip.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_blank(),
      legend.text=element_text(size = 14),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_rect(size = 0.3, linetype = "solid", color = "black")) +
    scale_fill_manual(values = colors_tsne) +
    scale_color_manual(values = colors_tsne) +
    guides(color = guide_legend(override.aes = list(size = 5)))
  
  pdf(file.path(outdir, paste0(prefix, "tSNEone", suffix, "_type2_subset.pdf")), width = 9, height = 7)                 
  print(ggp)
  dev.off()
  
  
}


### Plot all gates in one panel

legend <- get_legend(ggp_list[[1]])

ggp_list_new <- lapply(ggp_list, function(x){
  x + theme(legend.position="none")
})

ggp_list_new[["legend"]] <- legend


ggp <- plot_grid(plotlist = ggp_list_new, nrow = 2)

pdf(file.path(outdir, paste0(prefix, "tSNEone_type2_subset.pdf")), width = 12, height = 6)                 
print(ggp)
dev.off()




# ------------------------------------------------------------
# Make a table with automated clustering and manual merging cell counts
# ------------------------------------------------------------


mang_merge <- lapply(1:length(population_files), function(i){
  
  # i = 1
  man_fcs <- read.FCS(population_files[i], transformation = FALSE, truncate_max_range = FALSE)
  
  e <- exprs(man_fcs)
  
  # ------------------------------------------------------------
  # Identify cells to plot and create data to plot
  # ------------------------------------------------------------
  
  mang <- data.frame(e[, c("sample_id", "cell_id")])
  colnames(mang) <- c("sample_number", "cell_id")
  
  ### Add sample name
  mm <- match(mang$sample_number, md$sample_number)
  mang$sample_id <- md$shortname[mm]
  
  mang_split <- split(mang, factor(mang$sample_id, levels = md$shortname))
  sapply(mang_split, nrow)
  
  ### Update cell ids in the man gating
  mang_split_update <- lapply(1:length(mang_split), function(i){
    
    if(nrow(mang_split[[i]]) > 0){
      if(i > 1){
        mang_split[[i]]$cell_id <- mang_split[[i]]$cell_id + md$cell_counts_cumsum[i-1]
      }
    }
    return(mang_split[[i]])
  })
  
  mang_updated <- plyr::rbind.fill(mang_split_update)
  mang_updated$manual_gate <- name_filters[i]
  
  
  return(mang_updated)
  
})

mang_merge <- plyr::rbind.fill(mang_merge)

### Add manual merging to the clustering results

mm <- match(clustering$cell_id, mang_merge$cell_id)
clustering$manual_gate <- mang_merge$manual_gate[mm]


clustering_sub <- clustering[clustering$sample_id %in% samps2keep, ]
table(clustering_sub$manual_gate, useNA = "always")

clustering_sub$manual_gate[is.na(clustering_sub$manual_gate)] <- "unknown"


tbl <- table(clustering_sub$manual_gate, as.character(clustering_sub$label))
tbl <- as.data.frame.matrix(tbl)
tbl$manual_gate <- rownames(tbl)
tbl <- tbl %>% select("manual_gate", everything())
tbl

write.table(tbl , file = file.path(outdir, paste0(prefix, "table.xls")), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


# ------------------------------------------------------------
# Plot with frequencies of cells from manual gating
# ------------------------------------------------------------


colors_manual_gate <- c(colors_tsne, "grey")
names(colors_manual_gate) <- c(names(colors_tsne), "unknown")

freqs <- table(clustering_sub$manual_gate)
props <- freqs/sum(freqs) * 100


ggdf <- data.frame(manual_gate = names(props), proportion = as.numeric(props), sample_id = samps2keep)

ggdf$manual_gate <- factor(ggdf$manual_gate, levels = c(as.character(labels$label), "unknown"))


ylabels <- ggdf[order(ggdf$manual_gate, decreasing = TRUE), "proportion"]

### Bar plot 

ggp_prop <- ggplot(ggdf, aes(x = sample_id, y = proportion, fill = manual_gate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(y = ylabels/2 + c(0, cumsum(ylabels)[-length(proportion)]), label = percent(ylabels/100)), size = 4) +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"), 
    axis.title  = element_text(size = 14),
    title = element_text(size = 14),
    strip.text = element_text(size = 16, hjust = 0),
    strip.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.text=element_text(size = 14),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(size = 0.3, linetype = "solid", color = "black")) +
  scale_fill_manual(values = colors_manual_gate) +
  guides(color = guide_legend(override.aes = list(size = 5)))


ggp_list_new <- lapply(ggp_list, function(x){
  x + theme(legend.position="none")
})

ggp_list_new[["props"]] <- ggp_prop

ggp <- plot_grid(plotlist = ggp_list_new, nrow = 2)

pdf(file.path(outdir, paste0(prefix, "tSNEone_type2_subset_with_props.pdf")), width = 12, height = 6)                 
print(ggp)
dev.off()


### Pie plot

ggp_prop <- ggplot(ggdf, aes(x = sample_id, y = proportion, fill = manual_gate)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  geom_text(aes(y = ylabels/2 + c(0, cumsum(ylabels)[-length(proportion)]), label = percent(ylabels/100)), size = 4) +
  theme_bw() +
  theme(axis.text = element_blank(), 
    axis.title  = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    title = element_text(size = 14),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14)) +
  scale_fill_manual(values = colors_manual_gate) +
  guides(color = guide_legend(override.aes = list(size = 5)))


ggp_list_new[["props"]] <- ggp_prop + theme(legend.position="none")

legend <- get_legend(ggp_prop)

ggp_list_new[["legend"]] <- legend

ggp <- plot_grid(plotlist = ggp_list_new, nrow = 2)

pdf(file.path(outdir, paste0(prefix, "tSNEone_type2_subset_with_props_pie.pdf")), width = 14, height = 6)                 
print(ggp)
dev.off()













sessionInfo()











