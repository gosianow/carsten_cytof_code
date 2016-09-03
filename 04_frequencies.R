##############################################################################
## <<04_frequencies.R>>

# BioC 3.3
# Created 27 July 2016
# Updated 25 Aug 2016

##############################################################################
Sys.time()
##############################################################################

# Load packages
library(gdata)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2)
library(limma) # for strsplit2
library(RColorBrewer)
library(pheatmap)

##############################################################################
# Test arguments
##############################################################################

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
# freq_prefix='23_01_pca1_merging5_'
# freq_outdir='050_frequencies'
# path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.xlsx'
# path_clustering='030_heatmaps/23_01_pca1_merging5_clustering.xls'
# path_clustering_labels='030_heatmaps/23_01_pca1_merging5_clustering_labels.xls'
# path_fun_models='/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_models.R'

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
suffix <- ""
outdir <- freq_outdir

if(!file.exists(outdir)) 
  dir.create(outdir)

source(path_fun_models)

# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

md <- read.xls(path_metadata, stringsAsFactors=FALSE)

# add more info about samples
cond_split <- strsplit2(md$condition, "_")
colnames(cond_split) <- c("day", "response")

md[, c("day", "response")] <- cond_split
md$response <- factor(md$response, levels = c("NR", "R", "HD"))

### Colors 
colors <- unique(md[, c("condition", "color")])
colors$condition <- factor(colors$condition)
## replace _ with \n
levels(colors$condition) <- gsub("_", "\n", levels(colors$condition ))

color_groups <- colors$color
names(color_groups) <- colors$condition

color_samples <- md$color
names(color_samples) <- md$shortname


# ------------------------------------------------------------
# Load cluster data
# ------------------------------------------------------------

## clustering
clustering <- read.table(path_clustering, header = TRUE, sep = "\t", as.is = TRUE)

clust <- clustering[, "cluster"]
names(clust) <- clustering[, "cell_id"]


## clustering labels
labels <- read.table(path_clustering_labels, header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))

# ------------------------------------------------------------
### Colors for 20 clusters 
# ------------------------------------------------------------

# ggplot palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=60 , c=100)[1:n]
}

color_ramp <- c(colorRampPalette(brewer.pal(12,"Paired"))(12)[-c(11)],  gg_color_hue(max(1, nlevels(labels$label)-11)) )

tsne_colors <- color_ramp[1:nlevels(labels$label)]
names(tsne_colors) <- levels(labels$label)

# ---------------------------------------
# Calculate the cluster frequencies per sample
# ---------------------------------------

samp <- clustering[, "sample_id"]

# calculate frequencies
freq <- table(clust, samp)

prop <- t(t(freq) / colSums(freq)) * 100 # proportion of clusters in samples

# use labels as names of clusters
mlab <- match(rownames(freq), labels$cluster)


### Save the frequencies and proportions
prop_out <- data.frame(labels[mlab, c("cluster", "label")], as.data.frame.matrix(prop))

freq_out <- data.frame(labels[mlab, c("cluster", "label")], as.data.frame.matrix(freq))

write.table(prop_out, file=file.path(outdir, paste0(prefix, "frequencies.xls")), row.names=FALSE, quote=FALSE, sep="\t")
write.table(freq_out, file=file.path(outdir, paste0(prefix, "counts.xls")), row.names=FALSE, quote=FALSE, sep="\t")



# ------------------------------------------------------------
### Plot frequencies
# ------------------------------------------------------------

ggdf <- melt(prop_out, id.vars = c("cluster", "label"), value.name = "prop", variable.name = "samp")

## use labels as clusters
ggdf$cluster <- factor(ggdf$label, levels = labels$label)
ggdf <- ggdf[, c("cluster", "samp", "prop")]

## add group info
mm <- match(ggdf$samp, md$shortname)
ggdf$group <- factor(md$condition[mm])

## replace _ with \n
levels(ggdf$group) <- gsub("_", "\n", levels(ggdf$group))

## calculate mean and sd for the error bars on the plot
ggds <- ddply(ggdf, .(group, cluster), summarise, mean = mean(prop), sd = sd(prop))

# add more info about samples
ggdf$day <- strsplit2(ggdf$group, "\n")[, 1]
ggds$day <- strsplit2(ggds$group, "\n")[, 1]


clusters <- levels(ggdf$cluster)

# ------------------------------------
### plot each cluster as a separate page in the pdf file
ggp <- list()

for(i in 1:nlevels(ggdf$cluster)){
  # i = 1
  
  df <- ggdf[ggdf$cluster == clusters[i], , drop = FALSE]
  ds <- ggds[ggds$cluster == clusters[i], , drop = FALSE]
  
  ggp[[i]] <- ggplot(df, aes(x = group, y = prop, color = group)) +
    geom_jitter(size=2.5, shape = 17, width = 0.5, height = 0) +
    geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean, ymax=mean), color='black', width=0.4) +
    geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean-sd, ymax=mean+sd), color='black', width=0.25) +
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
      axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
      legend.position = "none") +
    scale_color_manual(values = color_groups)
  
}

pdf(file.path(outdir, paste0(prefix, "frequencies.pdf")), w=5, h=4, onefile=TRUE)
for(i in seq(length(ggp)))
  print(ggp[[i]])
dev.off()


# ------------------------------------
# plot all clusters in one pdf; colors per group; points

ggp <- ggplot(ggdf, aes(x = cluster, y = prop, color = group)) +
  # geom_jitter(size=2, shape = 16, width = 0.6, height = 0) +
  geom_point(size=2, shape = 16, position = position_jitterdodge(jitter.width = 3, jitter.height = 0, dodge.width = 0.7)) +
  geom_errorbar(data=ggds, aes(x=cluster, y=mean, ymin=mean, ymax=mean), width=0.4, position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.7)) +
  geom_errorbar(data=ggds, aes(x=cluster, y=mean, ymin=mean-sd, ymax=mean+sd), width=0.25, position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.7)) +
  theme_bw() +
  ylab("Frequency") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12, face="bold"), 
    axis.title.y = element_text(size=12, face="bold"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    axis.line.x = element_line(size = 0.5, linetype = "solid", color = "black"), 
    axis.line.y = element_line(size = 0.5, linetype = "solid", color = "black"),
    legend.title = element_blank(), legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1)) +
  scale_color_manual(values = color_groups) +
  facet_wrap(~ day, nrow = 1, scales = "free_x")

pdf(file.path(outdir, paste0(prefix, "frequencies_colors.pdf")), w=12, h=5)
print(ggp)
dev.off()


# ------------------------------------
# plot each cluster as a separate page in the pdf file; colors per group; barplot ordered by proportion

# ggp <- list()
# 
# for(i in 1:nlevels(ggdf$cluster)){
#   # i = 1
#   
#   df <- ggdf[ggdf$cluster == clusters[i], , drop = FALSE]
#   df$samp <- factor(df$samp, levels = df$samp[order(df$prop, decreasing = TRUE)])
#   
#   ggp[[i]] <- ggplot(df, aes(x = samp, y = prop, fill = group)) +
#     geom_bar(stat="identity") +
#     ggtitle(clusters[i]) +
#     theme_bw() +
#     ylab("Frequency") +
#     xlab("") +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12, face="bold"), 
#       axis.title.y = element_text(size=12, face="bold"), 
#       panel.grid.major = element_blank(), 
#       panel.grid.minor = element_blank(), 
#       panel.border = element_blank(), 
#       axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"), 
#       axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
#       legend.position = "none") +
#     scale_fill_manual(values = color_groups)
#   
# }
# 
# pdf(file.path(outdir, paste0(prefix, "frequencies_bar.pdf")), w=10, h=5, onefile=TRUE)
# for(i in seq(length(ggp)))
#   print(ggp[[i]])
# dev.off()



# ------------------------------------
# plot each cluster as a separate page in the pdf file; colors per group; barplot ordered by proportion; facet per day

ggp <- list()

for(i in 1:nlevels(ggdf$cluster)){
  # i = 1
  
  df <- ggdf[ggdf$cluster == clusters[i], , drop = FALSE]
  df$samp <- factor(df$samp, levels = df$samp[order(df$prop, decreasing = TRUE)])
  
  ggp[[i]] <- ggplot(df, aes(x = samp, y = prop, fill = group)) +
    geom_bar(stat="identity") +
    ggtitle(clusters[i]) +
    theme_bw() +
    ylab("Frequency") +
    xlab("") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12, face="bold"), 
      axis.title.y = element_text(size=12, face="bold"), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(), 
      axis.line.x = element_line(size = 0.5, linetype = "solid", color = "black"), 
      axis.line.y = element_line(size = 0.5, linetype = "solid", color = "black"),
      legend.position = "none") +
    scale_fill_manual(values = color_groups) +
    facet_wrap(~ day, nrow = 1, scales = "free_x")
  
}

pdf(file.path(outdir, paste0(prefix, "frequencies_bar_facet.pdf")), w=10, h=5, onefile=TRUE)
for(i in seq(length(ggp)))
  print(ggp[[i]])
dev.off()


# ------------------------------------
# Plot frequencies as barplots per sample

ggp <- ggplot(ggdf, aes(x = samp, y = prop, fill = cluster)) +
  geom_bar(stat="identity") +
  theme_bw() +
  ylab("Frequency") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=14, face="bold"), 
    axis.title.y = element_text(size=14, face="bold"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    axis.line.x = element_line(size = 0.5, linetype = "solid", color = "black"), 
    axis.line.y = element_line(size = 0.5, linetype = "solid", color = "black"),
    legend.key = element_blank()) +
  scale_fill_manual("", values = tsne_colors)


pdf(file.path(outdir, paste0(prefix, "frequencies_bar_sample.pdf")), w=10, h=5, onefile=TRUE)
print(ggp)
dev.off()



# ------------------------------------
# Plot frequencies as barplots per sample + facet by day

ggp <- ggplot(ggdf[!grepl("HD", ggdf$samp), ], aes(x = samp, y = prop, fill = cluster)) +
  geom_bar(stat="identity") +
  theme_bw() +
  ylab("Frequency") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=14, face="bold"), 
    axis.title.y = element_text(size=14, face="bold"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    axis.line.x = element_line(size = 0.5, linetype = "solid", color = "black"), 
    axis.line.y = element_line(size = 0.5, linetype = "solid", color = "black"),
    legend.key = element_blank()) +
  scale_fill_manual("", values = tsne_colors) +
  facet_wrap(~ day, nrow = 1, scales = "free_x")


pdf(file.path(outdir, paste0(prefix, "frequencies_bar_sample_facet.pdf")), w=10, h=5, onefile=TRUE)
print(ggp)
dev.off()



# ------------------------------------------------------------------------
# proportion of samples in clusters
# ------------------------------------------------------------------------

# prop2 <- freq / rowSums(freq) * 100
# rowSums(prop2)
# 
# prop2_out <- data.frame(labels[mlab, c("cluster", "label")], as.data.frame.matrix(prop2))
# 
# ggdf <- melt(prop2_out, id.vars = c("cluster", "label"), value.name = "prop", variable.name = "samp")
# 
# ## use labels as clusters
# ggdf$cluster <- factor(ggdf$label, levels = labels$label)
# ggdf <- ggdf[, c("cluster", "samp", "prop")]
# 
# ## add group info
# mm <- match(ggdf$samp, md$shortname)
# ggdf$group <- factor(md$condition[mm])
# 
# ## replace _ with \n
# levels(ggdf$group) <- gsub("_", "\n", levels(ggdf$group))
# 
# clusters <- levels(ggdf$cluster)
# 
# 
# # ------------------------------------
# # Plot frequencies as barplots per cluster
# 
# ggp <- ggplot(ggdf, aes(x = cluster, y = prop, fill = samp)) +
#   geom_bar(stat="identity") +
#   theme_bw() +
#   ylab("Frequency") +
#   xlab("") +
#   theme(axis.text.y = element_text(size=12, face="bold"), 
#     axis.title.y = element_text(size=12, face="bold"), 
#     panel.grid.major = element_blank(), 
#     panel.grid.minor = element_blank(), 
#     panel.border = element_blank(), 
#     axis.line.x = element_line(size = 0.5, linetype = "solid", color = "black"), 
#     axis.line.y = element_line(size = 0.5, linetype = "solid", color = "black"),
#     legend.key = element_blank()) +
#   scale_fill_manual("", values = color_samples) +
#   coord_flip()
# 
# 
# pdf(file.path(outdir, paste0(prefix, "frequencies_bar_cluster.pdf")), w=10, h=5, onefile=TRUE)
# print(ggp)
# dev.off()


# ------------------------------------------------------------
# Test for frequency differences between groups
# ------------------------------------------------------------

source(path_fun_models)


# -----------------------------
### Fit a logit GLM
# -----------------------------

fit_glm_logit_out <- fit_glm_logit(data = freq_out, md)

pvs_glm_logit <- data.frame(freq_out[, c("cluster", "label")], fit_glm_logit_out)

oo <- order(pvs_glm_logit$pval_responseR, decreasing = FALSE)
pvs_glm_logit <- pvs_glm_logit[oo, ]


## save the results
write.table(pvs_glm_logit, file=file.path(outdir, paste0(prefix, "frequencies_pvs_glm_logit", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")

table(pvs_glm_logit$adjp_responseR < 0.05, useNA = "always")



# ----------------------------------------
# Plot a heatmap of significant cases
# ----------------------------------------

which_top_pvs <- pvs_glm_logit$adjp_responseR < 0.05 & !is.na(pvs_glm_logit$adjp_responseR)
which(which_top_pvs)

if(sum(which_top_pvs) > 0){
  
  pvs_top <- pvs_glm_logit[which_top_pvs, c("cluster", "label", "adjp_responseR"), drop = FALSE]
  colnames(pvs_top) <- c("cluster", "label", "adjpval")
  
  expr_heat <- merge(pvs_top, prop_out, by = c("cluster", "label"), all.x = TRUE, sort = FALSE)
  
  # -----------------------------
  ### Plot one heatmap with R vs NR
  
  ## order the samples
  samples2plot <- md[md$response %in% c("NR", "R"), ]
  samples2plot <- samples2plot$shortname[order(samples2plot$response, samples2plot$day)]
  
  ## gap in the heatmap 
  gaps_col <- sum(grepl("_NR", samples2plot))
  
  ## expression scaled by row
  # expr <- expr_heat[, samples2plot, drop = FALSE]
  expr <- t(apply(expr_heat[, samples2plot, drop = FALSE], 1, function(x) (x-mean(x))/sd(x) ))
  th <- 2.5
  expr[expr > th] <- th
  expr[expr < -th] <- -th
  breaks = seq(from = -th, to = th, length.out = 101)
  legend_breaks = seq(from = -round(th), to = round(th), by = 1)
  
  labels_row <- paste0(expr_heat$label, " (", sprintf( "%.02e", expr_heat$adjpval), ")") 
  labels_col <- colnames(expr)
  
  pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "frequencies_pheatmap1", suffix, ".pdf")))
  
  
  # -----------------------------
  ### Plot two heatmaps with R vs NR for base and tx
  
  ## order the samples
  samples2plot <- md[md$response %in% c("NR", "R"), ]
  samples2plot <- samples2plot$shortname[order(samples2plot$response, samples2plot$day)]
  
  for(i in c("base", "tx")){
    
    samples2plot_sub <- samples2plot[grep(i, samples2plot)]
    
    ## gap in the heatmap 
    gaps_col<- sum(grepl("_NR", samples2plot_sub))
    
    ## expression scaled by row
    # expr <- expr_heat[, samples2plot_sub]
    expr <- t(apply(expr_heat[, samples2plot_sub], 1, function(x) (x-mean(x))/sd(x) ))
    th <- 2.5
    expr[expr > th] <- th
    expr[expr < -th] <- -th
    breaks = seq(from = -th, to = th, length.out = 101)
    legend_breaks = seq(from = -round(th), to = round(th), by = 1)
    
    labels_row <- paste0(expr_heat$label, " (", sprintf( "%.02e", expr_heat$adjpval), ")") 
    labels_col <- colnames(expr)
    
    pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "frequencies_pheatmap_", i, suffix, ".pdf")))
    
  }
  
  
  
}










sessionInfo()













################################
### 04_frequencies done!
################################