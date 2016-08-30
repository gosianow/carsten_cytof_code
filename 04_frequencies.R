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

##############################################################################
# Test arguments
##############################################################################

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
# freq_prefix='23_01_pca1_cl20_'
# freq_outdir='050_frequencies'
# path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.xlsx'
# path_clustering='030_heatmaps/23_01_pca1_cl20_clustering.xls'
# path_clustering_labels='030_heatmaps/23_01_pca1_cl20_clustering_labels.xls'
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

color_values <- colors$color
names(color_values) <- colors$condition

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


# ---------------------------------------
# Calculate the cluster frequencies per sample
# ---------------------------------------

samp <- clustering[, "sample_id"]

# calculate frequencies
freq <- table(cluster = clust, samp)

# use labels as names of clusters
mlab <- match(rownames(freq), labels$cluster)
rownames(freq) <- labels$label[mlab]

prop <- t(t(freq) / colSums(freq)) * 100


### Save the frequencies and proportions
prop_out <- data.frame(cluster = rownames(prop), as.data.frame.matrix(prop))
freq_out <- data.frame(cluster = rownames(freq), as.data.frame.matrix(freq))

write.table(prop_out, file=file.path(outdir, paste0(prefix, "frequencies.xls")), row.names=FALSE, quote=FALSE, sep="\t")
write.table(freq_out, file=file.path(outdir, paste0(prefix, "counts.xls")), row.names=FALSE, quote=FALSE, sep="\t")



# ------------------------------------------------------------
### Plot frequencies
# ------------------------------------------------------------

ggdf <- melt(prop, value.name = "prop")

## use labels as clusters
ggdf$cluster <- factor(ggdf$cluster, levels = labels$label)

## add group info
mm <- match(ggdf$samp, md$shortname)
ggdf$group <- factor(md$condition[mm])

## replace _ with \n
levels(ggdf$group) <- gsub("_", "\n", levels(ggdf$group))

## calculate mean and sd for the error bars on the plot
ggds <- ddply(ggdf, .(group, cluster), summarise, mean = mean(prop), sd = sd(prop))

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
    scale_color_manual(values = color_values)
  
}

pdf(file.path(outdir, paste0(prefix, "frequencies.pdf")), w=5, h=4, onefile=TRUE)
for(i in seq(length(ggp)))
  print(ggp[[i]])
dev.off()


# ------------------------------------
# plot all clusters in one pdf; colors per group; points

# add more info about samples
ggdf$day <- strsplit2(ggdf$group, "\n")[, 1]
ggds$day <- strsplit2(ggds$group, "\n")[, 1]


ggp <- ggplot(ggdf, aes(x = cluster, y = prop, color = group)) +
  geom_jitter(size=2, shape = 16, width = 0.6, height = 0) +
  geom_errorbar(data=ggds, aes(x=cluster, y=mean, ymin=mean, ymax=mean), width=0.4) +
  geom_errorbar(data=ggds, aes(x=cluster, y=mean, ymin=mean-sd, ymax=mean+sd), width=0.25) +
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
  scale_color_manual(values = color_values) +
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
#     scale_fill_manual(values = color_values)
#   
# }
# 
# pdf(file.path(outdir, paste0(prefix, "frequencies_bar.pdf")), w=10, h=5, onefile=TRUE)
# for(i in seq(length(ggp)))
#   print(ggp[[i]])
# dev.off()



# ------------------------------------
# plot each cluster as a separate page in the pdf file; colors per group; barplot ordered by proportion; facet per day

# add more info about samples
ggdf$day <- strsplit2(ggdf$group, "\n")[, 1]

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
    scale_fill_manual(values = color_values) +
    facet_wrap(~ day, nrow = 1, scales = "free_x")
  
}

pdf(file.path(outdir, paste0(prefix, "frequencies_bar_facet.pdf")), w=10, h=5, onefile=TRUE)
for(i in seq(length(ggp)))
  print(ggp[[i]])
dev.off()


# ------------------------------------------------------------
# Test for frequency differences between groups
# ------------------------------------------------------------

source(path_fun_models)


# -----------------------------
### Fit a logit GLM
# -----------------------------

pvs_glm_logit <- fit_glm_logit(data = freq_out, md)

out <- data.frame(cluster = rownames(freq_out), pvs_glm_logit)

## save the results
write.table(out, file=file.path(outdir, paste0(prefix, "frequencies_pvs_glm_logit", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")

table(pvs_glm_logit$adjp_responseR < 0.05)
table(pvs_glm_logit$adjp_daytx < 0.05)


# -----------------------------
### Fit a logit GLMM 
# -----------------------------

# pvs_glmm_logit <- fit_glmm_logit(data = freq_out, md)
# 
# out <- data.frame(cluster = rownames(freq_out), pvs_glmm_logit)
# 
# ## save the results
# write.table(out, file=file.path(outdir, paste0(prefix, "frequencies_pvs_glmm_logit", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
# 
# table(pvs_glmm_logit$adjp_responseR < 0.05)
# table(pvs_glmm_logit$adjp_daytx < 0.05)










sessionInfo()













################################
### 04_frequencies done!
################################