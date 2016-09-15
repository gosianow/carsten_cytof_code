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
library(gtools) # for logit


##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
freq_prefix='23_01_pca1_merging6_'
freq_outdir='050_frequencies'
path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.xlsx'
path_clustering='030_heatmaps/23_01_pca1_merging6_clustering.xls'
path_clustering_labels='030_heatmaps/23_01_pca1_merging6_clustering_labels.xls'
path_fun_models='/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_models.R'

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

## clustering labels
labels <- read.table(path_clustering_labels, header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))

## clustering
clustering <- read.table(path_clustering, header = TRUE, sep = "\t", as.is = TRUE)

## drop the "drop" cluster
if("drop" %in% labels$label){
  
  clust2drop <- labels$cluster[labels$label == "drop"]
  cells2drop <- clustering$cluster != clust2drop
  clustering <- clustering[cells2drop, , drop = FALSE]
  labels <- labels[labels$label != "drop", ,drop = FALSE]
  labels$label <- factor(labels$label)
  
}

clust <- clustering[, "cluster"]

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


ggdf$day <- factor(ggdf$day)
ggds$day <- factor(ggds$day)


# ------------------------------------
### plot each cluster as a separate page in the pdf file
ggp <- list()
clusters <- levels(ggdf$cluster)

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

pdf(file.path(outdir, paste0(prefix, "frequencies_plot.pdf")), w=5, h=4, onefile=TRUE)
for(i in seq(length(ggp)))
  print(ggp[[i]])
dev.off()


# ------------------------------------
# plot all clusters in one pdf; colors per group; points; separate pdf for base and tx

days <- levels(ggdf$day)

for(i in 1:nlevels(ggdf$day)){
  # i = 1
  
  df <- ggdf[ggdf$day == days[i], , drop = FALSE]
  ds <- ggds[ggds$day == days[i], , drop = FALSE]
  
  ggp <- ggplot(df, aes(x = cluster, y = prop, color = group)) +
    geom_point(size=2, shape = 16, alpha = 0.8, position = position_jitterdodge(jitter.width = 3, jitter.height = 0, dodge.width = 0.7)) +
    geom_errorbar(data=ds, aes(x=cluster, y=mean, ymin=mean, ymax=mean), width=0.4, position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.7), size = 1) +
    geom_errorbar(data=ds, aes(x=cluster, y=mean, ymin=mean-sd, ymax=mean+sd), width=0.25, position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.7), size = 1) +
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
      legend.title = element_blank(), legend.position = "right", legend.key = element_blank()) +
    guides(color = guide_legend(ncol = 1)) +
    scale_color_manual(values = color_groups)
  
  pdf(file.path(outdir, paste0(prefix, "frequencies_plot_", days[i] ,".pdf")), w=7, h=4)
  print(ggp)
  dev.off()
  
}


# ------------------------------------------------------------
# Test for frequency differences between groups
# ------------------------------------------------------------

source(path_fun_models)

# -----------------------------
### Fit a logit GLM with interactions + test contrasts with multcomp pckg
# -----------------------------

fit_glm_logit_interglht_out <- fit_glm_logit_interglht(data = freq_out, md)

pvs_glm_logit_interglht <- data.frame(freq_out[, c("cluster", "label")], fit_glm_logit_interglht_out[["pvals"]])
coeffs_glm_logit_interglht <- data.frame(freq_out[, c("cluster", "label")], fit_glm_logit_interglht_out[["coeffs"]])

oo <- order(pvs_glm_logit_interglht$pval_NRvsR, decreasing = FALSE)
pvs_glm_logit_interglht <- pvs_glm_logit_interglht[oo, ]
coeffs_glm_logit_interglht <- coeffs_glm_logit_interglht[oo, ]

## save the results
write.table(pvs_glm_logit_interglht, file=file.path(outdir, paste0(prefix, "frequencies_pvs_glm_logit_interglht", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
write.table(coeffs_glm_logit_interglht, file=file.path(outdir, paste0(prefix, "frequencies_coeffs_glm_logit_interglht", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")

table(pvs_glm_logit_interglht$adjp_NRvsR < 0.05, useNA = "always")
table(pvs_glm_logit_interglht$adjp_NRvsR_base < 0.05, useNA = "always")
table(pvs_glm_logit_interglht$adjp_NRvsR_tx < 0.05, useNA = "always")
table(pvs_glm_logit_interglht$adjp_NRvsR_basevstx < 0.05, useNA = "always")

table(pvs_glm_logit_interglht$pval_NRvsR < 0.05, useNA = "always")
table(pvs_glm_logit_interglht$pval_NRvsR_base < 0.05, useNA = "always")
table(pvs_glm_logit_interglht$pval_NRvsR_tx < 0.05, useNA = "always")
table(pvs_glm_logit_interglht$pval_NRvsR_basevstx < 0.05, useNA = "always")

# ----------------------------------------
# Plot a heatmap of significant cases
# ----------------------------------------

### normalize the expression for base and tx separately
expr_norm <- prop_out[, grep("cluster|label|_NR|_R", colnames(prop_out))]
th <- 2.5

for(i in c("base", "tx")){
  # i = "base"
  expr_norm[, grep(i, colnames(expr_norm))] <- t(apply(expr_norm[, grep(i, colnames(expr_norm)), drop = FALSE], 1, function(x){ x <- (x-mean(x))/sd(x); x[x > th] <- th; x[x < -th] <- -th; return(x)}))
}

breaks = seq(from = -th, to = th, length.out = 101)
legend_breaks = seq(from = -round(th), to = round(th), by = 1)

### add p-value info
expr_all <- merge(pvs_glm_logit_interglht, expr_norm, by = c("cluster", "label"), all.x = TRUE, sort = FALSE)


# -----------------------------
### Plot one heatmap with R vs NR

adjpval_name <- "adjp_NRvsR"

## group the expression by cluster
expr_all <- expr_all[order(expr_all[, adjpval_name]), , drop = FALSE]

which_top_pvs <- expr_all[, adjpval_name] < 0.05 & !is.na(expr_all[, adjpval_name])
which(which_top_pvs)

if(sum(which_top_pvs) > 0) {
  
  expr_heat <- expr_all[which_top_pvs, , drop = FALSE]
  
  # -----------------------------
  ## order the samples by NR and R
  
  samples2plot <- md[md$response %in% c("NR", "R"), ]
  samples2plot <- samples2plot$shortname[order(samples2plot$response, samples2plot$day)]
  
  ## gap in the heatmap 
  gaps_col <- sum(grepl("_NR", samples2plot))
  gaps_row <- NULL
  
  ## expression scaled by row
  expr <- expr_heat[ , samples2plot, drop = FALSE]
  
  labels_row <- paste0(expr_heat$label, " (", sprintf( "%.02e", expr_heat[, adjpval_name]), ")") 
  labels_col <- colnames(expr)
  
  pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "frequencies_pheatmap1", suffix, ".pdf")))
  
  
  # -----------------------------
  ## order the samples by base and tx
  samples2plot <- md[md$response %in% c("NR", "R"), ]
  samples2plot <- samples2plot$shortname[order(samples2plot$day, samples2plot$response)]
  
  ## gap in the heatmap 
  gaps_col <- c(max(grep("base_NR", samples2plot)), rep(max(grep("base", samples2plot)), 2), max(grep("tx_NR", samples2plot)))
  gaps_row <- NULL
  
  ## expression scaled by row
  expr <- expr_heat[ , samples2plot, drop = FALSE]
  
  labels_row <- paste0(expr_heat$label, " (", sprintf( "%.02e", expr_heat[, adjpval_name]), ")") 
  labels_col <- colnames(expr)
  
  pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "frequencies_pheatmap2", suffix, ".pdf")))
  
  
}


# -----------------------------
### Plot two heatmaps with R vs NR for base and tx

for(i in c("base", "tx")){
  # i = "base"
  
  adjpval_name <- paste0("adjp_NRvsR_", i)
  
  ## group the expression by cluster
  expr_all <- expr_all[order(expr_all[, adjpval_name]), , drop = FALSE]
  
  which_top_pvs <- expr_all[, adjpval_name] < 0.05 & !is.na(expr_all[, adjpval_name])
  which(which_top_pvs)
  
  if(sum(which_top_pvs) > 0){
    
    expr_heat <- expr_all[which_top_pvs, , drop = FALSE]
    
    # -----------------------------
    ## order the samples by NR and R
    
    samples2plot <- md[md$response %in% c("NR", "R"), ]
    samples2plot <- samples2plot$shortname[order(samples2plot$response, samples2plot$day)]
    samples2plot <- samples2plot[grep(i, samples2plot)]
    
    ## gap in the heatmap 
    gaps_col <- sum(grepl("_NR", samples2plot))
    gaps_row <- NULL
    
    ## expression scaled by row
    expr <- expr_heat[ , samples2plot, drop = FALSE]
    
    labels_row <- paste0(expr_heat$label, " (", sprintf( "%.02e", expr_heat[, adjpval_name]), ")") 
    labels_col <- colnames(expr)
    
    pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "frequencies_pheatmap_", i, suffix, ".pdf")))
    
    
  }
  
}

# -----------------------------
### Plot one heatmap with R vs NR + heatmap with p-values for NRvsR_base, NRvsR_tx and NRvsR_basevstx

adjpval_name <- c("adjp_NRvsR", "adjp_NRvsR_base", "adjp_NRvsR_tx", "adjp_NRvsR_basevstx")

## group the expression by cluster
expr_all <- expr_all[order(expr_all[, adjpval_name[4]], expr_all[, adjpval_name[3]], expr_all[, adjpval_name[2]], expr_all$label), , drop = FALSE]

which_top_pvs <- rowSums(expr_all[, adjpval_name] < 0.05) > 0 & rowSums(is.na(expr_all[, adjpval_name])) == 0
which(which_top_pvs)

if(sum(which_top_pvs) > 0) {
  
  expr_heat <- expr_all[which_top_pvs, , drop = FALSE]
  
  # -----------------------------
  ## order the samples by base and tx
  samples2plot <- md[md$response %in% c("NR", "R"), ]
  samples2plot <- samples2plot$shortname[order(samples2plot$day, samples2plot$response)]
  
  ## gap in the heatmap 
  gaps_col <- c(max(grep("base_NR", samples2plot)), rep(max(grep("base", samples2plot)), 2), max(grep("tx_NR", samples2plot)))
  gaps_row <- NULL
  
  ## expression 
  expr <- expr_heat[ , samples2plot, drop = FALSE]
  
  labels_row <- paste0(expr_heat$label) 
  labels_col <- colnames(expr)
  
  pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "frequencies_pheatmap3", suffix, ".pdf")))
  
  
  pvs <- expr_heat[, adjpval_name, drop = FALSE]
  labels_col <- colnames(pvs)
  gaps_col <- NULL
  
  
  pheatmap(pvs, cellwidth = 60, cellheight = 24, color = c("grey50", "grey90"), breaks = c(0, 0.05, 1), legend_breaks = c(0, 0.05, 1), legend = FALSE, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, display_numbers = TRUE, number_format = "%.02e", number_color = "black", fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "frequencies_pheatmap3pvs", suffix, ".pdf")))
  
  
}



# ----------------------------------------
# Plot coefficients NRvsR for base and tx (to show that they correlate)
# ----------------------------------------

adjpval_name <- "adjp_NRvsR_basevstx"

ggdf <- coeffs_glm_logit_interglht[, c("NRvsR_base", "NRvsR_tx")]

limmin <- min(ggdf, na.rm = TRUE)
limmax <- max(ggdf, na.rm = TRUE)

ggdf$interaction <- factor(pvs_glm_logit_interglht[, adjpval_name] < 0.05, levels = c("FALSE", "TRUE"))


ggp <- ggplot(data = ggdf, aes(x = NRvsR_base, y = NRvsR_tx, shape = interaction)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_abline(intercept = 0, slope = 1) +
  coord_cartesian(xlim = c(limmin, limmax), ylim = c(limmin, limmax)) +
  theme_bw() +
  theme(axis.text = element_text(size=14), 
    axis.title = element_text(size=14, face="bold"))

pdf(file.path(outdir, paste0(prefix, "frequencies_coeffs", suffix, ".pdf")), w=6, h=5, onefile=TRUE)
print(ggp)
dev.off()







sessionInfo()













################################
### 04_frequencies done!
################################