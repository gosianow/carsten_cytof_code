##############################################################################
## <<08_frequencies_merged_low_freq_out.R>>

# BioC 3.3
# Created 7 Nov 2016
# Updated 7 nov 2016

# Remove clusters with very low number of cells and re-adjust the p-values
# Use qvalue for the p-value adjustment (NOT anymore!!! - problems when there are too few test - the histograms are not uniform and smooth)

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
library(tools)
library(qvalue)

##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-merged_23_29/02_CD4'
freq_prefix='23_29_CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_cl100_out001_'
freq_outdir='10_cytokines_merged/03_frequencies_auto_2responses'

path_metadata=c('/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_02.xlsx','/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_29_02.xlsx')

path_counts=c('/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-merged_23_29/02_CD4/10_cytokines_merged/01_clustering/23_29_CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_cl100_data23_counts.xls','/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-merged_23_29/02_CD4/10_cytokines_merged/01_clustering/23_29_CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_cl100_data29_counts.xls')

data_name=c('data23','data29')

path_pvs='10_cytokines_merged/03_frequencies_auto_2responses/23_29_CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_cl100_frequencies_pvs_glmer_binomial_interglht.xls'

model2fit='glmer_binomial_interglht'
min_freq=0.01
pdf_hight=8

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

if(!file.exists(rwd)) 
  dir.create(rwd, recursive = TRUE)
setwd(rwd)

prefix <- freq_prefix
suffix <- ""
outdir <- freq_outdir

if(!file.exists(outdir)) 
  dir.create(outdir, recursive = TRUE)

if(!any(grepl("pdf_hight=", args))){
  pdf_hight=4
}

# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

md <- lapply(1:length(data_name), function(i){
  
  path <- path_metadata[i]
  md <- read.xls(path, stringsAsFactors=FALSE)
  md$data <- data_name[i]
  md
  
})

md <- rbind.fill(md)

rownames(md) <- md$shortname

# add more info about samples
cond_split <- strsplit2(md$condition, "_")
colnames(cond_split) <- c("day", "response")

md[, c("day", "response")] <- cond_split
md$response <- factor(md$response, levels = c("NR", "R", "HD"))
md$day <- factor(md$day, levels = c("base", "tx"))
md$patient_id <- factor(md$patient_id)

md$data <- factor(md$data, levels = data_name)

md$data_day <- interaction(md$data, md$day, lex.order = TRUE, drop = TRUE)


### Colors 
colors <- unique(md[, c("condition", "color")])
colors$condition <- factor(colors$condition)
## replace _ with \n
levels(colors$condition) <- gsub("_", "\n", levels(colors$condition ))

color_groups <- colors$color
names(color_groups) <- colors$condition

color_groupsb <- adjustcolor(color_groups, alpha = 0.3)
names(color_groupsb) <- colors$condition

color_samples <- md$color
names(color_samples) <- md$shortname

colors <- unique(md[, c("response", "color")])
color_response <- colors$color
names(color_response) <- colors$response


# ---------------------------------------
# Load cluster frequencies 
# ---------------------------------------

freq <- lapply(1:length(data_name), function(i){
  # i = 1
  path <- path_counts[i]
  freq <- read.table(path, header = TRUE, sep = "\t", as.is = TRUE)
  freq[, -which(colnames(freq) == "cluster")]
  
})

freq_out <- Reduce(function(...) merge(..., by = c("label"), all=TRUE, sort = FALSE), freq)

freq_out <- freq_out[freq_out$label != "drop", , drop = FALSE]
freq_out$cluster <- 1:nrow(freq_out)

freq_out <- freq_out[, c("cluster", "label", md$shortname)]

# if(any(!complete.cases(freq_out)))
#   stop(paste0("Files: ", paste(basename(path_counts), collapse = ", "), " have different clusters!!!"))
if(all(!complete.cases(freq_out)))
  stop("There is no common cluster in the merged data sets!")

prop_out <- freq_out
prop_out[md$shortname] <- t(t(freq_out[md$shortname]) / colSums(freq_out[md$shortname], na.rm = TRUE))


labels <- data.frame(cluster = freq_out$cluster, label = freq_out$label)
labels$label <- factor(labels$label, levels = unique(labels$label))
labels

# ---------------------------------------
# Keep only the clusters with high enough number of cells
# ---------------------------------------


clust2keep_freq <- rowSums(freq_out[md$shortname], na.rm = TRUE) /sum(freq_out[md$shortname], na.rm = TRUE) > min_freq

prop_out <- prop_out[clust2keep_freq, , drop = FALSE]
freq_out <- freq_out[clust2keep_freq, , drop = FALSE]

labels <- data.frame(cluster = freq_out$cluster, label = freq_out$label)
labels$label <- factor(labels$label, levels = unique(labels$label))
labels


# ------------------------------------------------------------
# Plot frequencies
# ------------------------------------------------------------

ggdf <- melt(prop_out[complete.cases(freq_out), , drop = FALSE], id.vars = c("cluster", "label"), value.name = "prop", variable.name = "samp")

## use labels as clusters
ggdf$cluster <- factor(ggdf$label, levels = labels$label)
ggdf <- ggdf[, c("cluster", "samp", "prop")]

## add group info
mm <- match(ggdf$samp, md$shortname)
ggdf$group <- factor(md$condition[mm])
## add data info
ggdf$data <- factor(md$data[mm], levels = data_name)

## replace _ with \n
levels(ggdf$group) <- gsub("_", "\n", levels(ggdf$group))

## calculate mean and sd for the error bars on the plot
ggds <- ddply(ggdf, .(group, cluster, data), summarise, mean = mean(prop), sd = sd(prop))

## add more info about samples
ggdf$day <- strsplit2(ggdf$group, "\n")[, 1]
ggds$day <- strsplit2(ggds$group, "\n")[, 1]

ggdf$day <- factor(ggdf$day)
ggds$day <- factor(ggds$day)


# ------------------------------------
# plot all clusters in one pdf; colors per response; boxplots + points


ggp <- ggplot(ggdf, aes(x = cluster, y = prop, color = group, shape = data, fill = group)) +
  geom_boxplot(width = 1, position = position_dodge(width = 0.9), outlier.colour = NA) +
  geom_point(size=1.5, alpha = 0.8, position = position_jitterdodge(jitter.width = 1.2, jitter.height = 0, dodge.width = 0.9)) +
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
  scale_color_manual(values = color_groups) +
  scale_fill_manual(values = color_groupsb) +
  facet_wrap(~ day)

pdf(file.path(outdir, paste0(prefix, "frequencies_plot.pdf")), w = nlevels(ggdf$cluster) + 3, h = pdf_hight)
print(ggp)
dev.off()



# ------------------------------------
# plot all clusters in one pdf; colors per response; separate pdf for base and tx; one boxplot + points

days <- levels(ggdf$day)

for(i in 1:nlevels(ggdf$day)){
  # i = 1
  
  df <- ggdf[ggdf$day == days[i], , drop = FALSE]
  
  ggp <- ggplot(df) +
    geom_boxplot(aes(x = cluster, y = prop, color = group, fill = group), width = 1, position = position_dodge(width = 0.9), outlier.colour = NA) +
    geom_point(aes(x = cluster, y = prop, color = group, shape = data), size=2, alpha = 0.8, position = position_jitterdodge(jitter.width = 1.2, jitter.height = 0, dodge.width = 0.9)) +
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
    scale_color_manual(values = color_groups) +
    scale_fill_manual(values = color_groupsb)
  
  pdf(file.path(outdir, paste0(prefix, "frequencies_plot_boxplotpoints_", days[i] ,".pdf")), w = nlevels(ggdf$cluster)/2 + 3, h = pdf_hight)
  print(ggp)
  dev.off()
  
}


# ------------------------------------------------------------
# Re-adjust the p-values
# ------------------------------------------------------------

k <- model2fit 

### p-value for sorting the output
pval_name1 <- "pval_NRvsR"
### p-value for plotting the pheatmap2
adjpval_name2 <- "adjp_NRvsR"
pval_name2 <- "pval_NRvsR"
### p-value for plotting the pheatmap3
adjpval_name_list <- c("adjp_NRvsR", "adjp_NRvsR_base", "adjp_NRvsR_tx", "adjp_NRvsR_basevstx")
pval_name_list <- c("pval_NRvsR", "pval_NRvsR_base", "pval_NRvsR_tx", "pval_NRvsR_basevstx")

# ----------------------------------------
# Read in the p-values and re-adjust them using qvalue()
# ----------------------------------------

# pvs <- read.table(path_pvs, header = TRUE, sep = "\t", as.is = TRUE)
# 
# pvs <- pvs[pvs$label %in% prop_out$label, , drop = FALSE]
# pval_colnames <- colnames(pvs)[grep("pval_", colnames(pvs))]
# 
# pdf(file.path(outdir, paste0(prefix, "frequencies_", k, "_pvs_hist", suffix, ".pdf")))
# 
# for(i in 1:length(pval_colnames)){
#   # i = 1
#   ### get q-value object
#   qobj <- qvalue(pvs[, pval_colnames[i]])
#   
#   hist(pvs[, pval_colnames[i]], breaks = 20)
#   
#   plot(qobj)
# 
#   print(hist(qobj))
#   
#   pvs[, gsub("pval_", "adjp_", pval_colnames[i])] <- qobj$qvalues
#   
# }
# 
# dev.off()
# 
# oo <- order(pvs[, pval_name], decreasing = FALSE)
# pvs <- pvs[oo, , drop = FALSE]
# 
# ## save the results
# write.table(pvs, file=file.path(outdir, paste0(prefix, "frequencies_pvs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")

# ----------------------------------------
# Read in the p-values and re-adjust them using qvalue()
# ----------------------------------------

pvs <- read.table(path_pvs, header = TRUE, sep = "\t", as.is = TRUE)

pvs <- pvs[pvs$label %in% prop_out$label, , drop = FALSE]

pval_colnames <- colnames(pvs)[grep("pval_", colnames(pvs))]


for(i in 1:length(pval_colnames)){
  # i = 1
  
  pvs[, gsub("pval_", "adjp_", pval_colnames[i])] <- p.adjust(pvs[, pval_colnames[i]], method = "BH")
  
}

oo <- order(pvs[, pval_name1], decreasing = FALSE)
pvs <- pvs[oo, , drop = FALSE]

## save the results
write.table(pvs, file=file.path(outdir, paste0(prefix, "frequencies_pvs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")


# ----------------------------------------
# Plot a heatmap of significant cases - transform proportions with arcsin-sqrt so the dispersion is the same for low and high props.
# ----------------------------------------

### normalize the expression
ass_freq_out <- prop_out
ass_freq_out[md$shortname] <- asin(sqrt(prop_out[md$shortname]))

ass_freq_out <- ass_freq_out[complete.cases(ass_freq_out), , drop = FALSE]

expr_norm <- ass_freq_out[, c("cluster", "label", md[md$response != "HD", "shortname"])]
th <- 2.5

data_days <- levels(md$data_day)

### Normalized to mean = 0 and sd = 1 per data and day
for(i in data_days){
  # i = "data23.base"
  expr_norm[, md[md$response != "HD" & md$data_day == i, "shortname"]] <- t(apply(expr_norm[, md[md$response != "HD" & md$data_day == i, "shortname"], drop = FALSE], 1, function(x){
    sdx <- sd(x, na.rm = TRUE)
    if(sdx == 0)
      x <- (x-mean(x, na.rm = TRUE))
    else
      x <- (x-mean(x, na.rm = TRUE))/sdx
    
    x[x > th] <- th
    x[x < -th] <- -th
    
    return(x)}))
}

breaks = seq(from = -th, to = th, length.out = 101)
legend_breaks = seq(from = -round(th), to = round(th), by = 1)

### add p-value info
expr_all <- merge(expr_norm, pvs, by = c("cluster", "label"), all.x = TRUE, sort = FALSE)


# -----------------------------
### Plot one heatmap with R vs NR

## group the expression by cluster
adjpval_name <- adjpval_name2
pval_name <- pval_name2
expr_all <- expr_all[order(expr_all[, pval_name]), , drop = FALSE]

which_top_pvs <- expr_all[, adjpval_name] < 0.05 & !is.na(expr_all[, adjpval_name])
which(which_top_pvs)

if(sum(which_top_pvs) > 0) {
  print("Plot pheatmap2")
  
  expr_heat <- expr_all[which_top_pvs, , drop = FALSE]
  
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
  
  annotation_col <- data.frame(response = factor(md[colnames(expr), "response"]))
  rownames(annotation_col) <- colnames(expr)
  
  annotation_colors <- list(response = color_response[levels(annotation_col$response)])
  
  pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_col = annotation_col, annotation_colors = annotation_colors, annotation_legend = FALSE, filename = file.path(outdir, paste0(prefix, "frequencies_", k, "_pheatmap2", suffix, ".pdf")))
  
  
}


# -----------------------------
### Plot two heatmaps with R vs NR for base and tx

for(i in levels(md$day)){
  # i = "base"
  print(paste0("Plot pheatmap_", i))
  
  adjpval_name <- paste0("adjp_NRvsR_", i)
  pval_name <- paste0("pval_NRvsR_", i)
  
  ## group the expression by cluster
  expr_all <- expr_all[order(expr_all[, pval_name]), , drop = FALSE]
  
  which_top_pvs <- expr_all[, adjpval_name] < 0.05 & !is.na(expr_all[, adjpval_name])
  which(which_top_pvs)
  
  if(sum(which_top_pvs) > 0){
    
    expr_heat <- expr_all[which_top_pvs, , drop = FALSE]
    
    # -----------------------------
    ## order the samples by NR and R
    
    samples2plot <- md[md$response %in% c("NR", "R"), ]
    samples2plot <- samples2plot[order(samples2plot$response, samples2plot$day), ]
    samples2plot <- samples2plot[grep(i, samples2plot$day), ]
    samples2plot <- samples2plot$shortname
    
    ## gap in the heatmap
    gaps_col <- sum(grepl("_NR", samples2plot))
    gaps_row <- NULL
    
    ## expression scaled by row
    expr <- expr_heat[ , samples2plot, drop = FALSE]
    
    labels_row <- paste0(expr_heat$label, " (", sprintf( "%.02e", expr_heat[, adjpval_name]), ")")
    labels_col <- colnames(expr)
    
    annotation_col <- data.frame(response = factor(md[colnames(expr), "response"]))
    rownames(annotation_col) <- colnames(expr)
    
    annotation_colors <- list(response = color_response[levels(annotation_col$response)])
    
    pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_col = annotation_col, annotation_colors = annotation_colors, annotation_legend = FALSE, filename = file.path(outdir, paste0(prefix, "frequencies_", k, "_pheatmap_", i, suffix, ".pdf")))
    
    
  }
  
}

# -----------------------------
### Plot one heatmap with R vs NR + heatmap with p-values for NRvsR_base, NRvsR_tx and NRvsR_basevstx

## group the expression by cluster and order by adjpval
for(i in length(adjpval_name_list):1){
  expr_all <- expr_all[order(expr_all[, pval_name_list[i]]), , drop = FALSE]
}


# which_top_pvs <- rowSums(expr_all[, adjpval_name_list] < 0.05, na.rm = TRUE) > 0 & rowSums(is.na(expr_all[, adjpval_name_list])) == 0
# which_top_pvs <- rowSums(is.na(expr_all[, adjpval_name_list])) < length(adjpval_name_list)
which_top_pvs <- rep(TRUE, nrow(expr_all))
which(which_top_pvs)

if(sum(which_top_pvs) > 0){
  print("Plot pheatmap3")
  
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
  
  annotation_col <- data.frame(response = factor(md[colnames(expr), "response"]))
  rownames(annotation_col) <- colnames(expr)
  
  annotation_colors <- list(response = color_response[levels(annotation_col$response)])
  
  pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_col = annotation_col, annotation_colors = annotation_colors, annotation_legend = FALSE, filename = file.path(outdir, paste0(prefix, "frequencies_", k, "_pheatmap3", suffix, ".pdf")))
  
  
  pvs_heat <- expr_heat[, adjpval_name_list, drop = FALSE]
  
  labels_col <- colnames(pvs_heat)
  gaps_col <- NULL
  
  
  pheatmap(pvs_heat, cellwidth = 60, cellheight = 24, color = c("grey50", "grey70", "grey90"), breaks = c(0, 0.05, 0.1, 1), legend_breaks = c(0, 0.05, 0.1, 1), legend = FALSE, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, display_numbers = TRUE, number_format = "%.02e", number_color = "black", fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "frequencies_", k, "_pheatmap3pvs", suffix, ".pdf")))
  
  
}










sessionInfo()




















################################################################
### 08_frequencies_merged_low_freq_out.R done!
################################################################