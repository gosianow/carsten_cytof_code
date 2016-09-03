##############################################################################
## <<04_expression.R>>

# BioC 3.3
# Created 22 Aug 2016
# Updated 29 Aug 2016

##############################################################################
Sys.time()
##############################################################################

# Load packages
library(gdata)
library(Repitools)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2)
library(RColorBrewer)
library(limma)
library(pheatmap)

##############################################################################
# Test arguments
##############################################################################

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
# expr_prefix='23_01_pca1_merging5_raw_'
# expr_outdir='080_expression'
# path_data='010_data/23_01_expr_raw.rds'
# path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.xlsx'
# path_clustering_observables='030_heatmaps/23_01_pca1_clustering_observables.xls'
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

prefix <- expr_prefix
outdir <- expr_outdir
suffix <- ""

if( !file.exists(outdir) ) 
  dir.create(outdir)

source(path_fun_models)

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

# read expression data
if(grepl(".txt", path_data)){
  expr <- read.table(path_data, header = TRUE, sep = "\t", as.is = TRUE)
}
if(grepl(".rds", path_data)){
  expr <- readRDS(path_data)
}

cell_id <- expr[, "cell_id"]
samp <- expr[, "sample_id"]
fcs_colnames <- colnames(expr)[!grepl("cell_id|sample_id", colnames(expr))]
e <- expr[, fcs_colnames]

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
# Load clustering data
# ------------------------------------------------------------

# clustering
clustering <- read.table(path_clustering, header = TRUE, sep = "\t", as.is = TRUE)

clust <- clustering[, "cluster"]
names(clust) <- clustering[, "cell_id"]

# clustering labels
labels <- read.table(path_clustering_labels, header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))
rownames(labels) <- labels$cluster

# clustering observables
clustering_observables <- read.table(path_clustering_observables, header = TRUE, sep = "\t", as.is = TRUE)
rownames(clustering_observables) <- clustering_observables$mass

clust_observ <- clustering_observables[clustering_observables$clustering_observable, "mass"]


# ------------------------------------------------------------
# get the isotope and antigen for fcs markers

m <- match(fcs_colnames, clustering_observables$mass)

fcs_panel <- data.frame(colnames = fcs_colnames, Isotope = clustering_observables$mass[m], Antigen = clustering_observables$marker[m], stringsAsFactors = FALSE)


# ------------------------------------------------------------

# Indeces of observables used for clustering 
scols <- which(fcs_colnames %in% clust_observ)

# Indeces of other observables
xcols <- which(!fcs_colnames %in% clust_observ)


# ordered by decreasing pca score
if("avg_score" %in% colnames(clustering_observables)){
  scols <- scols[order(clustering_observables[fcs_colnames[scols], "avg_score"], decreasing = TRUE)]
  xcols <- xcols[order(clustering_observables[fcs_colnames[xcols], "avg_score"], decreasing = TRUE)]
}


# ------------------------------------------------------------
# Get the median expression per cluster
# ------------------------------------------------------------

colnames(e) <- fcs_panel$Antigen

a <- aggregate(e, by = list(clust, samp), FUN = median)

mlab <- match(a$Group.1, labels$cluster)
a$label <- labels$label[mlab]

colnames(a)[1:2] <- c("cluster", "sample")

a <- a[, c("cluster", "label", "sample", fcs_panel$Antigen[c(scols, xcols)])]

### Save the median expression per cluster and sample
write.table(a, file.path(outdir, paste0(prefix, "expr_clust.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# ------------------------------------------------------------
# Get the median expression for cells from all the clusters
# ------------------------------------------------------------

colnames(e) <- fcs_panel$Antigen

aa <- aggregate(e, by = list(samp), FUN = median)

colnames(aa)[1] <- c("sample")

aa$cluster <- -1
aa$label <- "all"

aa <- aa[, c("cluster", "label", "sample", fcs_panel$Antigen[c(scols, xcols)])]

### Save the median expression per cluster and sample
write.table(aa, file.path(outdir, paste0(prefix, "expr_all.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# -----------------------------------------------------------------------------
### Plot expression per cluster
# -----------------------------------------------------------------------------

alist <- list()

a$label <- factor(a$label, levels = labels$label)

alist[[1]] <- a
alist[[2]] <- aa

out_names <- c("clust", "all")

pdf_width <- c(10, 4)
pdf_height <- c(8, 4)

for(j in 1:length(alist)){
  
  aaa <- alist[[j]]
  out_name <- out_names[j]
  
  ggdf <- melt(aaa, id.vars = c("cluster", "label", "sample"), value.name = "expr", variable.name = "marker")
  
  ## use labels as clusters
  ggdf$cluster <- ggdf$label
  
  ## add group info
  mm <- match(ggdf$samp, md$shortname)
  ggdf$group <- factor(md$condition[mm])
  
  ## replace _ with \n
  levels(ggdf$group) <- gsub("_", "\n", levels(ggdf$group))
  
  ## order markers as for heatmaps
  ggdf$marker <- factor(ggdf$marker, levels = fcs_panel$Antigen[c(scols, xcols)])
  
  ## calculate mean and sd for the error bars on the plot
  ggds <- ddply(ggdf, .(group, cluster, marker), summarise, mean = mean(expr), sd = sd(expr))
  
  markers <- levels(ggdf$marker)
  
  # ------------------------------------
  ## plot each cluster as a separate page in the pdf file
  
  ggp <- list()
  
  for(i in 1:nlevels(ggdf$marker)){
    # i = 1
    
    df <- ggdf[ggdf$marker == markers[i], , drop = FALSE]
    ds <- ggds[ggds$marker == markers[i], , drop = FALSE]
    
    ggp[[i]] <- ggplot(df, aes(x = group, y = expr, color = group)) +
      geom_jitter(size=2.5, shape = 16, width = 0.5, height = 0) +
      geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean, ymax=mean), colour='black', width=0.4) +
      geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean-sd, ymax=mean+sd), colour='black', width=0.25) +
      facet_wrap(~ cluster, scales = "free") +
      ggtitle(markers[i]) +
      theme_bw() +
      ylab("Expression") +
      xlab("") +
      theme(axis.text.x = element_text(size=12, face="bold"), 
        axis.title.y = element_text(size=12, face="bold"),
        legend.position = "none") +
      scale_color_manual(values = color_values)
    
  }
  
  pdf(file.path(outdir, paste0(prefix, "expr_", out_name, ".pdf")), w = pdf_width[j], h = pdf_height[j], onefile=TRUE)
  for(i in seq(length(ggp)))
    print(ggp[[i]])
  dev.off()
  
}


# -----------------------------------------------------------------------------
# Test for expression differences between groups per cluster
# -----------------------------------------------------------------------------

source(path_fun_models)

alist <- list()
alist[[1]] <- a
alist[[2]] <- aa

out_names <- c("clust", "all")


for(j in 1:length(alist)){
  # j = 1
  
  aaa <- alist[[j]]
  out_name <- out_names[j]
  
  ### Prepare the matrix with data (rows - markers X clusters; columns - samples)
  
  expr <- aaa
  
  exprm <- melt(expr, id.vars = c("cluster", "label", "sample"), value.name = "expr", variable.name = "marker")
  
  exprc <- dcast(exprm, cluster + label + marker ~ sample, value.var = "expr")
  
  
  # -----------------------------
  # Fit a normal GLM with interactions
  # -----------------------------
  
  fit_glm_norm_inter_out <- fit_glm_norm_inter(data = exprc, md)
  
  pvs_glm_norm_inter <- data.frame(exprc[, c("cluster", "label", "marker")], fit_glm_norm_inter_out[["pvals"]])
  coeffs_glm_norm_inter <- data.frame(exprc[, c("cluster", "label", "marker")], fit_glm_norm_inter_out[["coeffs"]])
  
  oo <- order(pvs_glm_norm_inter$pval_responseR, decreasing = FALSE)
  pvs_glm_norm_inter <- pvs_glm_norm_inter[oo, ]
  coeffs_glm_norm_inter <- coeffs_glm_norm_inter[oo, ]
  
  ## save the results
  write.table(pvs_glm_norm_inter, file = file.path(outdir, paste0(prefix, "expr_", out_name, "_pvs_glm_norm_inter", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  write.table(coeffs_glm_norm_inter, file = file.path(outdir, paste0(prefix, "expr_", out_name, "_coeffs_glm_norm_inter", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  
  table(pvs_glm_norm_inter$adjp_responseR < 0.05, useNA = "always")
  table(pvs_glm_norm_inter$adjp_responseR.daytx < 0.05, useNA = "always")
  
  
  # -----------------------------
  # Fit a normal GLM
  # -----------------------------
  
  fit_glm_norm_out <- fit_glm_norm(data = exprc, md)
  
  pvs_glm_norm <- data.frame(exprc[, c("cluster", "label", "marker")], fit_glm_norm_out[["pvals"]])
  coeffs_glm_norm <- data.frame(exprc[, c("cluster", "label", "marker")], fit_glm_norm_out[["coeffs"]])
  
  oo <- order(pvs_glm_norm$pval_responseR, decreasing = FALSE)
  pvs_glm_norm <- pvs_glm_norm[oo, ]
  coeffs_glm_norm <- coeffs_glm_norm[oo, ]
  
  ## save the results
  write.table(pvs_glm_norm, file = file.path(outdir, paste0(prefix, "expr_", out_name, "_pvs_glm_norm", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  write.table(coeffs_glm_norm, file = file.path(outdir, paste0(prefix, "expr_", out_name, "_coeffs_glm_norm", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  
  table(pvs_glm_norm$adjp_responseR < 0.05, useNA = "always")
  
  # -----------------------------
  # Fit a normal GLM per base and tx
  # -----------------------------
  
  fit_glm_norm_resp_base <- fit_glm_norm_resp(data = exprc, md = md[md$day == "base", ])
  
  fit_glm_norm_resp_tx <- fit_glm_norm_resp(data = exprc, md = md[md$day == "tx", ])
  
  pvs_base <- fit_glm_norm_resp_base[["pvals"]]
  pvs_tx <- fit_glm_norm_resp_tx[["pvals"]]
  
  table(pvs_base[, "adjp_responseR"] < 0.05, useNA = "always")
  table(pvs_tx[, "adjp_responseR"] < 0.05, useNA = "always")
  
  
  # -----------------------------
  # Fit a normal GLM with interactions + test contrasts with multcomp pckg
  # -----------------------------
  
  fit_glm_norm_interglht_out <- fit_glm_norm_interglht(data = exprc, md)
  
  pvs_glm_norm_interglht <- data.frame(exprc[, c("cluster", "label", "marker")], fit_glm_norm_interglht_out[["pvals"]])
  coeffs_glm_norm_interglht <- data.frame(exprc[, c("cluster", "label", "marker")], fit_glm_norm_interglht_out[["coeffs"]])
  
  oo <- order(pvs_glm_norm_interglht$pval_NRvsR, decreasing = FALSE)
  pvs_glm_norm_interglht <- pvs_glm_norm_interglht[oo, ]
  coeffs_glm_norm_interglht <- coeffs_glm_norm_interglht[oo, ]
  
  ## save the results
  write.table(pvs_glm_norm_interglht, file = file.path(outdir, paste0(prefix, "expr_", out_name, "_pvs_glm_norm_interglht", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  write.table(coeffs_glm_norm_interglht, file = file.path(outdir, paste0(prefix, "expr_", out_name, "_coeffs_glm_norm_interglht", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  
  
  table(pvs_glm_norm_interglht$adjp_NRvsR < 0.05, useNA = "always")
  table(pvs_glm_norm_interglht$adjp_NRvsR_base < 0.05, useNA = "always")
  table(pvs_glm_norm_interglht$adjp_NRvsR_tx < 0.05, useNA = "always")
  table(pvs_glm_norm_interglht$adjp_NRvsR_basevstx < 0.05, useNA = "always")
  
  # ----------------------------------------
  # Plot a heatmap of significant cases
  # ----------------------------------------
  
  which_top_pvs <- pvs_glm_norm_interglht$adjp_NRvsR < 0.05 & !is.na(pvs_glm_norm_interglht$adjp_NRvsR)
  which(which_top_pvs)
  
  if(sum(which_top_pvs) > 0){
    
    pvs_top <- pvs_glm_norm_interglht[which_top_pvs, c("cluster", "label", "marker", "adjp_NRvsR"), drop = FALSE]
    colnames(pvs_top) <- c("cluster", "label", "marker", "adjpval")
    
    expr_heat <- merge(pvs_top, exprc, by = c("cluster", "label", "marker"), all.x = TRUE, sort = FALSE)
    
    ## group the expression by cluster
    expr_heat <- expr_heat[order(expr_heat$label), , drop = FALSE]
    
    # -----------------------------
    ### Plot one heatmap with R vs NR
    
    ## order the samples
    samples2plot <- md[md$response %in% c("NR", "R"), ]
    samples2plot <- samples2plot$shortname[order(samples2plot$response, samples2plot$day)]
    
    ## gap in the heatmap 
    gaps_col <- sum(grepl("_NR", samples2plot))
    gaps_row <- unique(cumsum(table(expr_heat$label)))
    if(length(gaps_row) == 1) 
      gaps_row <- NULL
    
    ## expression scaled by row
    expr <- t(apply(expr_heat[, samples2plot, drop = FALSE], 1, function(x) (x-mean(x))/sd(x)))
    th <- 2.5
    expr[expr > th] <- th
    expr[expr < -th] <- -th
    breaks = seq(from = -th, to = th, length.out = 101)
    legend_breaks = seq(from = -round(th), to = round(th), by = 1)
    
    labels_row <- paste0(expr_heat$label, "/ ", expr_heat$marker, " (", sprintf( "%.02e", expr_heat$adjpval), ")") 
    labels_col <- colnames(expr)
    
    pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "expr_", out_name, "_pheatmap1", suffix, ".pdf")))
    
    
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
      expr <- t(apply(expr_heat[, samples2plot_sub], 1, function(x) (x-mean(x))/sd(x)))
      th <- 2.5
      expr[expr > th] <- th
      expr[expr < -th] <- -th
      breaks = seq(from = -th, to = th, length.out = 101)
      legend_breaks = seq(from = -round(th), to = round(th), by = 1)
      
      labels_row <- paste0(expr_heat$label, "/ ", expr_heat$marker, " (", sprintf( "%.02e", expr_heat$adjpval), ")") 
      labels_col <- colnames(expr)
      
      pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "expr_", out_name, "_pheatmap_", i, suffix, ".pdf")))
      
    }
    
    
    
  }
  
  # ----------------------------------------
  # Plot coefficients NRvsR for base and tx (to show that they correlate)
  # ----------------------------------------
  
  ggdf <- coeffs_glm_norm_interglht[, c("NRvsR_base", "NRvsR_tx")]
  
  limmin <- min(ggdf, na.rm = TRUE)
  limmax <- max(ggdf, na.rm = TRUE)
  
  ggp <- ggplot(data = ggdf, aes(x = NRvsR_base, y = NRvsR_tx)) +
    geom_point(size = 3, alpha = 0.75) +
    geom_abline(intercept = 0, slope = 1) +
    coord_cartesian(xlim = c(limmin, limmax), ylim = c(limmin, limmax)) +
    theme_bw() +
    theme(axis.text = element_text(size=14), 
      axis.title = element_text(size=14, face="bold"))
  
  pdf(file.path(outdir, paste0(prefix, "expr_", out_name, "_coeffs", suffix, ".pdf")), w=5, h=5, onefile=TRUE)
  print(ggp)
  dev.off()
  
  
}







sessionInfo()














################################
### 04_expression done!
################################