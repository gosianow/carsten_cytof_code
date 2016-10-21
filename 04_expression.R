##############################################################################
## <<04_expression.R>>

# BioC 3.3
# Created 22 Aug 2016
# Updated 13 Oct 2016

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
# expr_prefix='23_01_pca1_merging6_raw_'
# expr_outdir='080_expression'
# path_data='010_data/23_01_expr_raw.rds'
# path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.xlsx'
# path_clustering_observables='030_heatmaps/23_01_pca1_clustering_observables.xls'
# path_clustering='030_heatmaps/23_01_pca1_merging6_clustering.xls'
# path_clustering_labels='030_heatmaps/23_01_pca1_merging6_clustering_labels.xls'
# path_fun_models='/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_models.R'
# analysis_type='all'

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-03-23_01'
expr_prefix='0323_01_pca1_merging_raw_'
expr_outdir='080_expression'
path_data='010_data/0323_01_expr_raw.rds'
path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_0323_01.xlsx'
path_clustering_observables='030_heatmaps/0323_01_pca1_clustering_observables.xls'
path_clustering='030_heatmaps/0323_01_pca1_merging_clustering.xls'
path_clustering_labels='030_heatmaps/0323_01_pca1_merging_clustering_labels.xls'
path_fun_models='/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_models.R'
analysis_type='clust'


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

prefix <- expr_prefix
outdir <- expr_outdir
suffix <- ""

if( !file.exists(outdir) ) 
  dir.create(outdir)

source(path_fun_models)

if(!analysis_type %in% c("clust", "all"))
  stop("analysis_type must be 'all' or 'clust'!")

out_name <- analysis_type

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
md$response <- factor(md$response)
md$day <- factor(md$day, levels = c("base", "tx"))
md$day <- factor(md$day)
md$patient_id <- factor(md$patient_id)



### Colors 
colors <- unique(md[, c("condition", "color")])
colors$condition <- factor(colors$condition)
## replace _ with \n
levels(colors$condition) <- gsub("_", "\n", levels(colors$condition ))

color_groups <- colors$color
names(color_groups) <- colors$condition


# ------------------------------------------------------------
# Load clustering data
# ------------------------------------------------------------

# clustering observables
clustering_observables <- read.table(path_clustering_observables, header = TRUE, sep = "\t", as.is = TRUE)
rownames(clustering_observables) <- clustering_observables$mass

clust_observ <- clustering_observables[clustering_observables$clustering_observable, "mass"]


# clustering labels
labels <- read.table(path_clustering_labels, header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))
rownames(labels) <- labels$cluster

# clustering
clustering <- read.table(path_clustering, header = TRUE, sep = "\t", as.is = TRUE)

## drop the "drop" cluster
if("drop" %in% labels$label){
  
  clust2drop <- labels$cluster[labels$label == "drop"]
  cells2drop <- clustering$cluster != clust2drop
  clustering <- clustering[cells2drop, , drop = FALSE]
  labels <- labels[labels$label != "drop", ,drop = FALSE]
  labels$label <- factor(labels$label)
  e <- e[cells2drop, , drop = FALSE]
  
}


clust <- clustering[, "cluster"]

samp <- clustering[, "sample_id"]


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

if(analysis_type == "clust"){
  
  colnames(e) <- fcs_panel$Antigen
  
  a <- aggregate(e, by = list(clust, samp), FUN = median)
  
  mlab <- match(a$Group.1, labels$cluster)
  a$label <- labels$label[mlab]
  
  colnames(a)[1:2] <- c("cluster", "sample")
  
  a <- a[, c("cluster", "label", "sample", fcs_panel$Antigen[c(scols, xcols)])]
  a$label <- factor(a$label, levels = labels$label)
  
  ### Save the median expression per cluster and sample
  write.table(a, file.path(outdir, paste0(prefix, "expr_clust.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
}

if(analysis_type == "all"){
  
  colnames(e) <- fcs_panel$Antigen
  
  a <- aggregate(e, by = list(samp), FUN = median)
  
  colnames(a)[1] <- c("sample")
  
  a$cluster <- -1
  a$label <- "all"
  
  a <- a[, c("cluster", "label", "sample", fcs_panel$Antigen[c(scols, xcols)])]
  a$label <- factor(a$label)
  
  ### Save the median expression per cluster and sample
  write.table(a, file.path(outdir, paste0(prefix, "expr_all.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
}


# -----------------------------------------------------------------------------
### Plot expression per cluster
# -----------------------------------------------------------------------------


ggdf <- melt(a, id.vars = c("cluster", "label", "sample"), value.name = "expr", variable.name = "marker")

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

clusters <- levels(ggdf$cluster)

# ------------------------------------
## plot each cluster as a separate page in the pdf file

ggp <- list()

for(i in 1:nlevels(ggdf$cluster)){
  # i = 1
  
  df <- ggdf[ggdf$cluster == clusters[i], , drop = FALSE]
  ds <- ggds[ggds$cluster == clusters[i], , drop = FALSE]
  
  ggp[[i]] <- ggplot(df, aes(x = group, y = expr, color = group)) +
    geom_jitter(size=2.5, shape = 16, width = 0.5, height = 0) +
    geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean, ymax=mean), colour='black', width=0.4) +
    geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean-sd, ymax=mean+sd), colour='black', width=0.25) +
    facet_wrap(~ marker, scales = "free") +
    ggtitle(clusters[i]) +
    theme_bw() +
    ylab("Expression") +
    xlab("") +
    theme(axis.text.x = element_text(size=12, face="bold"), 
      axis.title.y = element_text(size=12, face="bold"),
      legend.position = "none") +
    scale_color_manual(values = color_groups)
  
}

pdf(file.path(outdir, paste0(prefix, "expr_", out_name, ".pdf")), w = 18, h = 12, onefile=TRUE)
for(i in seq(length(ggp)))
  print(ggp[[i]])
dev.off()



# -----------------------------------------------------------------------------
# Test for marker expression differences between groups overall and per cluster 
# -----------------------------------------------------------------------------

source(path_fun_models)

levels(md$day)
levels(md$response)


if(identical(levels(md$day), c("base", "tx")) && identical(levels(md$response), c("NR", "R", "HD"))){
  ## create formulas
  formula_lm <- y ~ response + day + response:day
  formula_lmer <- y ~ response + day + response:day + (1|patient_id)
  
  ## create contrasts
  contrast_names <- c("NRvsR", "NRvsR_base", "NRvsR_tx", "NRvsR_basevstx")
  k0 <- c(0, 1, 0, 0, 1/2, 0) # NR vs R
  k1 <- c(0, 1, 0, 0, 0, 0) # NR vs R in base
  k2 <- c(0, 1, 0, 0, 1, 0) # NR vs R in tx
  k3 <- c(0, 0, 0, 0, 1, 0) # whether NR vs R is different in base and tx
  K <- matrix(c(k0, k1, k2, k3), nrow = 4, byrow = TRUE)
  rownames(K) <- contrast_names
  
  ### p-value for sorting the output
  pval_name <- "pval_NRvsR"
  ### p-value for plotting the pheatmap2
  adjpval_name <- "adjp_NRvsR"
  ### p-value for plotting the pheatmap3
  adjpval_name_list <- c("adjp_NRvsR", "adjp_NRvsR_base", "adjp_NRvsR_tx", "adjp_NRvsR_basevstx")
  
  
}else if(identical(levels(md$day), "base") && identical(levels(md$response), c("NR", "R", "HD"))){
  ## create formulas
  formula_lm <- y ~ response
  formula_lmer <- y ~ response + (1|patient_id)
  
  ## create contrasts
  contrast_names <- c("NRvsR_base")
  k1 <- c(0, 1, 0) # NR vs R in base
  K <- matrix(k1, nrow = 1, byrow = TRUE)
  rownames(K) <- contrast_names
  
  ### p-value for sorting the output
  pval_name <- "pval_NRvsR_base"
  ### p-value for plotting the pheatmap2
  adjpval_name <- NULL
  ### p-value for plotting the pheatmap3
  adjpval_name_list <- "adjp_NRvsR_base"
  
}else{
  stop("Metadata does not fit to any the models that are specified !!!")  
}



### Prepare the matrix with data (rows - markers X clusters; columns - samples)

expr <- a
exprm <- melt(expr, id.vars = c("cluster", "label", "sample"), value.name = "expr", variable.name = "marker")
exprc <- dcast(exprm, cluster + label + marker ~ sample, value.var = "expr")

models2fit <- c("rlm_interglht", "lm_interglht", "test_wilcoxon", "lmer_interglht")


for(k in models2fit){
  # k = "lm_interglht"
  print(k)
  
  
  switch(k,
    lm_interglht = {
      
      # Fit a LM with interactions + test contrasts with multcomp pckg
      fit_out <- fit_lm_interglht(data = exprc, md, method = "lm", formula = formula_lm, K = K)
      
    }, 
    lmer_interglht = {
      
      # Fit a lmer with interactions + test contrasts with multcomp pckg
      fit_out <- fit_lmer_interglht(data = exprc, md, formula = formula_lmer, K = K)
      
    },
    lmrob_interglht = {
      ## Problems with running lmrob!!!
      fit_out <- fit_lm_interglht(data = exprc, md, method = "lmrob", formula = formula_lm, K = K)
    },
    rlm_interglht = {
      fit_out <- fit_lm_interglht(data = exprc, md, method = "rlm", formula = formula_lm, K = K)
    },
    test_wilcoxon = {
      fit_out <- test_wilcoxon(data = exprc, md)
    }
  
  )
  
  
  # ----------------------------------------
  # Extract p-values and coeffs
  # ----------------------------------------
  
  pvs <- data.frame(exprc[, c("cluster", "label", "marker")], fit_out[["pvals"]])
  coeffs <- data.frame(exprc[, c("cluster", "label", "marker")], fit_out[["coeffs"]])
  
  oo <- order(pvs[, pval_name], decreasing = FALSE)
  pvs <- pvs[oo, , drop = FALSE]
  coeffs <- coeffs[oo, , drop = FALSE]
  
  ## save the results
  write.table(pvs, file = file.path(outdir, paste0(prefix, "expr_", out_name, "_pvs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  write.table(coeffs, file = file.path(outdir, paste0(prefix, "expr_", out_name, "_coeffs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  
  
  # ----------------------------------------
  # Plot a heatmap of significant cases
  # ----------------------------------------
  
  ### normalize the expression
  expr_norm <- exprc[, c("cluster", "label", "marker", md[md$response != "HD", "shortname"])]
  th <- 2.5
  
  days <- levels(md$day)
  
  ### Normalized to mean = 0 and sd = 1 per day
  for(i in days){
    # i = "base"
    expr_norm[, md[md$response != "HD" & md$day == i, "shortname"]] <- t(apply(expr_norm[, md[md$response != "HD" & md$day == i, "shortname"], drop = FALSE], 1, function(x){ x <- (x-mean(x))/sd(x); x[x > th] <- th; x[x < -th] <- -th; return(x)}))
  }
  
  breaks = seq(from = -th, to = th, length.out = 101)
  legend_breaks = seq(from = -round(th), to = round(th), by = 1)
  
  ### add p-value info
  expr_all <- merge(pvs, expr_norm, by = c("cluster", "label", "marker"), all.x = TRUE, sort = FALSE)
  
  
  
  # -----------------------------
  ### Plot one heatmap with base and tx
  
  ## group the expression by cluster
  if(is.null(adjpval_name)){
    which_top_pvs <- FALSE
  }else{
    expr_all <- expr_all[order(expr_all$label, expr_all[, adjpval_name]), , drop = FALSE]
    
    which_top_pvs <- expr_all[, adjpval_name] < 0.05 & !is.na(expr_all[, adjpval_name])
    which(which_top_pvs)
  }
  
  
  if(sum(which_top_pvs) > 0) {
    print("Plot pheatmap2")
    
    expr_heat <- expr_all[which_top_pvs, , drop = FALSE]
    
    # -----------------------------
    ## order the samples by base and tx
    samples2plot <- md[md$response %in% c("NR", "R"), ]
    samples2plot <- samples2plot$shortname[order(samples2plot$day, samples2plot$response)]
    
    ## gap in the heatmap 
    gaps_col <- c(max(grep("base_NR", samples2plot)), rep(max(grep("base", samples2plot)), 2), max(grep("tx_NR", samples2plot)))
    gaps_row <- unique(cumsum(table(expr_heat$label)))
    gaps_row <- gaps_row[gaps_row > 0]
    if(length(gaps_row) == 1) 
      gaps_row <- NULL
    
    ## expression 
    expr <- expr_heat[ , samples2plot, drop = FALSE]
    
    labels_row <- paste0(expr_heat$label, "/ ", expr_heat$marker, " (", sprintf( "%.02e", expr_heat[, adjpval_name]), ")") 
    labels_col <- colnames(expr)
    
    pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "expr_", out_name, "_", k, "_pheatmap2", suffix, ".pdf")))
    
    
  }
  
  
  # -----------------------------
  ### Plot heatmaps with R vs NR for base and tx
  
  for(i in levels(md$day)){
    # i = "tx"
    
    adjpval_name <- paste0("adjp_NRvsR_", i)
    
    ## group the expression by cluster
    expr_all <- expr_all[order(expr_all$label, expr_all[, adjpval_name]), , drop = FALSE]
    
    which_top_pvs <- expr_all[, adjpval_name] < 0.05 & !is.na(expr_all[, adjpval_name])
    which(which_top_pvs)
    
    if(sum(which_top_pvs) > 0){
      print(paste0("Plot pheatmap_", i))
      
      expr_heat <- expr_all[which_top_pvs, , drop = FALSE]
      
      # -----------------------------
      ## order the samples by NR and R
      
      samples2plot <- md[md$response %in% c("NR", "R"), ]
      samples2plot <- samples2plot$shortname[order(samples2plot$response, samples2plot$day)]
      samples2plot <- samples2plot[grep(i, samples2plot)]
      
      ## gap in the heatmap 
      gaps_col <- sum(grepl("_NR", samples2plot))
      gaps_row <- unique(cumsum(table(expr_heat$label)))
      gaps_row <- gaps_row[gaps_row > 0]
      if(length(gaps_row) == 1) 
        gaps_row <- NULL
      
      ## expression
      expr <- expr_heat[ , samples2plot, drop = FALSE]
      
      labels_row <- paste0(expr_heat$label, "/ ", expr_heat$marker, " (", sprintf( "%.02e", expr_heat[, adjpval_name]), ")") 
      labels_col <- colnames(expr)
      
      pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "expr_", out_name, "_", k, "_pheatmap_", i, suffix, ".pdf")))
      
      
    }
    
  }
  
  # -----------------------------
  ### Plot one heatmap with R vs NR + heatmap with p-values for NRvsR_base, NRvsR_tx and NRvsR_basevstx
  
  ## group the expression by cluster and order by adjpval
  for(i in length(adjpval_name_list):1){
    expr_all <- expr_all[order(expr_all[, adjpval_name_list[i]]), , drop = FALSE]
  }
  
  expr_all <- expr_all[order(expr_all$label), , drop = FALSE]
  
  # which_top_pvs <- rowSums(expr_all[, adjpval_name_list, drop = FALSE] < 0.05, na.rm = TRUE) > 0 & rowSums(is.na(expr_all[, adjpval_name_list, drop = FALSE])) < length(adjpval_name_list)
  which_top_pvs <- rep(TRUE, nrow(expr_all))
  which(which_top_pvs)
  
  if(sum(which_top_pvs) > 0) {
    print("Plot pheatmap3")
    
    expr_heat <- expr_all[which_top_pvs, , drop = FALSE]
    
    # -----------------------------
    ## order the samples by base and tx
    samples2plot <- md[md$response %in% c("NR", "R"), ]
    samples2plot <- samples2plot$shortname[order(samples2plot$day, samples2plot$response)]
    
    ## gap in the heatmap 
    gaps_col <- c(max(grep("base_NR", samples2plot)), rep(max(grep("base", samples2plot)), 2), max(grep("tx_NR", samples2plot)))
    gaps_row <- unique(cumsum(table(expr_heat$label)))
    gaps_row <- gaps_row[gaps_row > 0]
    if(length(gaps_row) == 1) 
      gaps_row <- NULL
    
    ## expression 
    expr <- expr_heat[ , samples2plot, drop = FALSE]
    
    labels_row <- paste0(expr_heat$label, "/ ", expr_heat$marker) 
    labels_col <- colnames(expr)
    
    pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "expr_", out_name, "_", k, "_pheatmap3", suffix, ".pdf")))
    
    pvs_heat <- expr_heat[, adjpval_name_list, drop = FALSE]
    
    labels_col <- colnames(pvs_heat)
    gaps_col <- NULL
    
    pheatmap(pvs_heat, cellwidth = 60, cellheight = 24, color = c("grey50", "grey70", "grey90"), breaks = c(0, 0.05, 0.1, 1), legend_breaks = c(0, 0.05, 0.1, 1), legend = FALSE, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, display_numbers = TRUE, number_format = "%.02e", number_color = "black", fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "expr_", out_name, "_", k, "_pheatmap3pvs", suffix, ".pdf")))
    
    
  }
  
  
  
  # ----------------------------------------
  # Plot coefficients NRvsR for base and tx (to show that they correlate)
  # ----------------------------------------
  
  
  if("adjp_NRvsR_basevstx" %in% colnames(pvs)){
    
    adjpval_name <- "adjp_NRvsR_basevstx"
    
    ggdf <- coeffs[, c("NRvsR_base", "NRvsR_tx")]
    
    limmin <- min(ggdf, na.rm = TRUE)
    limmax <- max(ggdf, na.rm = TRUE)
    
    ggdf$interaction <- factor(pvs[, adjpval_name] < 0.05, levels = c("FALSE", "TRUE"))
    
    ggp <- ggplot(data = ggdf, aes(x = NRvsR_base, y = NRvsR_tx, shape = interaction)) +
      geom_point(size = 3, alpha = 0.75) +
      geom_abline(intercept = 0, slope = 1) +
      coord_cartesian(xlim = c(limmin, limmax), ylim = c(limmin, limmax)) +
      theme_bw() +
      theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face="bold"))
    
    pdf(file.path(outdir, paste0(prefix, "expr_", out_name, "_", k, "_coeffs", suffix, ".pdf")), w=6, h=5, onefile=TRUE)
    print(ggp)
    dev.off()
    
    
  }
  
  
  
  
} # models2fit







sessionInfo()














################################
### 04_expression done!
################################