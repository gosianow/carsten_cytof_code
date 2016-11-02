##############################################################################
## <<08_expression_merged.R>>

# BioC 3.3
# Created 13 Oct 2016
# Updated 20 Oct 2016

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
library(tools)

##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-merged_23_29/01'
expr_prefix='23m6_29m4_'
expr_outdir='08_expression_merged'

path_metadata=c('/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.xlsx','/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_29_01.xlsx')
path_expression=c('/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01/080_expression/23_01_pca1_merging6_raw_expr_all.xls','/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-29_01/080_expression/29_01_pca1_merging4_raw_expr_all.xls')
data_name=c('data23','data29')

path_fun_models='/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_models.R'
path_fun_formulas='/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_formulas_2datasets_3responses.R'

analysis_type='all'

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


if(!analysis_type %in% c("clust", "all"))
  stop("analysis_type must be 'all' or 'clust'!")

out_name <- analysis_type


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

color_samples <- md$color
names(color_samples) <- md$shortname

colors <- unique(md[, c("response", "color")])
color_response <- colors$color
names(color_response) <- colors$response


# ------------------------------------------------------------
# Load expression data
# ------------------------------------------------------------

a <- lapply(1:length(data_name), function(i){
  # i = 1
  path <- path_expression[i]
  a <- read.table(path, header = TRUE, sep = "\t", as.is = TRUE, check.names = FALSE)
  a[, -which(colnames(a) == "cluster")]
  
})

a <- rbind.fill(a)

## keep only those marksers that are common for all the merged datasets 
a <- a[, complete.cases(t(a)), drop = FALSE]

a <- a[a$label != "drop", , drop = FALSE]


## keep only these clusters that are present in all the merged datasets
nr_samples <- length(unique(a$sample))

labels_keep <- names(which(table(a$label) == nr_samples))

labels <- unique(a$label)
labels <- labels[labels %in% labels_keep]

labels <- data.frame(cluster = 1:length(labels), label = labels)
labels$label <- factor(labels$label, levels = unique(labels$label))


a <- a[a$label %in% labels_keep, , drop = FALSE]

mm <- match(a$label, labels$label)

a <- cbind(cluster = labels$cluster[mm], a)

markers_ordered <- colnames(a)[!colnames(a) %in% c("cluster", "label", "sample")]


# -----------------------------------------------------------------------------
### Plot expression per cluster
# -----------------------------------------------------------------------------


ggdf <- melt(a, id.vars = c("cluster", "label", "sample"), value.name = "expr", variable.name = "marker")

## use labels as clusters
ggdf$cluster <- factor(ggdf$label, levels = labels$label)


## add group info
mm <- match(ggdf$samp, md$shortname)
ggdf$group <- factor(md$condition[mm])
## add data info
ggdf$data <- factor(md$data[mm], levels = data_name)


## replace _ with \n
levels(ggdf$group) <- gsub("_", "\n", levels(ggdf$group))

## order markers as for heatmaps
ggdf$marker <- factor(ggdf$marker, levels = colnames(a)[!grepl("cluster|label|sample", colnames(a))])

## calculate mean and sd for the error bars on the plot
ggds <- ddply(ggdf, .(group, cluster, marker, data), summarise, mean = mean(expr), sd = sd(expr))

clusters <- levels(ggdf$cluster)

# ------------------------------------
## plot each cluster as a separate page in the pdf file

ggp <- list()

for(i in 1:nlevels(ggdf$cluster)){
  # i = 1
  
  df <- ggdf[ggdf$cluster == clusters[i], , drop = FALSE]
  ds <- ggds[ggds$cluster == clusters[i], , drop = FALSE]
  
  ggp[[i]] <- ggplot(df, aes(x = group, y = expr, color = group, shape = data)) +
    geom_point(size=2.5, position = position_jitterdodge(jitter.width = 3, jitter.height = 0, dodge.width = 0.7)) +
    geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean, ymax=mean), color='black', width=0.4, position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.7)) +
    geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean-sd, ymax=mean+sd), color='black', width=0.25, position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.7)) +
    facet_wrap(~ marker, scales = "free") +
    ggtitle(clusters[i]) +
    theme_bw() +
    ylab("Expression") +
    xlab("") +
    theme(axis.text.x = element_text(size=10, face="bold"),
      axis.title.y = element_text(size=12, face="bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
      axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
      legend.position = "right") +
    scale_color_manual(values = color_groups)
  
}

pdf(file.path(outdir, paste0(prefix, "expr_", out_name, ".pdf")), w = 18, h = 12, onefile=TRUE)
for(i in seq(length(ggp)))
  print(ggp[[i]])
dev.off()


# -----------------------------------------------------------------------------
# Prepare the matrix with data (rows - markers X clusters; columns - samples)
# Plot a heatmap with clustered samples
# -----------------------------------------------------------------------------


expr <- a
exprm <- melt(expr, id.vars = c("cluster", "label", "sample"), value.name = "expr", variable.name = "marker")
exprc <- dcast(exprm, cluster + label + marker ~ sample, value.var = "expr")


### normalize the expression
expr_norm <- exprc[, c("cluster", "label", "marker", md[md$response != "HD", "shortname"])]
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

# ### Normalize to mean = 0 per data and day
# for(i in data_day){
#   # i = "data23.base"
#   expr_norm[, md[md$response != "HD" & md$data_day == i, "shortname"]] <- t(apply(expr_norm[, md[md$response != "HD" & md$data_day == i, "shortname"], drop = FALSE], 1, function(x){ x <- x-mean(x); return(x)}))
# }
# ### Normalize to sd = 1 for all samples
# expr_norm[, md[md$response != "HD", "shortname"]] <- t(apply(expr_norm[, md[md$response != "HD", "shortname"], drop = FALSE], 1, function(x){ x <- x/sd(x); x[x > th] <- th; x[x < -th] <- -th; return(x)}))


breaks = seq(from = -th, to = th, length.out = 101)
legend_breaks = seq(from = -round(th), to = round(th), by = 1)

if(analysis_type == "all"){
  
  expr_heat <- expr_norm
  rownames(expr_heat) <- expr_heat$marker
  
  expr_heat <- expr_heat[markers_ordered, ]
  expr <- expr_heat[, md[md$response != "HD", "shortname"]]
  
  labels_row <- paste0(expr_heat$label, "/ ", expr_heat$marker) 
  labels_col <- colnames(expr)
  
  annotation_col <- data.frame(response = factor(md[md$response != "HD", "response"]))
  rownames(annotation_col) <- md[md$response != "HD", "shortname"]
  
  annotation_colors <- list(response = color_response[levels(annotation_col$response)])
  
  cluster_cols <- hclust(dist(t(expr)), method = "ward.D2")
  cluster_rows <- hclust(dist(expr), method = "ward.D2")
  
  pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = cluster_cols, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_col = annotation_col, annotation_colors = annotation_colors, filename = file.path(outdir, paste0(prefix, "expr_", out_name,  "_pheatmap_colclust", suffix, ".pdf")))
  
}




# -----------------------------------------------------------------------------
# Test for marker expression differences between groups overall and per cluster 
# -----------------------------------------------------------------------------

### Load functions fitting models
source(path_fun_models)
### Load formulas that are fit in the models - this function may change the md object!!!
source(path_fun_formulas)

levels(md$data)
levels(md$day)
levels(md$response)


### Fit all the models

models2fit <- c("lm_interglht", "lmer_interglht", "rlm_interglht")

for(k in models2fit){
  # k = "lmer_interglht"
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
  
  # table(pvs$adjp_NRvsR < 0.05, useNA = "always")
  # table(pvs$adjp_NRvsR_base < 0.05, useNA = "always")
  # table(pvs$adjp_NRvsR_tx < 0.05, useNA = "always")
  # table(pvs$adjp_NRvsR_basevstx < 0.05, useNA = "always")
  #
  # table(pvs$pval_NRvsR < 0.05, useNA = "always")
  # table(pvs$pval_NRvsR_base < 0.05, useNA = "always")
  # table(pvs$pval_NRvsR_tx < 0.05, useNA = "always")
  # table(pvs$pval_NRvsR_basevstx < 0.05, useNA = "always")
  
  # ----------------------------------------
  # Plot a heatmap of significant cases
  # ----------------------------------------
  
  ### add p-value info
  expr_all <- merge(pvs, expr_norm, by = c("cluster", "label", "marker"), all.x = TRUE, sort = FALSE)
  
  # -----------------------------
  ### Plot one heatmap with base and tx
  
  ## group the expression by cluster
  adjpval_name <- adjpval_name2
  expr_all <- expr_all[order(expr_all$label, expr_all[, adjpval_name]), , drop = FALSE]
  
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
  ### Plot two heatmaps with R vs NR for base and tx
  
  for(i in levels(md$day)){
    # i = "tx"
    print(paste0("Plot pheatmap_", i))
    
    adjpval_name <- paste0("adjp_NRvsR_", i)
    
    ## group the expression by cluster
    expr_all <- expr_all[order(expr_all$label, expr_all[, adjpval_name]), , drop = FALSE]
    
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
  
  which_top_pvs <- rowSums(expr_all[, adjpval_name_list, drop = FALSE] < 0.05, na.rm = TRUE) > 0 & rowSums(is.na(expr_all[, adjpval_name_list, drop = FALSE])) < length(adjpval_name_list)
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
### 08_expression_merged.R done!
################################