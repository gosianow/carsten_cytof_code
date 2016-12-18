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
library(ComplexHeatmap)


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
path_fun_plot_heatmaps <- "/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_plot_heatmaps_for_sign_expr.R"
path_fun_plot_expression <- "/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_plot_expression.R"
path_marker_exclusion='23m4_29m2_expr_marker_exclusion.txt'
analysis_type='all'

### Optional arguments
FDR_cutoff=0.05
suffix="_top005"

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

prefix <- expr_prefix
outdir <- expr_outdir


if( !file.exists(outdir) ) 
  dir.create(outdir)


if(!analysis_type %in% c("clust", "all"))
  stop("analysis_type must be 'all' or 'clust'!")

out_name <- analysis_type


if(!any(grepl("FDR_cutoff=", args))){
  FDR_cutoff=0.05
}

if(!any(grepl("suffix=", args))){
  suffix="_top005"
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



# ------------------------------------------------------------
# Load expression data
# ------------------------------------------------------------

a <- lapply(1:length(data_name), function(i){
  # i = 2
  path <- path_expression[i]
  a <- read.table(path, header = TRUE, sep = "\t", as.is = TRUE, check.names = FALSE)
  a <- a[, -which(colnames(a) == "cluster")]
  a
})


## keep only those markers that are common for all the merged datasets 
marker_list <- unlist(lapply(a, colnames))
marker_list <- marker_list[!grepl("label|sample", marker_list)]

marker_table <- table(marker_list)
marker_overlap <- names(marker_table)[marker_table == length(data_name)]

a <- rbind.fill(a)

a <- a[, colnames(a) %in% c(marker_overlap, "label", "sample"), drop = FALSE]

## drop the "drop" cluster
a <- a[a$label != "drop", , drop = FALSE]
a <- a[a$sample %in% md$shortname, , drop = FALSE]

## keep only those clusters that are present in all the merged datasets
labels_keep <- names(which(table(a$label) == nrow(md)))

labels <- unique(a$label)
labels <- labels[labels %in% labels_keep]

labels <- data.frame(cluster = 1:length(labels), label = labels)
labels$label <- factor(labels$label, levels = unique(labels$label))
labels

a <- a[a$label %in% labels_keep, , drop = FALSE]

mm <- match(a$label, labels$label)

a <- cbind(cluster = labels$cluster[mm], a)

markers_ordered <- colnames(a)[!colnames(a) %in% c("cluster", "label", "sample")]

### Save the median expression per cluster and sample
write.table(a, file.path(outdir, paste0(prefix, "expr_clust_", out_name, ".xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# ---------------------------------------
# Keep only those samples that have enough cells in any of the clusters
# Clusters with at least one NA will be skipped in the expression analysis (see 00_models.R)
# ---------------------------------------

table_sample <- table(a[complete.cases(a), "sample"])

keep_samps <- names(table_sample)[table_sample > 0]

a <- a[a$sample %in% keep_samps, , drop = FALSE]

md <- md[md$shortname %in% keep_samps, , drop = FALSE]

## drop unused levels
md$response <- factor(md$response)
md$day <- factor(md$day)
md$patient_id <- factor(md$patient_id)


# ------------------------------------------------------------
# load marker exclusion for plotting on the heatmaps
# ------------------------------------------------------------

marker_exclusion <- NULL

if(file.exists(file.path(path_marker_exclusion))){
  
  marker_exclusion <- read.table(file.path(path_marker_exclusion), header = TRUE, sep = "\t", as.is = TRUE)
  marker_exclusion <- marker_exclusion[, 1]
  
  if(!all(marker_exclusion %in% markers_ordered))
    stop("Marker exclusion is wrong")
  
}


# -----------------------------------------------------------------------------
### Plot expression per cluster
# -----------------------------------------------------------------------------

source(path_fun_plot_expression)


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


plot_expression(ggdf = ggdf, ggds = ggds, color_groups = color_groups, outdir = outdir, prefix = prefix, prefix2 = out_name)


# -----------------------------------------------------------------------------
# Prepare the matrix with data for heatmaps (rows - markers X clusters; columns - samples)
# Plot a heatmap with clustered columns
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
  print(i)
  
  expr_norm[, md[md$response != "HD" & md$data_day == i, "shortname"]] <- t(apply(expr_norm[, md[md$response != "HD" & md$data_day == i, "shortname"], drop = FALSE], 1, function(x){ 
    
    if(sum(!is.na(x)) == 0)
      return(x)
    
    if(sum(!is.na(x)) < 2)
      return(x-mean(x, na.rm = TRUE))
    
    sdx <- sd(x, na.rm = TRUE)
    if(sdx == 0)
      x <- x-mean(x, na.rm = TRUE)
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


### Plot a heatmap with clustered columns and all the rows
if(analysis_type == "all"){
  
  expr_heat <- expr_norm
  rownames(expr_heat) <- expr_heat$marker
  
  expr_heat <- expr_heat[markers_ordered, ]
  expr <- expr_heat[, md[md$response != "HD", "shortname"]]
  
  labels_row <- paste0(expr_heat$marker) 
  labels_col <- colnames(expr)
  
  annotation_col <- data.frame(response = factor(md[md$response != "HD", "response"]))
  rownames(annotation_col) <- md[md$response != "HD", "shortname"]
  
  annotation_colors <- list(response = color_response[levels(annotation_col$response)])
  
  cluster_cols <- hclust(dist(t(expr)), method = "ward.D2")
  cluster_rows <- hclust(dist(expr), method = "ward.D2")
  
  # Using pheatmap
  
  pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = cluster_cols, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_col = annotation_col, annotation_colors = annotation_colors, annotation_legend = TRUE, filename = file.path(outdir, paste0(prefix, "expr_", out_name,  "_pheatmap_colclust", ".pdf")))
  
  
  ### Using ComplexHeatmap
  
  ha <-  HeatmapAnnotation(df = annotation_col, col = list(response = color_response[levels(annotation_col$response)]))
  
  ht1 <- Heatmap(expr, name = "", col = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), cluster_columns = cluster_cols, cluster_rows = cluster_rows, column_dend_reorder = FALSE, row_dend_reorder = FALSE, top_annotation = ha, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"))
  
  pdf(file.path(outdir, paste0(prefix, "expr_", out_name,  "_ComplexHeatmap_colclust", ".pdf")), width = 10, height = 7)
  draw(ht1)
  dev.off()
  
  ## Plot only those markers that are not excluded
  if(!is.null(marker_exclusion)){
    
    expr_sub <- expr[!rownames(expr) %in% marker_exclusion, , drop = FALSE]
    labels_row_sub <- rownames(expr_sub)
    
    cluster_cols <- hclust(dist(t(expr_sub)), method = "ward.D2")
    cluster_rows <- hclust(dist(expr_sub), method = "ward.D2")
    
    # Using pheatmap
    pheatmap(expr_sub, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = cluster_cols, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row_sub, fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_col = annotation_col, annotation_colors = annotation_colors, annotation_legend = TRUE, filename = file.path(outdir, paste0(prefix, "expr_", out_name,  "_pheatmap_colclust_ex", ".pdf")))
    
    
    ### Using ComplexHeatmap
    ha <-  HeatmapAnnotation(df = annotation_col, col = list(response = color_response[levels(annotation_col$response)]))
    ht1 <- Heatmap(expr_sub, name = "", col = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), cluster_columns = cluster_cols, cluster_rows = cluster_rows, column_dend_reorder = FALSE, row_dend_reorder = FALSE, top_annotation = ha, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"))
    
    pdf(file.path(outdir, paste0(prefix, "expr_", out_name,  "_ComplexHeatmap_colclust_ex", ".pdf")), width = 10, height = 7)
    draw(ht1)
    dev.off()
    
    
  }
  
}




# -----------------------------------------------------------------------------
# Test for marker expression differences between groups overall and per cluster
# -----------------------------------------------------------------------------
## The model functions do not anlyse a cluster with NAs; 
## For merged data it means such cluster was not present in all the datasets
## For expression data clusters with no cells are skipped


### Load functions fitting models
source(path_fun_models)
### Load formulas that are fit in the models - this function may change the md object!!!
source(path_fun_formulas)

source(path_fun_plot_heatmaps)


levels(md$data)
levels(md$day)
levels(md$response)


### Fit all the models

# models2fit <- c("lm_interglht", "lmer_interglht", "rlm_interglht")
models2fit <- c("lmer_interglht")


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
  
  oo <- order(pvs[, pval_name1], decreasing = FALSE)
  pvs <- pvs[oo, , drop = FALSE]
  coeffs <- coeffs[oo, , drop = FALSE]
  
  ## save the results
  write.table(pvs, file = file.path(outdir, paste0(prefix, "expr_", out_name, "_pvs_", k, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  write.table(coeffs, file = file.path(outdir, paste0(prefix, "expr_", out_name, "_coeffs_", k, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  
  
  # ----------------------------------------
  # Plot a heatmap of significant cases
  # ----------------------------------------
  
  ### add p-value info
  expr_all <- merge(pvs, expr_norm, by = c("cluster", "label", "marker"), all.x = TRUE, sort = FALSE)
  
  prefix2 <- paste0(out_name, "_", k, "_")
  
  plot_heatmaps_for_sign_expr(expr_all = expr_all, md = md, FDR_cutoff = FDR_cutoff, pval_name2 = pval_name2, adjpval_name2 = adjpval_name2, pval_name_list = pval_name_list, adjpval_name_list = adjpval_name_list, breaks = breaks, legend_breaks = legend_breaks, outdir = outdir, prefix = prefix, prefix2 = prefix2, suffix = suffix)
  
  
  ## Plot only those markers that are not excluded
  if(!is.null(marker_exclusion)){
    
    expr_all <- expr_all[!expr_all$marker %in% marker_exclusion, , drop = FALSE]
    
    plot_heatmaps_for_sign_expr(expr_all = expr_all, md = md, FDR_cutoff = FDR_cutoff, pval_name2 = pval_name2, adjpval_name2 = adjpval_name2, pval_name_list = pval_name_list, adjpval_name_list = adjpval_name_list, breaks = breaks, legend_breaks = legend_breaks, outdir = outdir, prefix = prefix, prefix2 = prefix2, suffix = paste0(suffix, "_ex"))
    
  }
  
  
  
} # models2fit







sessionInfo()














################################
### 08_expression_merged.R done!
################################