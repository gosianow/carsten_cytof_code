

Sys.time()

# Load packages
library(gdata)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2)
library(RColorBrewer)
library(limma)
library(pheatmap)
library(ComplexHeatmap)

##############################################################################
# Test arguments
##############################################################################


prefix='23_01_pca1_cl20_all_'
outdir='../carsten_cytof/PD1_project/CK_2016-06-23_01/080_expression_auto'
path_metadata='../carsten_cytof/PD1_project/CK_metadata/metadata_23_01.xlsx'
path_expression='../carsten_cytof/PD1_project/CK_2016-06-23_01/080_expression_auto/23_01_pca1_cl20_all_expr.xls'
path_fun_models='00_models.R'
path_fun_formulas='00_formulas_1dataset_3responses.R'
path_fun_plot_heatmaps="00_plot_heatmaps_for_sign_expr.R"
path_marker_exclusion='../carsten_cytof/PD1_project/CK_2016-06-23_01/010_helpfiles/23_01_pca1_cl20_all_marker_exclusion.txt'
FDR_cutoff='05'


### Can be NULL
# path_marker_exclusion


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

if( !file.exists(outdir) ) 
  dir.create(outdir, recursive = TRUE)

suffix <- paste0("_top", FDR_cutoff)
FDR_cutoff <- as.numeric(paste0("0.", FDR_cutoff))
FDR_cutoff

# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

md <- lapply(1:length(path_metadata), function(i){
  md <- read.xls(path_metadata[i], stringsAsFactors = FALSE)
  md
})

md <- plyr::rbind.fill(md)
rownames(md) <- md$shortname

### Factor arrangment
md$response <- factor(md$response, levels = c("NR", "R", "HD"))
md$response <- factor(md$response)
md$day <- factor(md$day, levels = c("base", "tx"))
md$day <- factor(md$day)
md$patient_id <- factor(md$patient_id)
md$data <- factor(md$data)
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
# Load median expression data
# ------------------------------------------------------------


a <- lapply(1:length(path_expression), function(i){
  # i = 1
  a <- read.table(path_expression[i], header = TRUE, sep = "\t", as.is = TRUE, check.names = FALSE)
  a <- a[, -which(colnames(a) == "cluster")]
  a
})


## keep only those markers that are common for all the merged datasets 
marker_list <- unlist(lapply(a, colnames))
marker_list <- marker_list[!grepl("label|sample", marker_list)]

marker_table <- table(marker_list)
marker_overlap <- names(marker_table)[marker_table == length(path_expression)]

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
head(a)

markers_ordered <- colnames(a)[!colnames(a) %in% c("cluster", "label", "sample")]
markers_ordered


### Skipp samples that have NA for all the clusters
### Clusters with at least one NA will be skipped in the expression analysis (see 00_models.R)

table_sample <- table(a[complete.cases(a), "sample"])

keep_samps <- names(table_sample)[table_sample > 0]

a <- a[a$sample %in% keep_samps, , drop = FALSE]

md <- md[md$shortname %in% keep_samps, , drop = FALSE]


# Drop unused levels
md$response <- factor(md$response)
md$day <- factor(md$day)
md$data <- factor(md$data)
md$patient_id <- factor(md$patient_id)



# ------------------------------------------------------------
# Load marker exclusion for plotting on the heatmaps
# ------------------------------------------------------------

marker_exclusion <- NULL

if(!is.null(path_marker_exclusion)){
  if(file.exists(path_marker_exclusion)){
    
    marker_exclusion <- read.table(path_marker_exclusion, header = TRUE, sep = "\t", as.is = TRUE)
    marker_exclusion <- marker_exclusion[, 1]
    
    if(!all(marker_exclusion %in% markers_ordered))
      stop("Marker exclusion is wrong")
    
  }
}


# -----------------------------------------------------------------------------
# Normalize the expression
# -----------------------------------------------------------------------------
### Prepare the matrix with data (rows - markers x clusters; columns - samples)

expr <- a
exprm <- melt(expr, id.vars = c("cluster", "label", "sample"), value.name = "expr", variable.name = "marker")
exprc <- dcast(exprm, cluster + label + marker ~ sample, value.var = "expr")

expr_norm <- exprc[, c("cluster", "label", "marker", md[md$response != "HD", "shortname"])]
th <- 2.5

md$data_day <- factor(md$data_day)
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

breaks = seq(from = -th, to = th, length.out = 101)
legend_breaks = seq(from = -round(th), to = round(th), by = 1)



# -----------------------------------------------------------------------------
# Plot a heatmap with clustered columns and all the rows
# -----------------------------------------------------------------------------


if(all(a$label == "all")){
  
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
  
  pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = cluster_cols, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_col = annotation_col, annotation_colors = annotation_colors, annotation_legend = TRUE, filename = file.path(outdir, paste0(prefix, "expr_pheatmap_colclust_both_asis.pdf")))
  
  
  ### Using ComplexHeatmap
  
  ha <-  HeatmapAnnotation(df = annotation_col, col = list(response = color_response[levels(annotation_col$response)]))
  
  ht1 <- Heatmap(expr, name = "", col = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), cluster_columns = cluster_cols, cluster_rows = cluster_rows, column_dend_reorder = FALSE, row_dend_reorder = FALSE, top_annotation = ha, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"))
  
  pdf(file.path(outdir, paste0(prefix, "expr_ComplexHeatmap_colclust_both_asis.pdf")), width = 10, height = 7)
  try(draw(ht1), silent = TRUE)
  dev.off()
  
  
  ### Reorder the dendrogram
  ht1 <- Heatmap(expr, name = "", col = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), cluster_columns = cluster_cols, cluster_rows = cluster_rows, column_dend_reorder = TRUE, row_dend_reorder = FALSE, top_annotation = ha, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"))
  
  pdf(file.path(outdir, paste0(prefix, "expr_ComplexHeatmap_colclust_both_reord.pdf")), width = 10, height = 7)
  try(draw(ht1), silent = TRUE)
  dev.off()
  
  
  ## Plot only those markers that are not excluded
  if(!is.null(marker_exclusion)){
    
    expr_sub <- expr[!rownames(expr) %in% marker_exclusion, , drop = FALSE]
    labels_row_sub <- rownames(expr_sub)
    
    cluster_cols <- hclust(dist(t(expr_sub)), method = "ward.D2")
    cluster_rows <- hclust(dist(expr_sub), method = "ward.D2")
    
    # Using pheatmap
    pheatmap(expr_sub, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = cluster_cols, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row_sub, fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_col = annotation_col, annotation_colors = annotation_colors, annotation_legend = TRUE, filename = file.path(outdir, paste0(prefix, "expr_pheatmap_colclust_both_asis_ex.pdf")))
    
    
    ### Using ComplexHeatmap
    ha <-  HeatmapAnnotation(df = annotation_col, col = list(response = color_response[levels(annotation_col$response)]))
    ht1 <- Heatmap(expr_sub, name = "", col = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), cluster_columns = cluster_cols, cluster_rows = cluster_rows, column_dend_reorder = FALSE, row_dend_reorder = FALSE, top_annotation = ha, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"))
    
    pdf(file.path(outdir, paste0(prefix, "expr_ComplexHeatmap_colclust_both_asis_ex.pdf")), width = 10, height = 7)
    try(draw(ht1), silent = TRUE)
    dev.off()
    
    ### Reorder the dendrogram
    ht1 <- Heatmap(expr_sub, name = "", col = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), cluster_columns = cluster_cols, cluster_rows = cluster_rows, column_dend_reorder = TRUE, row_dend_reorder = FALSE, top_annotation = ha, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"))
    
    pdf(file.path(outdir, paste0(prefix, "expr_ComplexHeatmap_colclust_both_reord_ex.pdf")), width = 10, height = 7)
    try(draw(ht1), silent = TRUE)
    dev.off()
    
  }
  
  ### ----------------------------------------------------------------
  ### Plot using only the base or tx samples
  ### ----------------------------------------------------------------
  
  for(i in c("base", "tx")){
    
    expr_subday <- expr_heat[, md[md$response != "HD" & md$day == i, "shortname"]]
    labels_col_base <- colnames(expr_subday)
    
    annotation_col <- data.frame(response = factor(md[md$response != "HD" & md$day == i, "response"]))
    rownames(annotation_col) <- md[md$response != "HD" & md$day == i, "shortname"]
    
    annotation_colors <- list(response = color_response[levels(annotation_col$response)])
    
    cluster_cols <- hclust(dist(t(expr_subday)), method = "ward.D2")
    cluster_rows <- hclust(dist(expr_subday), method = "ward.D2")
    
    # Using pheatmap
    
    pheatmap(expr_subday, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = cluster_cols, cluster_rows = cluster_rows, labels_col = labels_col_base, labels_row = labels_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_col = annotation_col, annotation_colors = annotation_colors, annotation_legend = TRUE, filename = file.path(outdir, paste0(prefix, "expr_pheatmap_colclust_", i, "_asis.pdf")))
    
    
    ### Using ComplexHeatmap
    
    ha <-  HeatmapAnnotation(df = annotation_col, col = list(response = color_response[levels(annotation_col$response)]))
    
    ht1 <- Heatmap(expr_subday, name = "", col = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), cluster_columns = cluster_cols, cluster_rows = cluster_rows, column_dend_reorder = FALSE, row_dend_reorder = FALSE, top_annotation = ha, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"))
    
    pdf(file.path(outdir, paste0(prefix, "expr_ComplexHeatmap_colclust_", i, "_asis.pdf")), width = 7, height = 7)
    try(draw(ht1), silent = TRUE)
    dev.off()
    
    ### Reorder the dendrogram
    ht1 <- Heatmap(expr_subday, name = "", col = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), cluster_columns = cluster_cols, cluster_rows = cluster_rows, column_dend_reorder = TRUE, row_dend_reorder = FALSE, top_annotation = ha, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"))
    
    pdf(file.path(outdir, paste0(prefix, "expr_ComplexHeatmap_colclust_", i, "_reord.pdf")), width = 7, height = 7)
    try(draw(ht1), silent = TRUE)
    dev.off()
    
    ## Plot only those markers that are not excluded
    if(!is.null(marker_exclusion)){
      
      expr_sub <- expr_subday[!rownames(expr_subday) %in% marker_exclusion, , drop = FALSE]
      labels_row_sub <- rownames(expr_sub)
      
      cluster_cols <- hclust(dist(t(expr_sub)), method = "ward.D2")
      cluster_rows <- hclust(dist(expr_sub), method = "ward.D2")
      
      # Using pheatmap
      pheatmap(expr_sub, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = cluster_cols, cluster_rows = cluster_rows, labels_col = labels_col_base, labels_row = labels_row_sub, fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_col = annotation_col, annotation_colors = annotation_colors, annotation_legend = TRUE, filename = file.path(outdir, paste0(prefix, "expr_pheatmap_colclust_", i, "_asis_ex.pdf")))
      
      
      ### Using ComplexHeatmap
      ha <-  HeatmapAnnotation(df = annotation_col, col = list(response = color_response[levels(annotation_col$response)]))
      ht1 <- Heatmap(expr_sub, name = "", col = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), cluster_columns = cluster_cols, cluster_rows = cluster_rows, column_dend_reorder = FALSE, row_dend_reorder = FALSE, top_annotation = ha, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"))
      
      pdf(file.path(outdir, paste0(prefix, "expr_ComplexHeatmap_colclust_", i, "_asis_ex.pdf")), width = 7, height = 7)
      try(draw(ht1), silent = TRUE)
      dev.off()
      
      ### Reorder the dendrogram
      ht1 <- Heatmap(expr_sub, name = "", col = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), cluster_columns = cluster_cols, cluster_rows = cluster_rows, column_dend_reorder = TRUE, row_dend_reorder = FALSE, top_annotation = ha, heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"))
      
      pdf(file.path(outdir, paste0(prefix, "expr_ComplexHeatmap_colclust_", i, "_reord_ex.pdf")), width = 7, height = 7)
      try(draw(ht1), silent = TRUE)
      dev.off()
      
    }
    
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


# models2fit <- c("lm_interglht", "lmer_interglht")
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
  write.table(pvs, file = file.path(outdir, paste0(prefix, "expr_pvs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  write.table(coeffs, file = file.path(outdir, paste0(prefix, "expr_coeffs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  
  
  # ----------------------------------------
  # Plot a heatmap of significant cases
  # ----------------------------------------
  
  ### add p-value info
  expr_all <- merge(pvs, expr_norm, by = c("cluster", "label", "marker"), all.x = TRUE, sort = FALSE)
  
  prefix2 <- paste0(k, "_")
  
  plot_heatmaps_for_sign_expr(expr_all = expr_all, md = md, FDR_cutoff = FDR_cutoff, pval_name2 = pval_name2, adjpval_name2 = adjpval_name2, pval_name_list = pval_name_list, adjpval_name_list = adjpval_name_list, breaks = breaks, legend_breaks = legend_breaks, color_response = color_response, outdir = outdir, prefix = prefix, prefix2 = prefix2, suffix = suffix)
  
  
  ## Plot only those markers that are not excluded
  if(!is.null(marker_exclusion)){
    
    expr_all <- expr_all[!expr_all$marker %in% marker_exclusion, , drop = FALSE]
    
    plot_heatmaps_for_sign_expr(expr_all = expr_all, md = md, FDR_cutoff = FDR_cutoff, pval_name2 = pval_name2, adjpval_name2 = adjpval_name2, pval_name_list = pval_name_list, adjpval_name_list = adjpval_name_list, breaks = breaks, legend_breaks = legend_breaks, color_response = color_response, outdir = outdir, prefix = prefix, prefix2 = prefix2, suffix = paste0(suffix, "_ex"))
    
  }
  
  
}











sessionInfo()





