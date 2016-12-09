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

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
expr_prefix='23_01_pca1_merging6_raw_'
expr_outdir='080_expression'
path_data='010_data/23_01_expr_raw.rds'
path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.xlsx'
path_clustering_observables='030_heatmaps/23_01_pca1_clustering_observables.xls'
path_clustering='030_heatmaps/23_01_pca1_merging6_clustering.xls'
path_clustering_labels='030_heatmaps/23_01_pca1_merging6_clustering_labels.xls'
path_fun_models='/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_models.R'
path_fun_formulas='/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_formulas_1dataset_3responses.R'
analysis_type='all'

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-29_03'
expr_prefix='29_03_pca1_merging2_raw_'
expr_outdir='080_expression'
path_data='010_data/29_03_expr_raw.rds'
path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_29_03.xlsx'
path_clustering_observables='030_heatmaps/29_03_pca1_clustering_observables.xls'
path_clustering='030_heatmaps/29_03_pca1_merging2_clustering.xls'
path_clustering_labels='030_heatmaps/29_03_pca1_merging2_clustering_labels.xls'
path_fun_models='/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_models.R'
path_fun_formulas='/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_formulas_1dataset_3responses.R'
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

path_fun_plot_heatmaps <- "/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_plot_heatmaps_for_sign_expr.R"
source(path_fun_plot_heatmaps)

setwd(rwd)

prefix <- expr_prefix
outdir <- expr_outdir
suffix <- ""

if( !file.exists(outdir) ) 
  dir.create(outdir, recursive = TRUE)


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
# Get the median expression per cluster, and if sample has not enough cells, set expression to NA
# ------------------------------------------------------------


if(analysis_type == "clust"){
  
  min_cells <- 20
  
  table_samp <- aggregate(e[, 1, drop = FALSE], by = list(clust, samp), FUN = length, drop = FALSE)
  
  keep_samps <- table_samp[, 3] > min_cells
  
  colnames(e) <- fcs_panel$Antigen
  
  a <- aggregate(e, by = list(clust, samp), FUN = median, drop = FALSE)
  
  mlab <- match(a$Group.1, labels$cluster)
  a$label <- labels$label[mlab]
  
  colnames(a)[1:2] <- c("cluster", "sample")
  
  a <- a[, c("cluster", "label", "sample", fcs_panel$Antigen[c(scols, xcols)])]
  a$label <- factor(a$label, levels = labels$label)
  
  a[!keep_samps, fcs_panel$Antigen[c(scols, xcols)]] <- NA
  
  ### Save the median expression per cluster and sample
  write.table(a, file.path(outdir, paste0(prefix, "expr_clust.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  ### Skipp samples that have NA for all the clusters
  ### Clusters with at least one NA will be skipped in the expression analysis (see 00_models.R)
  
  table_sample <- table(a[complete.cases(a), "sample"])
  
  keep_samps <- names(table_sample)[table_sample > 0]
  
  a <- a[a$sample %in% keep_samps, , drop = FALSE]
  
  md <- md[md$shortname %in% keep_samps, , drop = FALSE]
  
  ## drop unused levels
  md$response <- factor(md$response)
  md$day <- factor(md$day)
  md$patient_id <- factor(md$patient_id)
  
}

if(analysis_type == "all"){
  
  min_cells <- 50
  
  table_samp <- table(samp)
  
  keep_samps <- names(table_samp)[which(table_samp > min_cells)]
  
  
  colnames(e) <- fcs_panel$Antigen
  
  a <- aggregate(e, by = list(samp), FUN = median)
  
  colnames(a)[1] <- c("sample")
  
  a$cluster <- -1
  a$label <- "all"
  
  a <- a[, c("cluster", "label", "sample", fcs_panel$Antigen[c(scols, xcols)])]
  a$label <- factor(a$label)
  
  a[!a$sample %in% keep_samps, fcs_panel$Antigen[c(scols, xcols)]] <- NA
  
  ### Save the median expression per cluster and sample
  write.table(a, file.path(outdir, paste0(prefix, "expr_all.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  ### Keep only those samples that have enough cells 
  
  a <- a[a$sample %in% keep_samps, , drop = FALSE]
  
  md <- md[md$shortname %in% keep_samps, , drop = FALSE]
  
  ## drop unused levels
  md$response <- factor(md$response)
  md$day <- factor(md$day)
  md$patient_id <- factor(md$patient_id)
  
  
}



# -----------------------------------------------------------------------------
# Plot expression per cluster
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
### normalize the expression
# -----------------------------------------------------------------------------

### Prepare the matrix with data (rows - markers X clusters; columns - samples)

expr <- a
exprm <- melt(expr, id.vars = c("cluster", "label", "sample"), value.name = "expr", variable.name = "marker")
exprc <- dcast(exprm, cluster + label + marker ~ sample, value.var = "expr")


expr_norm <- exprc[, c("cluster", "label", "marker", md[md$response != "HD", "shortname"])]
th <- 2.5

days <- levels(md$day)

### Normalized to mean = 0 and sd = 1 per day
for(i in days){
  # i = "base"
  expr_norm[, md[md$response != "HD" & md$day == i, "shortname"]] <- t(apply(expr_norm[, md[md$response != "HD" & md$day == i, "shortname"], drop = FALSE], 1, function(x){ 
    
    if(sum(!is.na(x)) == 0)
      return(x)
    
    if(sum(!is.na(x)) < 2)
      return(x-mean(x, na.rm = TRUE))
    
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

levels(md$day)
levels(md$response)


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
  write.table(pvs, file = file.path(outdir, paste0(prefix, "expr_", out_name, "_pvs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  write.table(coeffs, file = file.path(outdir, paste0(prefix, "expr_", out_name, "_coeffs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  
  
  # ----------------------------------------
  # Plot a heatmap of significant cases
  # ----------------------------------------
  
  ### add p-value info
  expr_all <- merge(pvs, expr_norm, by = c("cluster", "label", "marker"), all.x = TRUE, sort = FALSE)
  
  plot_heatmaps_for_sign_expr()
  
  
} # models2fit







sessionInfo()














################################
### 04_expression done!
################################