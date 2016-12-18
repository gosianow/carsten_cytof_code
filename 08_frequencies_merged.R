##############################################################################
## <<08_frequencies_merged.R>>

# BioC 3.3
# Created 22 Sep 2016
# Updated 20 Oct 2016

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

##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-merged_23_29/01'
freq_prefix='23m6_29m4_'
freq_outdir='08_frequencies_merged_2responses'
path_metadata=c('/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.xlsx','/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_29_01.xlsx')
path_counts=c('/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01/050_frequencies/23_01_pca1_merging6_counts.xls','/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-29_01/050_frequencies/29_01_pca1_merging4_counts.xls')
data_name=c('data23','data29')
path_fun_models='/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_models.R'
path_fun_formulas='/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_formulas_2datasets_2responses.R'
path_fun_plot_frequencies <- "/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_plot_frequencies.R"
path_fun_plot_heatmaps <- "/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_plot_heatmaps_for_sign_freqs.R"

### Optional arguments
pdf_hight=8
plot_only=FALSE
FDR_cutoff=0.05
suffix="_top005"
# For selecting and ordering clusters to plot
path_cluster_selection='23CD4m5_29CD4m5_frequencies_cluster_selection.txt' 

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
outdir <- freq_outdir

if(!file.exists(outdir)) 
  dir.create(outdir, recursive = TRUE)


if(!any(grepl("pdf_hight=", args))){
  pdf_hight=4
}

if(!any(grepl("FDR_cutoff=", args))){
  FDR_cutoff=0.05
}

if(!any(grepl("suffix=", args))){
  suffix="_top005"
}

if(!any(grepl("plot_only=", args))){
  plot_only=FALSE
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

if(all(!complete.cases(freq_out)))
  stop("There is no common cluster in the merged data sets!")

prop_out <- freq_out
prop_out[md$shortname] <- t(t(freq_out[md$shortname]) / colSums(freq_out[md$shortname], na.rm = TRUE)) * 100


labels <- data.frame(cluster = freq_out$cluster, label = freq_out$label)
labels$label <- factor(labels$label, levels = unique(labels$label))
labels


write.table(prop_out, file=file.path(outdir, paste0(prefix, "frequencies.xls")), row.names=FALSE, quote=FALSE, sep="\t")
write.table(freq_out, file=file.path(outdir, paste0(prefix, "counts.xls")), row.names=FALSE, quote=FALSE, sep="\t")

# ---------------------------------------
# Keep only those samples that have enough cells 
# ---------------------------------------

min_cells <- 0

table_samp <- colSums(freq_out[md$shortname], na.rm = TRUE)
names(table_samp) <- md$shortname

keep_samps <- names(table_samp)[which(table_samp > min_cells)]


prop_out <- prop_out[, colnames(prop_out) %in% c("cluster", "label", keep_samps), drop = FALSE]
freq_out <- freq_out[, colnames(freq_out) %in% c("cluster", "label", keep_samps), drop = FALSE]

table(md$shortname %in% keep_samps)

md <- md[md$shortname %in% keep_samps, , drop = FALSE]


## drop unused levels
md$response <- factor(md$response)
md$day <- factor(md$day)
md$patient_id <- factor(md$patient_id)
md$data <- factor(md$data)
md$data_day <- factor(md$data_day)



# ------------------------------------------------------------
# Plot frequencies
# ------------------------------------------------------------

source(path_fun_plot_frequencies)


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

## add more info about samples
ggdf$day <- strsplit2(ggdf$group, "\n")[, 1]
ggdf$day <- factor(ggdf$day)



### Cluster selection
### Use only the clusters that are specified in the cluster_selection.txt file

if(any(grepl("path_cluster_selection=", args))){
  
  if(file.exists(path_cluster_selection)){
    
    cluster_selection <- read.table(path_cluster_selection, header = FALSE, sep = "\t", as.is = TRUE)[, 1]
    
    if(!all(cluster_selection %in% levels(ggdf$cluster)))
      stop("Cluster selection is wrong!!!")
    
    ggdf <- ggdf[ggdf$cluster %in% cluster_selection, , drop = FALSE]
    ggdf$cluster <- factor(ggdf$cluster, levels = cluster_selection)
    
  }else{
    
    ggdf$prop <- NA
    
  }
  
}


plot_frequencies(ggdf = ggdf, color_groups = color_groups, color_groupsb = color_groupsb, outdir = outdir, prefix = prefix, pdf_hight = pdf_hight)


# -----------------------------------------------------------------------------
# Prepare the matrix with data for heatmaps (rows - clusters; columns - samples)
# Plot a heatmap with clustered columns
# -----------------------------------------------------------------------------
# Transform proportions with arcsin-sqrt so the dispersion is the same for low and high props.


### normalize the expression
ass_freq_out <- freq_out
ass_freq_out[md$shortname] <- asin(sqrt((t(t(freq_out[md$shortname]) / colSums(freq_out[md$shortname], na.rm = TRUE)))))

ass_freq_out <- ass_freq_out[complete.cases(ass_freq_out), , drop = FALSE]

expr_norm <- ass_freq_out[, c("cluster", "label", md[md$response != "HD", "shortname"])]
th <- 2.5

data_days <- levels(md$data_day)

### Normalized to mean = 0 and sd = 1 per data and day
for(i in data_days){
  # i = "data23.base"
  expr_norm[, md[md$response != "HD" & md$data_day == i, "shortname"]] <- t(apply(expr_norm[, md[md$response != "HD" & md$data_day == i, "shortname"], drop = FALSE], 1, function(x){
    
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


### Plot a heatmap with clustered columns and all the rows

expr_heat <- expr_norm
rownames(expr_heat) <- expr_heat$label

expr <- expr_heat[, md[md$response != "HD", "shortname"]]

labels_row <- paste0(expr_heat$label) 
labels_col <- colnames(expr)

annotation_col <- data.frame(response = factor(md[md$response != "HD", "response"]))
rownames(annotation_col) <- md[md$response != "HD", "shortname"]

annotation_colors <- list(response = color_response[levels(annotation_col$response)])

cluster_cols <- hclust(dist(t(expr)), method = "ward.D2")
cluster_rows <- hclust(dist(expr), method = "ward.D2")

# Using pheatmap
pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = cluster_cols, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_col = annotation_col, annotation_colors = annotation_colors, annotation_legend = TRUE, filename = file.path(outdir, paste0(prefix, "frequencies_pheatmap_colclust", ".pdf")))





# ------------------------------------------------------------
# Test for frequency differences between groups
# ------------------------------------------------------------
## The model functions do not anlyse a cluster with NAs; 
## For merged data it means such cluster was not present in all the datasets
## For expression data clusters with no cells are skipped

if(!plot_only){
  
  ### Load functions fitting models
  source(path_fun_models)
  ### Load formulas that are fit in the models - this function may change the md object!!!
  source(path_fun_formulas)
  
  source(path_fun_plot_heatmaps)
  
  
  levels(md$data)
  levels(md$day)
  levels(md$response)
  
  # models2fit <- c("glm_binomial_interglht", "glm_quasibinomial_interglht", "glmer_binomial_interglht", "lmer_logit_interglht", "lmer_arcsinesqrt_interglht", "lm_logit_interglht", "lm_arcsinesqrt_interglht")
  models2fit <- c("glmer_binomial_interglht")
  
  for(k in models2fit){
    # k = "glmer_logit_interglht"
    print(k)
    
    switch(k,
      glm_binomial_interglht = {
        # Fit a GLM binomial with interactions + test contrasts with multcomp pckg
        fit_out <- fit_glm_interglht(data = freq_out, md, family = "binomial", formula = formula_glm_binomial, K = K)
        
      },
      glm_quasibinomial_interglht = {
        # Fit a GLM quasibinomial with interactions + test contrasts with multcomp pckg
        fit_out <- fit_glm_interglht(data = freq_out, md, family = "quasibinomial", formula = formula_glm_binomial, K = K)
        
      },
      glmer_binomial_interglht = {
        # Fit a GLMM binomial with interactions + test contrasts with multcomp pckg
        fit_out <- fit_glmer_interglht(data = freq_out, md, family = "binomial", formula = formula_glmer_binomial, K = K)
        
      },
      lmer_logit_interglht = {
        
        logit_freq_out <- freq_out
        logit_freq_out[md$shortname] <- logit(t(t(freq_out[md$shortname]) / colSums(freq_out[md$shortname], na.rm = TRUE)))
        ## Be carefull about Inf and -Inf for prop = 0, 1
        fit_out <- fit_lmer_interglht(data = logit_freq_out, md, formula = formula_lmer, K = K)
        
      },
      lm_logit_interglht = {
        
        logit_freq_out <- freq_out
        logit_freq_out[md$shortname] <- logit(t(t(freq_out[md$shortname]) / colSums(freq_out[md$shortname], na.rm = TRUE)))
        ## Be carefull about Inf and -Inf for prop = 0, 1
        fit_out <- fit_lm_interglht(data = logit_freq_out, md, formula = formula_lm, K = K)
        
      },
      lmer_arcsinesqrt_interglht = {
        
        ass_freq_out <- freq_out
        ass_freq_out[md$shortname] <- asin(sqrt((t(t(freq_out[md$shortname]) / colSums(freq_out[md$shortname], na.rm = TRUE)))))
        
        fit_out <- fit_lmer_interglht(data = ass_freq_out, md, formula = formula_lmer, K = K)
        
      },
      lm_arcsinesqrt_interglht = {
        
        ass_freq_out <- freq_out
        ass_freq_out[md$shortname] <- asin(sqrt((t(t(freq_out[md$shortname]) / colSums(freq_out[md$shortname], na.rm = TRUE)))))
        
        fit_out <- fit_lm_interglht(data = ass_freq_out, md, formula = formula_lm, K = K)
        
      },
      
      glmmadmb_fixed_beta_interglht = {
        fit_out <- fit_glm_interglht(data = freq_out, md, family = "beta", formula = formula_glm_beta, K = K)
      },
      
      glmmadmb_fixed_betabinomial_interglht = {
        fit_out <- fit_glm_interglht(data = freq_out, md, family = "betabinomial", formula = formula_glm_binomial, K = K)
      },
      
      test_wilcoxon = {
        
        n01_freq_out <- freq_out
        n01_freq_out[md$shortname] <- t(t(freq_out[md$shortname]) / colSums(freq_out[md$shortname]))
        
        fit_out <- test_wilcoxon(data = n01_freq_out, md)
      }
      
    )
    
    # ----------------------------------------
    # Extract p-values and coeffs
    # ----------------------------------------
    
    pvs <- data.frame(freq_out[, c("cluster", "label"), drop = FALSE], fit_out[["pvals"]])
    coeffs <- data.frame(freq_out[, c("cluster", "label"), drop = FALSE], fit_out[["coeffs"]])
    
    oo <- order(pvs[, pval_name1], decreasing = FALSE)
    pvs <- pvs[oo, , drop = FALSE]
    coeffs <- coeffs[oo, , drop = FALSE]
    
    ## save the results
    write.table(pvs, file=file.path(outdir, paste0(prefix, "frequencies_pvs_", k, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
    write.table(coeffs, file=file.path(outdir, paste0(prefix, "frequencies_coeffs_", k, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
    
    
    # ----------------------------------------
    # Plot a heatmap of significant cases 
    # ----------------------------------------
    
    ### add p-value info
    expr_all <- merge(expr_norm, pvs, by = c("cluster", "label"), all.x = TRUE, sort = FALSE)
    
    prefix2 <- paste0(k, "_")
    
    plot_heatmaps_for_sign_freqs(expr_all = expr_all, md = md, FDR_cutoff = FDR_cutoff, pval_name2 = pval_name2, adjpval_name2 = adjpval_name2, pval_name_list = pval_name_list, adjpval_name_list = adjpval_name_list, breaks = breaks, legend_breaks = legend_breaks, outdir = outdir, prefix = prefix, prefix2 = prefix2, suffix = suffix)
    
    
  }
  
  
}





sessionInfo()













################################
### 08_frequencies_merged.R done!
################################