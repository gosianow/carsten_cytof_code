##############################################################################
## <<08_frequencies_merged_bimatrix>>

# BioC 3.3
# Created 1 Dec 2016

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
library(UpSetR)

##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-merged_23_29/02_CD4'
freq_prefix='23_29_CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_cd69_positive_'
freq_outdir='10_cd69_merged_overall_freqs/03_frequencies_auto_2responses'
path_data=c('/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_02_CD4_merging2/070_cd69_bimatrix/01_clustering/23CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_cd69_positive_bimatrix.txt','/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-29_02_CD4_merging/070_cd69_bimatrix/01_clustering/29CD4_02CD4_pca1_merging3_Tmem_cytCM_raw2_cd69_positive_bimatrix.txt')
path_clustering_observables=c('/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_02_CD4_merging2/070_cd69_bimatrix/01_clustering/23CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_cd69_positive_clustering_observables.xls','/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-29_02_CD4_merging/070_cd69_bimatrix/01_clustering/29CD4_02CD4_pca1_merging3_Tmem_cytCM_raw2_cd69_positive_clustering_observables.xls')
path_metadata=c('/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_02.xlsx','/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_29_02.xlsx')
data_name=c('data23','data29')
path_fun_models='/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_models.R'
path_fun_formulas='/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_formulas_2datasets_2responses.R'
path_fun_plot_frequencies <- "/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_plot_frequencies.R"
path_fun_plot_heatmaps <- "/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_plot_heatmaps_for_sign_freqs.R"


### Optional arguments
pdf_hight=4
plot_only=FALSE
FDR_cutoff=0.1
suffix='_top01'


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

### Optional arguments
if(!any(grepl("pdf_hight=", args))){
  pdf_hight=4
}

if(!any(grepl("plot_only=", args))){
  plot_only=FALSE
}

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
# Load data - bimatrix
# ------------------------------------------------------------


expr <- lapply(1:length(data_name), function(i){
  
  if(grepl(".txt", path_data[i])){
    expr <- read.table(path_data[i], header = TRUE, sep = "\t", as.is = TRUE)
  }
  if(grepl(".rds", path_data[i])){
    expr <- readRDS(path_data[i])
  }
  
  expr$data_id <- data_name[i]
  return(expr)
})

expr <- rbind.fill(expr)

### cell_id is a paste of cell_id and sample_id because cell ids are the same for data 23 and 29 
cell_id <- paste0(expr[, "cell_id"])
samp <- expr[, "sample_id"]
fcs_colnames <- colnames(expr)[!grepl("cell_id|sample_id|data_id", colnames(expr))]
e <- expr[, fcs_colnames]

data_id <- paste0(expr[, "data_id"])

# ------------------------------------------------------------
# Load clustering_observables
# ------------------------------------------------------------
# Use an intersection of observables used in each dataset

clustering_observables <- lapply(1:length(data_name), function(i){
  
  clustering_observables <- read.table(path_clustering_observables[i], header = TRUE, sep = "\t", as.is = TRUE)
  return(clustering_observables)
  
})

clustering_observables <- Reduce(function(...) merge(..., by = c("mass", "marker"), all=TRUE, sort = FALSE), clustering_observables)

clustering_observables$clustering_observable <- apply(clustering_observables[, grep("clustering_observable", colnames(clustering_observables))], 1, all)

clust_observ <- clustering_observables[clustering_observables$clustering_observable, "mass"]


# selected columns for clustering 

scols <- which(fcs_colnames %in% clust_observ)

ef <- as.matrix(e[, scols])


# ------------------------------------------------------------
# Upsetr plots
# ------------------------------------------------------------

bidf <- data.frame(ef[, clustering_observables[clustering_observables$clustering_observable, "mass"]], row.names = 1:nrow(ef), check.names = FALSE)
colnames(bidf) <- clustering_observables[clustering_observables$clustering_observable, "marker"]

pdf(file.path(outdir, paste0(prefix, "upsetr.pdf")), w = 16, h = 6)
upset(bidf, sets = colnames(bidf), nintersects = 50, order.by = "freq")
dev.off()


# --------------------------------------------------------------------------
# Calculate the positive marker frequencies
# --------------------------------------------------------------------------

clust_observ <- colnames(bidf)

freq_list <- lapply(clust_observ, function(i){
  # i <- "PD-1"
  
  freq_tmp <- table(bidf[, i], samp)
  rownames(freq_tmp) <- paste0(i, c("-", "+"))
  
  freq_tmp
  
})

prop_list <- lapply(freq_list, function(x){
  
  prop_tmp <- t(t(x[2, , drop = FALSE]) / colSums(x, na.rm = TRUE)) * 100
  
})

prop <- do.call(rbind, prop_list)
prop_out <- data.frame(cluster = rownames(prop), label = rownames(prop), prop)


freq <- do.call(rbind, freq_list)
freq_out <- data.frame(cluster = rownames(freq), label = rownames(freq), freq)


write.table(prop_out, file=file.path(outdir, paste0(prefix, "frequencies.xls")), row.names=FALSE, quote=FALSE, sep="\t")
write.table(freq_out, file=file.path(outdir, paste0(prefix, "counts.xls")), row.names=FALSE, quote=FALSE, sep="\t")



freq_out_list <- lapply(freq_list, function(x){
  
  data.frame(cluster = rownames(x), label = rownames(x), as.data.frame.matrix(x))
  
})


# ------------------------------------------------------------
# Plot frequencies
# ------------------------------------------------------------

source(path_fun_plot_frequencies)


ggdf <- melt(prop_out, id.vars = c("cluster", "label"), value.name = "prop", variable.name = "samp")

## use labels as clusters
ggdf$cluster <- factor(ggdf$label, levels = rownames(prop_out))
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


plot_frequencies(ggdf = ggdf, color_groups = color_groups, color_groupsb = color_groupsb, outdir = outdir, prefix = prefix, pdf_hight = pdf_hight)



# -----------------------------------------------------------------------------
# Prepare the matrix with data for heatmaps (rows - clusters; columns - samples)
# Plot a heatmap with clustered columns
# -----------------------------------------------------------------------------
# Transform proportions with arcsin-sqrt so the dispersion is the same for low and high props.


### normalize the expression

prop_list <- lapply(freq_list, function(x){
  
  prop_tmp <- asin(sqrt(t(t(x[2, , drop = FALSE]) / colSums(x, na.rm = TRUE))))
  
})

prop <- do.call(rbind, prop_list)
ass_freq_out <- data.frame(cluster = rownames(prop), label = rownames(prop), prop)


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
  
  models2fit <- c("glmer_binomial_interglht")
  
  
  for(k in models2fit){
    # k = "glmer_binomial_interglht"
    print(k)
    
    
    out_list <- lapply(1:length(freq_out_list), function(i){
      # i = 1
      
      freq_out <- freq_out_list[[i]]
      
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
      
      return(list(pvs = pvs, coeffs = coeffs))
      
    })
    
    ## Extracs p-values
    pvs <- lapply(out_list, function(x){
      x[["pvs"]][2, , drop = FALSE]
    })
    
    pvs <- do.call(rbind, pvs)
    
    ## Readjust the p-values
    pval_colnames <- colnames(pvs)[grep("pval_", colnames(pvs))]
    
    for(i in 1:length(pval_colnames)){
      # i = 1
      pvs[, gsub("pval_", "adjp_", pval_colnames[i])] <- p.adjust(pvs[, pval_colnames[i]], method = "BH")
    }
    
    ## Extracs coefficients
    coeffs <- lapply(out_list, function(x){
      x[["coeffs"]][2, , drop = FALSE]
    })
    
    coeffs <- do.call(rbind, coeffs)
    
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














################################################################
### 08_frequencies_merged_bimatrix done!
################################################################