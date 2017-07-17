

Sys.time()

# Load packages
library(gdata)
library(ggplot2)
library(reshape2)
library(limma) # for strsplit2
library(RColorBrewer)
library(pheatmap)
library(gtools) # for logit
library(plyr) # for rbind.fill

##############################################################################
# Test arguments
##############################################################################


prefix='01CD4_23CD4merging5_29CD4merging5_'
outdir='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/01_CD4_noHD10/050_frequencies/HDvsRandNR'
path_metadata=c('../carsten_cytof/PD1_project/CK_metadata/metadata_23_01.xlsx','../carsten_cytof/PD1_project/CK_metadata/metadata_29_01_noHD10.xlsx')
path_counts=c('../carsten_cytof/PD1_project/CK_2016-06-23_01_CD4_mergingNEW2/050_frequencies/23CD4_01CD4_pca1_merging5_counts.xls','../carsten_cytof/PD1_project/CK_2016-06-29_01_CD4_merging2/050_frequencies/29CD4_01CD4_pca1_merging5_counts.xls')
path_fun_models='./00_models.R'
FDR_cutoff='10'


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

if(!file.exists(outdir)) 
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
# Load cluster frequencies
# ------------------------------------------------------------


freq <- lapply(1:length(path_counts), function(i){
  # i = 1
  freq <- read.table(path_counts[i], header = TRUE, sep = "\t", as.is = TRUE)
})

freq_out <- Reduce(function(...) merge(..., by = c("cluster", "label"), all=TRUE, sort = FALSE), freq)

## Drop the 'drop' cluster
freq_out <- freq_out[freq_out$label != "drop", , drop = FALSE]

if(!all(complete.cases(freq_out))){
  stop("There are some clusters that are not common in the merged data sets or have different cluster number!")
}


## keep only those samples that are also in the metadata file
freq_out <- freq_out[, c("cluster", "label", md$shortname)]

freq_out <- freq_out[order(freq_out$cluster), , drop = FALSE]

labels <- data.frame(cluster = freq_out$cluster, label = freq_out$label)
labels$label <- factor(labels$label, levels = unique(labels$label))
labels


# ---------------------------------------
# Keep only those samples that have enough cells
# Samples with 0 cells are eliminated from the analysis
# ---------------------------------------

min_cells <- 0

table_samp <- colSums(freq_out[md$shortname], na.rm = TRUE)
names(table_samp) <- md$shortname

keep_samps <- names(table_samp)[which(table_samp > min_cells)]

freq_out <- freq_out[, colnames(freq_out) %in% c("cluster", "label", keep_samps), drop = FALSE]

md <- md[md$shortname %in% keep_samps, , drop = FALSE]

## drop unused levels
md$response <- factor(md$response)
md$day <- factor(md$day)
md$patient_id <- factor(md$patient_id)



# ------------------------------------------------------------
# Test for frequency differences between groups
# ------------------------------------------------------------
## The model functions do not anlyse a cluster with NAs; 
## For merged data it means such cluster was not present in all the datasets
## For expression data clusters with no cells are skipped


### Load functions fitting models
source(path_fun_models)


### Load formulas that are fit in the models - this function may change the md object!!!

md$response2 <- ifelse(md$response == "HD", "HD", "cancer")
md$response2 <- factor(md$response2, levels = c("HD", "cancer")) 


model.matrix( ~ response2 + data_day + response2:data_day, data = md)


if(identical(levels(md$data), c("data23", "data29")) && identical(levels(md$day), c("base", "tx")) && identical(levels(md$response2), c("HD", "cancer")) && identical(levels(md$data_day), c("data23.base", "data23.tx", "data29.base", "data29.tx"))){
  
  ## create formulas
  formula_glmer_binomial <- y/total ~ response2 + data_day + response2:data_day + (1|patient_id)
  
  ## create contrasts
  contrast_names <- c("HDvsCancer", "HDvsCancer_base", "HDvsCancer_tx", "HDvsCancer_basevstx")
  k1 <- c(0, 1, 0, 0, 0, 0, 1/2, 0) # HDvsCancer in base
  k2 <- c(0, 1, 0, 0, 0, 1/2, 0, 1/2) # HDvsCancer in tx
  k0 <- (k1 + k2) / 2 # HDvsCancer
  k3 <- c(0, 0, 0, 0, 0, 1/2, 0, 1/2) # whether HDvsCancer is different in base and tx
  K <- matrix(c(k0, k1, k2, k3), nrow = 4, byrow = TRUE)
  rownames(K) <- contrast_names
  
  ### p-value for sorting the output
  pval_name1 <- "pval_HDvsCancer"
  ### p-value for plotting the pheatmap2
  adjpval_name2 <- "adjp_HDvsCancer"
  pval_name2 <- "pval_HDvsCancer"
  ### p-value for plotting the pheatmap3
  adjpval_name_list <- c("adjp_HDvsCancer", "adjp_HDvsCancer_base", "adjp_HDvsCancer_tx", "adjp_HDvsCancer_basevstx")
  pval_name_list <- c("pval_HDvsCancer", "pval_HDvsCancer_base", "pval_HDvsCancer_tx", "pval_HDvsCancer_basevstx")
  
}else{
  
  stop("Metadata does not fit to any the models that are specified !!!")
  
}







k <- "glmer_binomial_interglht"


# Fit a GLMM binomial with interactions + test contrasts with multcomp pckg
fit_out <- fit_glmer_interglht(data = freq_out, md, family = "binomial", formula = formula_glmer_binomial, K = K)


# ----------------------------------------
# Extract p-values and coeffs
# ----------------------------------------

pvs <- data.frame(freq_out[, c("cluster", "label")], fit_out[["pvals"]])
coeffs <- data.frame(freq_out[, c("cluster", "label")], fit_out[["coeffs"]])

oo <- order(pvs[, pval_name1], decreasing = FALSE)
pvs <- pvs[oo, , drop = FALSE]
coeffs <- coeffs[oo, , drop = FALSE]

## save the results
write.table(pvs, file=file.path(outdir, paste0(prefix, "frequencies_pvs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
write.table(coeffs, file=file.path(outdir, paste0(prefix, "frequencies_coeffs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")














sessionInfo()




