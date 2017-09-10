

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

prefix='23CD8allall_29CD8allall_02CD8v2_cl49_top10_glmer_binomial_interglht_'
outdir='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/02v2_CD8/090_cytokine_bimatrix_frequencies_clustering/3responses_both'
path_metadata=c('../carsten_cytof/PD1_project/CK_metadata/metadata_23_02.xlsx','../carsten_cytof/PD1_project/CK_metadata/metadata_29_02.xlsx')
path_frequencies='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/02v2_CD8/090_cytokine_bimatrix_frequencies_clustering/3responses_both/23CD8allall_29CD8allall_02CD8v2_cl49_frequencies.xls'
path_pvs='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/02v2_CD8/090_cytokine_bimatrix_frequencies_clustering/3responses_both/23CD8allall_29CD8allall_02CD8v2_cl49_frequencies_pvs_glmer_binomial_interglht_top10.xls'
path_fun_plot_frequencies='./00_plot_frequencies.R'
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
md$response <- factor(md$response, levels = c("HD", "NR", "R")) ### For plotting HD is a reference!!!
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

prop <- lapply(1:length(path_frequencies), function(i){
  # i = 1
  prop <- read.table(path_frequencies[i], header = TRUE, sep = "\t", as.is = TRUE)
  print(prop[, c("cluster", "label")])
  return(prop)
})


prop_out <- Reduce(function(...) merge(..., by = c("cluster", "label"), all=TRUE, sort = FALSE), prop)
prop_out[, c("cluster", "label")]

## drop the 'drop' cluster
prop_out <- prop_out[prop_out$label != "drop", , drop = FALSE]

if(!all(complete.cases(prop_out))){
  stop("There are some clusters that are not common in the merged data sets or have different cluster number!")
}

## keep only those samples that are also in the metadata file
prop_out <- prop_out[, c("cluster", "label", md$shortname)]

prop_out <- prop_out[order(prop_out$cluster), , drop = FALSE]

labels <- data.frame(cluster = prop_out$cluster, label = prop_out$label)
labels$label <- factor(labels$label, levels = unique(labels$label))
labels

# ------------------------------------------------------------
# Load pvalues
# ------------------------------------------------------------

pvs <- read.table(path_pvs, header = TRUE, sep = "\t", as.is = TRUE)


comparisons <- colnames(pvs)[grep("adjp_", colnames(pvs))]



# ------------------------------------------------------------
# Colors for clusters
# ------------------------------------------------------------

# ggplot palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=60 , c=100)[1:n]
}

# color blind palette

colors_muted <- c("#DC050C", "#E8601C", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6", "#F1932D", "#F6C141", "#F7EE55", "#4EB265", "#90C987", "#CAEDAB")
color_ramp <- c(colors_muted, gg_color_hue(max(1, nlevels(labels$label) - length(colors_muted))))

colors_clusters <- color_ramp[1:nlevels(labels$label)]
names(colors_clusters) <- levels(labels$label)


# ------------------------------------------------------------
# Plot frequencies
# ------------------------------------------------------------

source(path_fun_plot_frequencies)


for(i in 1:length(comparisons)){
  # i = 1
  
  comparison <- comparisons[i]
  print(comparison)
  comparison_prefix <- paste0(gsub("adjp_", "", comparison), "_")
  
  pvs_sign <- pvs[pvs[, comparison] < FDR_cutoff & !is.na(pvs[, comparison]), , drop = FALSE]
  print(pvs_sign)
  
  if(nrow(pvs_sign) == 0){
    
    ggdf <- NULL
    
  }else{
    
    ### Order by the significance
    pvs_sign <- pvs_sign[order(pvs_sign[, gsub("adjp_", "pval_", comparison)]), , drop = FALSE]
    
    labels_sign <- as.character(pvs_sign$label)
    
    prop_sign <- prop_out[prop_out$label %in% labels_sign, , drop = FALSE]
    
    ggdf <- melt(prop_sign, id.vars = c("cluster", "label"), value.name = "prop", variable.name = "samp")
    
    ## use labels as clusters
    ggdf$cluster <- factor(ggdf$label, levels = labels_sign)
    
    # add more info about samples
    mm <- match(ggdf$samp, md$shortname)
    ggdf$group <- factor(md$condition[mm])
    ggdf$day <- factor(md$day[mm])
    ggdf$data <- factor(md$data[mm])
    ggdf$response <- factor(md$response[mm])
    
    ## replace _ with \n
    levels(ggdf$group) <- gsub("_", "\n", levels(ggdf$group))
    
  }
  
  prefix_sign <- paste0(prefix, comparison_prefix)
  
  plot_frequencies(ggdf = ggdf, color_groups = color_groups, color_groupsb = color_groupsb, colors_clusters = colors_clusters, color_response = color_response, outdir = outdir, prefix = prefix_sign)
  
  
}

















sessionInfo()




