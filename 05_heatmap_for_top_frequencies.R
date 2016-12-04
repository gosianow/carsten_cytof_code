##############################################################################
## <<05_heatmap_for_top_frequencies.R>>

# BioC 3.3
# Created 14 Nov 2016
# Updated 14 nov 2016

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
library(qvalue)

##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-merged_23_29/02_CD4'
freq_prefix='23_29_CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_cl40_'
freq_outdir='10_cytokines_merged/03_frequencies_auto_2responses'
path_pvs='10_cytokines_merged/03_frequencies_auto_2responses/23_29_CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_cl40_frequencies_pvs_glmer_binomial_interglht.xls'
path_clusters='10_cytokines_merged/01_clustering/23_29_CD4_02CD4_pca1_merging2_Tmem_cytCM_raw2_cl40_clusters.xls'
model2fit='glmer_binomial_interglht'
pvs_list=c('NRvsR','NRvsR_base','NRvsR_tx')
adjp_top=NA
adjp_cut=0.05
pheatmap_palette='Greys'
pheatmap_palette_rev=FALSE

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

prefix <- paste0(freq_prefix, "frequencies_", model2fit, "_")
suffix <- ""
outdir <- freq_outdir

if(!file.exists(outdir)) 
  dir.create(outdir, recursive = TRUE)



# ----------------------------------------
# Read in the p-values and clusters data
# ----------------------------------------

pvs <- read.table(path_pvs, header = TRUE, sep = "\t", as.is = TRUE)

clusters <- read.table(path_clusters, header = TRUE, sep = "\t", as.is = TRUE)

## merge only by label because cluster may be different
datam <- merge(pvs, clusters, by = "label", all.x = TRUE)

markers <- colnames(clusters)[!grepl("cluster|label|frequencies|counts", colnames(clusters))]



# ------------------------------------------------------------
# Plot pheatmaps of expression
# ------------------------------------------------------------

if(pheatmap_palette_rev){
  color <- colorRampPalette(rev(brewer.pal(n = 8, name = pheatmap_palette)))(100)
}else{
  color <- colorRampPalette(brewer.pal(n = 8, name = pheatmap_palette))(100)
}


for(i in 1:length(pvs_list)){
  # i = 1
  
  datam <- datam[order(datam[, paste0("pval_", pvs_list[i])], decreasing = FALSE), , drop = FALSE]
  
  if(!is.na(adjp_cut)){
    keep_pvs <- datam[, paste0("adjp_", pvs_list[i])] <= adjp_cut & !is.na(datam[, paste0("adjp_", pvs_list[i])])
    pdf_name <- gsub("\\.", "", as.character(adjp_cut))
    if(sum(keep_pvs, na.rm = TRUE) == 0){
      pdf(file.path(outdir, paste0(prefix, "top", pdf_name, "_pheatmap_", pvs_list[i],".pdf")))
      plot(1, type="n", axes=F, xlab="", ylab="")
      dev.off()
      next
    }
    
  }
  
  if(!is.na(adjp_top)){
    keep_pvs <- intersect(1:min(adjp_top, nrow(datam)), which(!is.na(datam[, paste0("adjp_", pvs_list[i])])))
    pdf_name <- adjp_top
  }
  
  
  expr <- as.matrix(datam[, markers, drop = FALSE])
  rownames(expr) <- datam$label
  
  expr <- expr[keep_pvs, , drop = FALSE]
  
  cluster_rows <- FALSE
  
  labels_row <- paste0(datam$label, " (", sprintf( "%.02e", datam[, paste0("adjp_", pvs_list[i])]), ")")[keep_pvs]
  labels_col <- gsub(pattern = "\\.", replacement = "-" , colnames(expr)) 
  
  
  pheatmap(expr, color = color, cellwidth = 24, cellheight = 24, cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, display_numbers = FALSE, number_color = "black", fontsize_number = 8,  fontsize_row = 10, fontsize_col = 10, fontsize = 6, filename = file.path(outdir, paste0(prefix, "top", pdf_name, "_pheatmap_", pvs_list[i],".pdf")))
  
  ### Save the labels that are in the top list
  
  write.table(datam$label[keep_pvs], file = file.path(outdir, paste0(prefix, "top", pdf_name, "_selection_", pvs_list[i],".txt")), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  
}














sessionInfo()




















################################################################
### 05_heatmap_for_top_frequencies.R done!
################################################################