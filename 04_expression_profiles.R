

Sys.time()

# Load packages
library(flowCore)
library(gdata)
library(ggplot2)
library(reshape2)
library(limma) # for strsplit2



##############################################################################
# Test arguments
##############################################################################


prefix='23CD4_01CD4_pca1_merging5_lmer_interglht_'
outdir='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/01_CD4_noHD10/080_expression/3responses_both/expression_profiles'
path_data='../carsten_cytof/PD1_project/CK_2016-06-23_01_CD4_mergingNEW2/010_data/23CD4_01CD4_expr_raw.rds'
path_clustering_observables='../carsten_cytof/PD1_project/CK_2016-06-23_01_CD4_mergingNEW2/030_heatmaps/23CD4_01CD4_pca1_clustering_observables.xls'
path_clustering='../carsten_cytof/PD1_project/CK_2016-06-23_01_CD4_mergingNEW2/030_heatmaps/23CD4_01CD4_pca1_merging5_clustering.xls'
path_clustering_labels='../carsten_cytof/PD1_project/CK_2016-06-23_01_CD4_mergingNEW2/030_heatmaps/23CD4_01CD4_pca1_merging5_clustering_labels.xls'
path_pvs='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/01_CD4_noHD10/080_expression/3responses_both/01CD4_23CD4merging5_29CD4merging5_all_expr_pvs_lmer_interglht_top10.xls'
path_metadata='../carsten_cytof/PD1_project/CK_metadata/metadata_23_01.xlsx'
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
suffix
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

head(md)

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
# Load pvalues
# ------------------------------------------------------------

pvs <- read.table(path_pvs, header = TRUE, sep = "\t", as.is = TRUE)

comparisons <- colnames(pvs)[grep("adjp_", colnames(pvs))]
comparisons

# ------------------------------------------------------------
# Load expression data
# ------------------------------------------------------------

expr <- readRDS(path_data)

fcs_colnames <- colnames(expr)[!grepl("cell_id|sample_id", colnames(expr))]
e <- expr[, fcs_colnames]

head(e)

# ------------------------------------------------------------
# Load clustering data
# ------------------------------------------------------------

# clustering observables
clustering_observables <- read.table(path_clustering_observables, header = TRUE, sep = "\t", as.is = TRUE)
rownames(clustering_observables) <- clustering_observables$mass
clustering_observables

clust_observ <- clustering_observables[clustering_observables$clustering_observable, "mass"]
clust_observ


# clustering labels
labels <- read.table(path_clustering_labels, header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))
rownames(labels) <- labels$cluster
labels


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
# Plot expression distribution stratified by cluster and sample for the significant markers 
# ------------------------------------------------------------


df <- data.frame(samp = samp, clust = clust, e, check.names = FALSE)
dfm <- melt(df, id.var = c("samp", "clust"), variable.name = "mass", value.name = "expression")

mm <- match(dfm$samp, md$shortname)
dfm$day <- md$day[mm]
dfm$day <- factor(dfm$day, levels = levels(md$day))
dfm$response <- md$response[mm]
dfm$response <- factor(dfm$response, levels = levels(md$response))


mm <- match(dfm$mass, clustering_observables$mass)
dfm$marker <- clustering_observables$marker[mm]

mm <- match(dfm$clust, labels$cluster)
dfm$label <- labels$label[mm]
dfm$label <- factor(dfm$label, levels = levels(labels$label))

head(dfm)


comparisons

for(i in 1:length(comparisons)){
  # i = 1
  
  comparison <- comparisons[i]
  print(comparison)
  comparison_suffix <- paste0(gsub("adjp_", "", comparison), suffix)
  comparison_suffix
  
  pvs_sign <- pvs[pvs[, comparison] < FDR_cutoff & !is.na(pvs[, comparison]), , drop = FALSE]
  print(pvs_sign)
  
  
  if(nrow(pvs_sign) == 0){
    
    
    
  }else{
    
    stopifnot(all(pvs_sign$marker %in% dfm$marker))
    
    pvs_sign <- pvs_sign[order(pvs_sign[, comparison], decreasing = FALSE), , drop = FALSE]
    
    ggdf <- dfm[dfm$marker %in% pvs_sign$marker, , drop = FALSE]
    ggdf$marker <- factor(ggdf$marker, levels = pvs_sign$marker)
    
    ww <- 2 * nlevels(ggdf$marker)
    hh <- 2 * nlevels(ggdf$label)
    
    ggp <- ggplot(ggdf, aes(x = expression, color = samp, linetype = day)) +
      geom_density(adjust = 1) +
      theme_bw() +
      theme(axis.text = element_text(size = 10), 
        axis.title  = element_text(size = 10),
        strip.text = element_text(size = 8, hjust = 0),
        strip.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(size = 0.3, linetype = "solid", color = "black")) +
      guides(color = guide_legend(nrow = 4)) +
      scale_color_manual(values = color_samples) +
      facet_grid(label ~ marker, scales = "free")
    
    
    pdf(file.path(outdir, paste0(prefix, "density_", comparison_suffix,"_all.pdf")), width = ww, height = hh)
    print(ggp)
    dev.off()
    
    ### Plots per marker
    for(j in 1:nrow(pvs_sign)){
      # j = 1
      
      marker <- pvs_sign$marker[j]
      print(marker)
      ggdf <- dfm[dfm$marker %in% marker & dfm$day %in% "base" & !dfm$response %in% "HD", , drop = FALSE]
      
      ww <- 10
      hh <- 7
      
      ggp <- ggplot(ggdf, aes(x = expression, color = label, linetype = day)) +
        geom_density(adjust = 1, position = "stack") +
        theme_bw() +
        theme(axis.text = element_text(size = 10), 
          axis.title  = element_text(size = 10),
          strip.text = element_text(size = 8, hjust = 0),
          strip.background = element_blank(),
          legend.key = element_blank(),
          legend.title = element_blank(),
          legend.position = "right",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_rect(size = 0.3, linetype = "solid", color = "black")) +
        guides(color = guide_legend(nrow = 4)) +
        # scale_color_manual(values = color_samples) +
        facet_wrap(~ samp, scales = "free")
      
      
      pdf(file.path(outdir, paste0(prefix, "density_", comparison_suffix,"_", marker ,"_stack.pdf")), width = ww, height = hh)
      print(ggp)
      dev.off()
      
      
    }
    
    
    
  }
  
  
}






































sessionInfo()






