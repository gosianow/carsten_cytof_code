##############################################################################
## <<06_cytokines.R>>

# BioC 3.3
# Created 5 Aug 2016
# Updated 22 Aug 2016

##############################################################################
Sys.time()
##############################################################################

# Load packages
library(flowCore)
library(gdata)
library(Repitools)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2)
library(coop)
library(RColorBrewer)
library(UpSetR)
library(limma)
library(FlowSOM)
library(ConsensusClusterPlus)
library(pheatmap)


##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_02_CD4_merging2'
cytokines_prefix='23CD4_02CD4_pca1_merging_Tmem_cytCM_raw_'
path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_02.xlsx'
path_cytokines_cutoffs='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_panels/panel2CD4_cytokines_CM.xlsx'
path_clustering='23CD4_02CD4_pca1_merging_clustering.xls'
path_clustering_labels='23CD4_02CD4_pca1_merging_clustering_labels.xls'
clsubset = c('CM','EM','TE')
cutoff_colname='positive_cutoff_raw_base'
data2analyse='raw'
nmetaclusts=40
path_marker_selection='23CD4_02CD4_pca1_merging_Tmem_marker_selection.txt'
path_rtsne_out='23CD4_02CD4_pca1_rtsne_out_norm.rda'
path_rtsne_data='23CD4_02CD4_pca1_rtsne_data_norm.xls'
pdf_width=15
pdf_height=10
tsnep_suffix='_norm'

##############################################################################
# Read in the arguments
##############################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(args)

##############################################################################

setwd(rwd)

prefix <- cytokines_prefix
suffix <- ""

prefix_clust <- paste0("cl", nmetaclusts, "_")

if(!data2analyse %in% c("raw", "norm"))
  stop("data2analyse must be 'raw' or 'norm'!")


# ------------------------------------------------------------
# define directories
# ------------------------------------------------------------

fcsDir <- "010_cleanfcs"; if( !file.exists(fcsDir) ) dir.create(fcsDir)
hmDir <- "030_heatmaps"; if( !file.exists(hmDir) ) dir.create(hmDir)
sneDir <- "040_tsnemaps"; if( !file.exists(sneDir) ) dir.create(sneDir)
cyDir <- "060_cytokines"; if( !file.exists(cyDir) ) dir.create(cyDir)

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

# read metadata
md <- read.xls(path_metadata, stringsAsFactors=FALSE)

# define FCS file names
f <- file.path(fcsDir, md$filename)
names(f) <- md$shortname

# read raw FCS files in
fcs <- lapply(f, read.FCS)

fcs_colnames <- colnames(fcs[[1]])

### Create sample info
samp <- rep(names(fcs), sapply(fcs, nrow))

# ------------------------------------------------------------
# Load more data
# ------------------------------------------------------------

## positive cutoffs for cytokines
cytokines_cutoffs <- read.xls(path_cytokines_cutoffs, stringsAsFactors=FALSE)

if(!cutoff_colname %in% colnames(cytokines_cutoffs))
  stop("There are no such column with cutoffs!")


# -------------------------------------
# get the isotope and antigen for fcs markers

m <- match(fcs_colnames, cytokines_cutoffs$fcs_colname)

fcs_panel <- data.frame(colnames = fcs_colnames, Isotope = cytokines_cutoffs$Isotope[m], Antigen = cytokines_cutoffs$Antigen[m], stringsAsFactors = FALSE)


# -------------------------------------

### Indeces of observables used for positive-negative analysis

pncols <- which(fcs_panel$colnames %in% cytokines_cutoffs[!is.na(cytokines_cutoffs[, cutoff_colname]), "fcs_colname"])

## Order them in the way as in panel
mm <- match(cytokines_cutoffs[!is.na(cytokines_cutoffs[, cutoff_colname]), "Antigen"], fcs_panel$Antigen[pncols])
pncols <- pncols[mm]


# ------------------------------------------------------------

## clustering results
clust <- read.table(file.path(hmDir, path_clustering), header = TRUE, sep = "\t", as.is = TRUE)
clust <- clust[, "cluster"]

## cluster labels
labels <- read.table(file.path(hmDir, path_clustering_labels), header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))


if(!all(clsubset %in% labels$label))
  stop("Cluster labels are wrong!")


# load marker selection for plotting on the heatmaps
marker_selection <- NULL

if(file.exists(file.path(path_marker_selection))){
  
  marker_selection <- read.table(file.path(path_marker_selection), header = TRUE, sep = "\t", as.is = TRUE)
  marker_selection <- marker_selection[, 1]
  
  if(!all(marker_selection %in% fcs_panel$Antigen[pncols]))
    stop("Marker selection is wrong")
  
}


# ------------------------------------------------------------
# Load tSNE data
# ------------------------------------------------------------

load(file.path(sneDir, path_rtsne_out))

rtsne_data <- read.table(file.path(sneDir, path_rtsne_data), header = TRUE, sep = "\t", as.is = TRUE)


# ------------------------------------------------------------
# Get marker expression
# ------------------------------------------------------------


# arc-sin-h the columns specific 
fcsT <- lapply(fcs, function(u) {
  e <- exprs(u)
  e[ , pncols] <- asinh( e[ , pncols] / 5 )
  exprs(u) <- e
  u
})


# Re-extract the data as list of tables

epn <- lapply(fcsT, function(u) {
  exprs(u)[,pncols]
})

epn <- do.call("rbind",epn)
colnames(epn) <- fcs_panel$Antigen[pncols]


### Normalize the data to 01

epnl <- epn
rng <- apply(epnl,2,quantile,p=c(.01,.99))
for(i in 1:ncol(epnl)) {
  epnl[,i] <- (epnl[,i]-rng[1,i])/(rng[2,i]-rng[1,i])
}
epnl[epnl<0] <- 0
epnl[epnl>1] <- 1

### Keep only the data from the specified clusters

cells2keep <- clust %in% labels[labels$label %in% clsubset, "cluster"]

samp <- rep(names(fcs), sapply(fcs, nrow))[cells2keep]


# Use raw data
if(data2analyse == "raw")
  e <- epn[cells2keep, ]

# Use 01 normalized data
if(data2analyse == "norm")
  e <- epnl[cells2keep, ]


# ------------------------------------------------------------
# Create the bimatrix - TRUE when cells are positive expressed for a given marker
# ------------------------------------------------------------

## get the corresponding cutoffs
mm <- match(colnames(e), cytokines_cutoffs$Antigen)

e_cut <- cytokines_cutoffs[mm, cutoff_colname]
names(e_cut) <- colnames(e)

e_cut

## create the bimatrix
bimatrix <- t(t(e) > e_cut)
bimatrix <- apply(bimatrix, 2, as.numeric)

## remove cells that are always negative
cells2keep2 <- rowSums(bimatrix) > 0

bm <- bimatrix[cells2keep2, ]
samp <- samp[cells2keep2]

# ------------------------------------------------------------
# Upsetr plots
# ------------------------------------------------------------

bidf <- data.frame(bimatrix, row.names = 1:nrow(bimatrix))

pdf(file.path(cyDir, paste0(prefix, "upsetr", suffix, ".pdf")), w = 16, h = 6)
upset(bidf, sets = colnames(bidf), nintersects = 50, order.by = "freq")
dev.off()


# ------------------------------------------------------------
# Cell clustering with FlowSOM
# ------------------------------------------------------------

# Number of clusters
rand_seed <- 1234

# SOM
set.seed(rand_seed)
fsom <- FlowSOM::SOM(bm)


# consensus clustering that is reproducible with seed
data <- fsom$codes
k <- nmetaclusts

pdf(file.path(cyDir, paste0(prefix, prefix_clust, "ConsensusClusterPlus", suffix, ".pdf")), width = 7, height = 7)

results <- ConsensusClusterPlus::ConsensusClusterPlus(t(data),
  maxK = k, reps = 100, pItem = 0.9, pFeature = 1, title = tempdir(),
  plot = NULL, verbose = FALSE, clusterAlg = "hc", distance = "euclidean", seed = rand_seed)

dev.off()


# get cluster ids
fsom_mc <- results[[k]]$consensusClass
clust <- fsom_mc[fsom$mapping[,1]]


which_cells_kept <- seq(length(cells2keep))[cells2keep][cells2keep2]

### Save clustering results
clust_out <- data.frame(cluster = clust, cell_id = which_cells_kept, bm, stringsAsFactors = FALSE)

write.table(clust_out, file = file.path(cyDir, paste0(prefix, prefix_clust, "clustering", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")



### Save cluster frequencies and the mean expression
a <- aggregate(bm, by=list(clust), FUN = mean)

# get cluster frequencies
freq_clust <- table(clust)

clusters_out <- data.frame(cluster = names(freq_clust), label = names(freq_clust), counts = as.numeric(freq_clust), frequencies = as.numeric(freq_clust)/sum(freq_clust) * 100, a[, -1])

write.table(clusters_out, file.path(cyDir, paste0(prefix, prefix_clust, "clusters", suffix, ".xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# ----------------------------
# Plot heatmaps
# ----------------------------

expr <- a[, -1]
rownames(expr) <- a[, 1]

cluster_rows <- hclust(dist(expr), method = "average")
cluster_cols <- hclust(dist(t(expr)), method = "average")

labels_row <- paste0(rownames(expr), " (", round(as.numeric(freq_clust)/sum(freq_clust)*100, 2), "%)") 
labels_col <- colnames(expr)

pheatmap(expr, color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(100), cluster_cols = cluster_cols, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, breaks = seq(from = 0, to = 1, length.out = 101), legend_breaks = seq(from = 0, to = 1, by = 0.2), display_numbers = TRUE, number_color = "black", fontsize_number = 7,  fontsize_row = 10, fontsize_col = 10, fontsize = 7, filename = file.path(cyDir, paste0(prefix, prefix_clust, "pheatmap_all_row_clust", suffix, ".pdf")), width = 10, height = 7)


pheatmap(expr, color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(100), cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, breaks = seq(from = 0, to = 1, length.out = 101), legend_breaks = seq(from = 0, to = 1, by = 0.2), display_numbers = TRUE, number_color = "black", fontsize_number = 7, fontsize_row = 10, fontsize_col = 10, fontsize = 7, filename = file.path(cyDir, paste0(prefix, prefix_clust, "pheatmap_all", suffix, ".pdf")), width = 10, height = 7)



## Plot only the selected markers
if(!is.null(marker_selection)){
  
  expr_sub <- expr[, marker_selection, drop = FALSE]
  labels_col_sub <- colnames(expr_sub)
  
  pheatmap(expr_sub, color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(100), cluster_cols = FALSE, cluster_rows = cluster_rows, labels_col = labels_col_sub, labels_row = labels_row, breaks = seq(from = 0, to = 1, length.out = 101), legend_breaks = seq(from = 0, to = 1, by = 0.2), display_numbers = TRUE, number_color = "black", fontsize_number = 7,  fontsize_row = 10, fontsize_col = 10, fontsize = 7, filename = file.path(cyDir, paste0(prefix, prefix_clust, "pheatmap_s1_row_clust", suffix, ".pdf")), width = 10, height = 7)
  
  
  pheatmap(expr_sub, color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(100), cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col_sub, labels_row = labels_row, breaks = seq(from = 0, to = 1, length.out = 101), legend_breaks = seq(from = 0, to = 1, by = 0.2), display_numbers = TRUE, number_color = "black", fontsize_number = 7, fontsize_row = 10, fontsize_col = 10, fontsize = 7, filename = file.path(cyDir, paste0(prefix, prefix_clust, "pheatmap_s1", suffix, ".pdf")), width = 10, height = 7)
  
  
}




# ------------------------------------------------------------
# tSNE plots
# ------------------------------------------------------------

# get clustering for cells that were used in tSNE and belong to the cluster subset 

which_cells_kept <- seq(length(cells2keep))[cells2keep][cells2keep2]

clust_tsne <- clust[which_cells_kept %in% rtsne_data$cell_index]
e_tsne <- epnl[rtsne_data$cell_index[rtsne_data$cell_index %in% which_cells_kept], ]

ggdf <- data.frame(tSNE1 = rtsne_out$Y[rtsne_data$cell_index %in% which_cells_kept, 1], tSNE2 = rtsne_out$Y[rtsne_data$cell_index %in% which_cells_kept, 2], cluster = clust_tsne, sample = rtsne_data$sample_name[rtsne_data$cell_index %in% which_cells_kept], e_tsne)

## order clusters like in the heatmap
ggdf$cluster <- factor(ggdf$cluster, levels = cluster_rows$labels[cluster_rows$order])

# add group info
mm <- match(ggdf$sample, md$shortname)
ggdf$group <- factor(md$condition[mm])


### Plot of tsne - all cells, all clusters

## facet per group
ggp <- ggplot(ggdf,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
  geom_point(size=1) +
  facet_wrap(~ group) +
  labs(x = "tSNE 1", y="tSNE 2") + 
  xlim(range(rtsne_out$Y[, 1])) +
  ylim(range(rtsne_out$Y[, 2])) +
  theme_bw() +
  theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold")) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

pdf(file.path(cyDir, paste0(prefix, prefix_clust, "tSNE", tsnep_suffix, ".pdf")), width = pdf_width, height = pdf_height)           
print(ggp)
dev.off()


## one plot 
ggp <- ggplot(ggdf,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
  geom_point(size = 1) +
  labs(x = "tSNE 1", y="tSNE 2")+ 
  xlim(range(rtsne_out$Y[, 1])) +
  ylim(range(rtsne_out$Y[, 2])) +
  theme_bw() +
  theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold")) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

pdf(file.path(cyDir, paste0(prefix, prefix_clust, "tSNEone", tsnep_suffix, ".pdf")), width = 9, height = 7)                 
print(ggp)
dev.off()


### Plot of tsne - PD.1 expression

## facet per group
ggp <- ggplot(ggdf,  aes(x = tSNE1, y = tSNE2, color = PD.1)) +
  geom_point(size=1) +
  facet_wrap(~ group) +
  labs(x = "tSNE 1", y="tSNE 2") + 
  xlim(range(rtsne_out$Y[, 1])) +
  ylim(range(rtsne_out$Y[, 2])) +
  scale_colour_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(100)) +
  theme_bw() +
  theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold")) 

pdf(file.path(cyDir, paste0(prefix, "tSNE_PD1", tsnep_suffix, ".pdf")), width = pdf_width, height = pdf_height)                 
print(ggp)
dev.off()


## one plot 
ggp <- ggplot(ggdf,  aes(x = tSNE1, y = tSNE2, color = PD.1)) +
  geom_point(size = 1) +
  labs(x = "tSNE 1", y="tSNE 2")+ 
  xlim(range(rtsne_out$Y[, 1])) +
  ylim(range(rtsne_out$Y[, 2])) +
  scale_colour_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(100)) +
  theme_bw() +
  theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold")) 

pdf(file.path(cyDir, paste0(prefix, "tSNEone_PD1", tsnep_suffix, ".pdf")), width = 9, height = 7)                 
print(ggp)
dev.off()

# ------------------------------------------------------------
# Plots of frequencies for clusters
# ------------------------------------------------------------

cluster <- clust

## Freqs per sample
tabl2 <- table(cluster, samp)

prop <- t(t(tabl2) / colSums(tabl2)) * 100


ggdf <- melt(prop, value.name = "prop")


## add group info
mm <- match(ggdf$samp, md$shortname)
ggdf$group <- factor(md$condition[mm])

## merge base_HD and tx_HD into one level - HD
# new_levels <- levels(ggdf$group)
# new_levels[grep("HD", new_levels)] <- "HD"
# levels(ggdf$group) <- new_levels

## replace _ with \n
levels(ggdf$group) <- gsub("_", "\n", levels(ggdf$group))

## calculate mean and sd for the error bars on the plot
ggds <- ddply(ggdf, .(group, cluster), summarise, mean = mean(prop), sd = sd(prop))


## plot each cluster as a separate page in the pdf file
ggp <- list()

ggdf$cluster <- factor(ggdf$cluster)
clusters <- levels(ggdf$cluster)

for(i in 1:nlevels(ggdf$cluster)){
  
  df <- ggdf[ggdf$cluster == clusters[i], , drop = FALSE]
  ds <- ggds[ggds$cluster == clusters[i], , drop = FALSE]
  
  ggp[[i]] <- ggplot(df, aes(x = group, y = prop)) +
    geom_jitter(size=2.5, shape = 17, width = 0.5, height = 0) +
    geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean, ymax=mean), colour='black', width=0.4) +
    geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean-sd, ymax=mean+sd), colour='black', width=0.25) +
    ggtitle(clusters[i]) +
    theme_bw() +
    ylab("Frequency") +
    xlab("") +
    theme(axis.text.x = element_text(size=12, face="bold"), 
      axis.title.y = element_text(size=12, face="bold"), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(), 
      axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"), 
      axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))
  
}


pdf(file.path(cyDir, paste0(prefix, prefix_clust, "frequencies", suffix, ".pdf")), w=5, h=4, onefile=TRUE)
for(i in seq(length(ggp)))
  print(ggp[[i]])
dev.off()



# # ------------------------------------------------------------
# # Test for frequency differences between groups
# # ------------------------------------------------------------
# 
# ## Freqs per sample
# tabl2 <- table(clust, samp)
# 
# freq2 <- data.frame(cluster = rownames(tabl2), as.data.frame.matrix(tabl2), stringsAsFactors = FALSE)
# 
# prop2 <- data.frame(cluster = rownames(tabl2), as.data.frame.matrix(t(t(tabl2) / colSums(tabl2))) * 100, stringsAsFactors = FALSE)
# 
# ### Save the frequencies and proportions
# write.table(prop2, file=file.path(cyDir, paste0(prefix, "frequencies", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
# write.table(freq2, file=file.path(cyDir,paste0(prefix, "counts", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
# 
# 
# # add more info about samples
# cond_split <- strsplit2(md$condition, "_")
# colnames(cond_split) <- c("day", "response")
# 
# md[, c("day", "response")] <- cond_split
# md$response <- factor(md$response, levels = c("NR", "R", "HD"))
# 
# 
# # -----------------------------
# ### Run two-way ANOVA
# 
# pvs_anova <- t(apply(prop2[, md$shortname], 1, function(y){
#   # y <- prop2[1, md$shortname]
#   
#   ## there must be at least 10 proportions greater than 0
#   if(sum(y > 0) < 10)
#     return(rep(NA, 3))
#   
#   data_tmp <- data.frame(y = as.numeric(y), md[, c("day", "response")])
#   
#   res_tmp <- aov(y ~ day * response, data = data_tmp)
#   sum_tmp <- summary(res_tmp)
#   
#   sum_tmp[[1]][1:3, "Pr(>F)"]
#   
# }))
# 
# movars <- c("day", "response", "day:response")
# colnames(pvs_anova) <- paste0("pval_", movars)
# pvs_anova <- data.frame(pvs_anova)
# 
# 
# ## get adjusted p-values
# 
# adjp_anova <- data.frame(apply(pvs_anova, 2, p.adjust, method = "BH"))
# colnames(adjp_anova) <- paste0("adjp_", movars)
# 
# ## save the results
# pvs_anova_out <- data.frame(cluster = rownames(pvs_anova), pvs_anova, adjp_anova)
# 
# write.table(pvs_anova_out, file=file.path(cyDir, paste0(prefix, "pvs_anova", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
# 
# 
# table(adjp_anova$adjp_response < 0.05)
# table(adjp_anova$adjp_day < 0.05)
# 
# 
# 
# # -----------------------------
# ### Fit a GLM
# 
# 
# pvs_glm <- t(apply(prop2[, md$shortname], 1, function(y){
#   # y <- prop2[1, md$shortname]
#   
#   ## there must be at least 10 proportions greater than 0
#   if(sum(y > 0) < 10)
#     return(rep(NA, 4))
#   
#   data_tmp <- data.frame(y = as.numeric(y), md[, c("day", "response")])
#   
#   res_tmp <- glm(y ~ response + day, data = data_tmp)
#   
#   sum_tmp <- summary(res_tmp)
#   
#   out <- as.numeric(sum_tmp$coefficients[, "Pr(>|t|)"])
#   
#   return(out)
#   
# }))
# 
# data_tmp <- data.frame(y = as.numeric(prop2[1, md$shortname]), md[, c("day", "response")])
# momat <- model.matrix(y ~ response + day, data = data_tmp)
# movars <- colnames(momat)
# 
# colnames(pvs_glm) <- paste0("pval_", movars)
# pvs_glm <- data.frame(pvs_glm)
# 
# 
# ## get adjusted p-values
# 
# adjp_glm <- data.frame(apply(pvs_glm, 2, p.adjust, method = "BH"))
# colnames(adjp_glm) <- paste0("adjp_", movars)
# 
# ## save the results
# pvs_glm_out <- data.frame(cluster = rownames(pvs_glm), pvs_glm, adjp_glm)
# 
# write.table(pvs_glm_out, file=file.path(cyDir, paste0(prefix, "pvs_glm", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
# 
# 
# table(adjp_glm$adjp_responseR < 0.05)
# table(adjp_glm$adjp_daytx < 0.05)
# 
# # -----------------------------
# ### Fit a GLM with inteactions
# 
# 
# pvs_glm <- t(apply(prop2[, md$shortname], 1, function(y){
#   # y <- prop2[1, md$shortname]
#   
#   ## there must be at least 10 proportions greater than 0
#   if(sum(y > 0) < 10)
#     return(rep(NA, 4))
#   
#   data_tmp <- data.frame(y = as.numeric(y), md[, c("day", "response")])
#   
#   res_tmp <- glm(y ~ response + day + response:day, data = data_tmp)
#   
#   sum_tmp <- summary(res_tmp)
#   
#   out <- as.numeric(sum_tmp$coefficients[, "Pr(>|t|)"])
#   
#   return(out)
#   
# }))
# 
# data_tmp <- data.frame(y = as.numeric(prop2[1, md$shortname]), md[, c("day", "response")])
# momat <- model.matrix(y ~ response + day + response:day, data = data_tmp)
# movars <- colnames(momat)
# 
# colnames(pvs_glm) <- paste0("pval_", movars)
# pvs_glm <- data.frame(pvs_glm)
# 
# 
# ## get adjusted p-values
# 
# adjp_glm <- data.frame(apply(pvs_glm, 2, p.adjust, method = "BH"))
# colnames(adjp_glm) <- paste0("adjp_", movars)
# 
# ## save the results
# pvs_glm_out <- data.frame(cluster = rownames(pvs_glm), pvs_glm, adjp_glm)
# 
# write.table(pvs_glm_out, file=file.path(cyDir, paste0(prefix, "pvs_glminter", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
# 
# 
# table(adjp_glm$adjp_responseR < 0.05)
# table(adjp_glm$adjp_daytx < 0.05)
# 




























################################
### 06_cytokines done!
################################