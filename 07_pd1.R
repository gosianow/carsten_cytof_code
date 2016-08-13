##############################################################################
## <<07_pd1.R>>

# BioC 3.3
# Created 11 Aug 2016
# Updated 11 Aug 2016

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

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_02_CD4_merging2'
# prefix='pnlCD4_pca1_merging_CD4_pd1Tmem_'
# path_panel='panel_CD4.xlsx'
# path_cytokines_cutoffs='panel_CD4_cytokines_CM_RAW.xlsx'
# path_clustering='pnlCD4_pca1_merging_CD4_clustering.xls'
# path_clustering_labels='pnlCD4_pca1_merging_CD4_clustering_labels.xls'
# clusters2analyse=c('CM','EM','TE')
# cutoff_colname='positive_cutoff_raw'
# data2analyse='raw'
# suffix='_20cl_raw'
# nmetaclusts=20
# path_metadata

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

if(!data2analyse %in% c("raw", "norm"))
  stop("data2analyse must be 'raw' or 'norm'!")


# ------------------------------------------------------------
# define directories
# ------------------------------------------------------------

fcsDir <- "010_cleanfcs"; if( !file.exists(fcsDir) ) dir.create(fcsDir)
pcaDir <- "020_pcascores"; if( !file.exists(pcaDir) ) dir.create(pcaDir)
hmDir <- "030_heatmaps"; if( !file.exists(hmDir) ) dir.create(hmDir)
pd1Dir <- "070_pd1"; if( !file.exists(pd1Dir) ) dir.create(pd1Dir)

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


# read panel, pick which columns to use
panel <- read.xls(path_panel,stringsAsFactors=FALSE)

# get isotope mass of columns in fcs files.. to match against the panel
panel_mass <- as.numeric(gsub("[[:alpha:]]", "", colnames(fcs[[1]])))

# cols - get fcs columns that are in the panel with transform = 1
cols <- which(panel_mass %in% panel$Isotope[panel$transform==1])

# Antigen - get the antigen name
m <- match(panel_mass, panel$Isotope)

fcs_panel <- data.frame(colnames = colnames(fcs[[1]]), Isotope = panel_mass, cols = panel_mass %in% panel$Isotope[panel$transform==1], Antigen = panel$Antigen[m], stringsAsFactors = FALSE)

fcs_panel$Antigen[is.na(fcs_panel$Antigen)] <- ""


# arc-sin-h the columns specific 
fcsT <- lapply(fcs, function(u) {
  e <- exprs(u)
  e[,cols] <- asinh( e[,cols] / 5 )
  exprs(u) <- e
  u
})


# ------------------------

# add more info about samples
cond_split <- strsplit2(md$condition, "_")
colnames(cond_split) <- c("day", "response")

md[, c("day", "response")] <- cond_split
md$response <- factor(md$response, levels = c("NR", "R", "HD"))


# ------------------------------------------------------------
# Load more data
# ------------------------------------------------------------

## positive cutoffs for cytokines
cytokines_cutoffs <- read.xls(path_cytokines_cutoffs, stringsAsFactors=FALSE)

if(!cutoff_colname %in% colnames(cytokines_cutoffs))
  stop("There are no such column with cutoffs!")

cytokines_cutoffs <- merge(panel, cytokines_cutoffs[, c("Isotope", cutoff_colname)], by = "Isotope", all.x = TRUE, sort = FALSE)

## clustering results
clust <- read.table(file.path(hmDir, path_clustering), header = TRUE, sep = "\t", as.is = TRUE)
clust <- clust[, 1]

## cluster labels
labels <- read.table(file.path(hmDir, path_clustering_labels), header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))


if(any(clusters2analyse == "all")){
  clusters2analyse <- labels$label
}else{
  if(!all(clusters2analyse %in% labels$label))
    stop("Cluster labels are wrong!")
}

# ------------------------------------------------------------

### Indeces of observables used for positive-negative analysis

pncols <- which(fcs_panel$Isotope %in% cytokines_cutoffs[!is.na(cytokines_cutoffs[, cutoff_colname]), "Isotope"])

## Order them in the way as in panel
mm <- match(cytokines_cutoffs[!is.na(cytokines_cutoffs[, cutoff_colname]), "Antigen"], fcs_panel$Antigen[pncols])
pncols <- pncols[mm]

if(!all(pncols %in% cols))
  stop("Cytokine positive cutoffs can be defined only for the transform=1 observables!")

### Index of PD-1

pd1col <- which(fcs_panel$Antigen == "PD-1")

if(!pd1col %in% pncols)
  stop("There is no positive cutoff defined for PD-1")


# ------------------------------------------------------------
# Get marker expression
# ------------------------------------------------------------

### Re-extract the data as list of tables

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

cells2keep <- clust %in% labels[labels$label %in% clusters2analyse, "cluster"]

epn <- epn[cells2keep, ]
epnl <- epnl[cells2keep, ]


# Use raw data
if(data2analyse == "raw")
  e <- epn

# Use 01 normalized data
if(data2analyse == "norm")
  e <- epnl

samp <- rep(names(fcs), sapply(fcs, nrow))[cells2keep]



# ------------------------------------------------------------
# Plot distribution of observables
# ------------------------------------------------------------


df <- data.frame(samp = samp, e)
dfm <- melt(df, id.var = "samp")

# add group info
mm <- match(dfm$samp, md$shortname)
dfm$group <- factor(md$condition[mm])

ggp <- ggplot(dfm, aes(x=value)) + 
  geom_density(adjust = 1, fill = "black", alpha = 0.3) + 
  facet_wrap(~ variable, nrow = 3, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf(file.path(pd1Dir, paste0(prefix, "distrosmer", suffix,".pdf")), w = ncol(e), h = 10)
print(ggp)
dev.off()


## create colors per group not per sample
mm <- match(levels(dfm$samp), md$shortname)
groups <- factor(md$condition[mm])
color_values <- colorRampPalette(brewer.pal(12,"Paired"))(12)[c(1,3,5,2,4,6)]
color_values <- color_values[as.numeric(groups)]
names(color_values) <- levels(dfm$samp)

ggp <- ggplot(dfm, aes(x=value, color = samp)) + 
  geom_density(adjust = 1) + 
  facet_wrap(~ variable, nrow = 3, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title = element_blank(), legend.position = "bottom") +
  guides(color = guide_legend(nrow = 3, byrow = TRUE)) +
  scale_color_manual(values = color_values)

pdf(file.path(pd1Dir, paste0(prefix, "distrosgrp", suffix,".pdf")), w = ncol(e), h = 10)
print(ggp)
dev.off()



# ------------------------------------------------------------
# Get frequencies of positive PD-1
# ------------------------------------------------------------

epd1 <- e[, "PD-1"]

pd1cut <- cytokines_cutoffs[cytokines_cutoffs$Antigen == "PD-1", cutoff_colname]

bivec <- epd1 > pd1cut

## Freqs of positive PD-1 per sample
tabl1 <- table(bivec, samp)

freq1 <- data.frame(cluster = rownames(tabl1), as.data.frame.matrix(tabl1), stringsAsFactors = FALSE)

prop1 <- data.frame(cluster = rownames(tabl1), as.data.frame.matrix(t(t(tabl1) / colSums(tabl1))) * 100, stringsAsFactors = FALSE)


### Save the frequencies and proportions
write.table(prop1, file=file.path(pd1Dir, paste0(prefix, "frequencies_pd1pos", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
write.table(freq1, file=file.path(pd1Dir,paste0(prefix, "counts_pd1pos", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")


# ------------------------------------------------------------
# Plot PD-1 pos frequencies
# ------------------------------------------------------------

ggdf <- melt(prop1, value.name = "prop", variable.name = "samp", id.vars = "cluster")

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
    geom_jitter(size=2.5, shape = 17, width = 0.5) +
    geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean, ymax=mean), colour='black', width=0.4) +
    geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean-sd, ymax=mean+sd), colour='black', width=0.25) +
    ggtitle(clusters[i]) +
    theme_bw() +
    ylab("Frequency") +
    xlab("") +
    # ylim(c(0, max(df$prop))) +
    theme(axis.text.x = element_text(size=12, face="bold"), 
      axis.title.y = element_text(size=12, face="bold"), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(), 
      axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"), 
      axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))
  
}


pdf(file.path(pd1Dir, paste0(prefix, "frequencies_pd1pos", suffix, ".pdf")), w=5, h=4, onefile=TRUE)
for(i in seq(length(ggp)))
  print(ggp[[i]])
dev.off()



# ------------------------------------------------------------
# Test for frequency differences between groups
# ------------------------------------------------------------

data <- prop1

# -----------------------------
### Run two-way ANOVA
# -----------------------------

pvs_anova <- t(apply(data[, md$shortname], 1, function(y){
  # y <- data[1, md$shortname]
  
  ## there must be at least 10 proportions greater than 0
  if(sum(y > 0) < 10)
    return(rep(NA, 3))
  
  data_tmp <- data.frame(y = as.numeric(y), md[, c("day", "response")])
  
  res_tmp <- aov(y ~ day * response, data = data_tmp)
  sum_tmp <- summary(res_tmp)
  
  sum_tmp[[1]][1:3, "Pr(>F)"]
  
}))

movars <- c("day", "response", "day:response")
colnames(pvs_anova) <- paste0("pval_", movars)
pvs_anova <- data.frame(pvs_anova)


## get adjusted p-values

adjp_anova <- data.frame(apply(pvs_anova, 2, p.adjust, method = "BH"))
colnames(adjp_anova) <- paste0("adjp_", movars)

## save the results
pvs_anova_out <- data.frame(cluster = rownames(pvs_anova), pvs_anova, adjp_anova)

write.table(pvs_anova_out, file=file.path(pd1Dir, paste0(prefix, "pvs_anova", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")


table(adjp_anova$adjp_response < 0.05)
table(adjp_anova$adjp_day < 0.05)



# -----------------------------
### Fit a GLM
# -----------------------------

pvs_glm <- t(apply(data[, md$shortname], 1, function(y){
  # y <- data[1, md$shortname]
  
  ## there must be at least 10 proportions greater than 0
  if(sum(y > 0) < 10)
    return(rep(NA, 4))
  
  data_tmp <- data.frame(y = as.numeric(y), md[, c("day", "response")])
  
  res_tmp <- glm(y ~ response + day, data = data_tmp)
  
  sum_tmp <- summary(res_tmp)
  
  out <- as.numeric(sum_tmp$coefficients[, "Pr(>|t|)"])
  
  return(out)
  
}))

data_tmp <- data.frame(y = as.numeric(data[1, md$shortname]), md[, c("day", "response")])
momat <- model.matrix(y ~ response + day, data = data_tmp)
movars <- colnames(momat)

colnames(pvs_glm) <- paste0("pval_", movars)
pvs_glm <- data.frame(pvs_glm)


## get adjusted p-values

adjp_glm <- data.frame(apply(pvs_glm, 2, p.adjust, method = "BH"))
colnames(adjp_glm) <- paste0("adjp_", movars)

## save the results
pvs_glm_out <- data.frame(cluster = rownames(pvs_glm), pvs_glm, adjp_glm)

write.table(pvs_glm_out, file=file.path(pd1Dir, paste0(prefix, "pvs_glm", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")


table(adjp_glm$adjp_responseR < 0.05)
table(adjp_glm$adjp_daytx < 0.05)

# -----------------------------
### Fit a GLM with inteactions
# -----------------------------

pvs_glm <- t(apply(data[, md$shortname], 1, function(y){
  # y <- data[1, md$shortname]
  
  ## there must be at least 10 proportions greater than 0
  if(sum(y > 0) < 10)
    return(rep(NA, 4))
  
  data_tmp <- data.frame(y = as.numeric(y), md[, c("day", "response")])
  
  res_tmp <- glm(y ~ response + day + response:day, data = data_tmp)
  
  sum_tmp <- summary(res_tmp)
  
  out <- as.numeric(sum_tmp$coefficients[, "Pr(>|t|)"])
  
  return(out)
  
}))

data_tmp <- data.frame(y = as.numeric(data[1, md$shortname]), md[, c("day", "response")])
momat <- model.matrix(y ~ response + day + response:day, data = data_tmp)
movars <- colnames(momat)

colnames(pvs_glm) <- paste0("pval_", movars)
pvs_glm <- data.frame(pvs_glm)


## get adjusted p-values

adjp_glm <- data.frame(apply(pvs_glm, 2, p.adjust, method = "BH"))
colnames(adjp_glm) <- paste0("adjp_", movars)

## save the results
pvs_glm_out <- data.frame(cluster = rownames(pvs_glm), pvs_glm, adjp_glm)

write.table(pvs_glm_out, file=file.path(pd1Dir, paste0(prefix, "pvs_glminter", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")


table(adjp_glm$adjp_responseR < 0.05)
table(adjp_glm$adjp_daytx < 0.05)




# ------------------------------------------------------------
# Continue with stratification if there are other observables to consider positive
# ------------------------------------------------------------

if(ncol(e) > 1){
  
  ### Treat PD-1 positive cells as 100% 
  
  e <- e[epd1 > pd1cut, ]
  samp <- samp[epd1 > pd1cut]
  
  # ------------------------------------------------------------
  # Create the bimatrix - TRUE when cells are positive expressed for a given marker
  # ------------------------------------------------------------
  
  ## get the corresponding cutoffs
  mm <- match(colnames(e), cytokines_cutoffs$Antigen)
  
  cytcut <- cytokines_cutoffs[mm, cutoff_colname]
  names(cytcut) <- colnames(e)
  
  ## create the bimatrix
  bimatrix <- t(t(e) > cytcut)
  bimatrix <- apply(bimatrix, 2, as.numeric)
  
  ## remove cells that are always negative
  bm <- bimatrix[rowSums(bimatrix) > 0, ]
  samp <- samp[rowSums(bimatrix) > 0]
  
  
  # ------------------------------------------------------------
  # Upsetr plots
  # ------------------------------------------------------------
  
  bidf <- data.frame(bimatrix, row.names = 1:nrow(bimatrix))
  
  pdf(file.path(pd1Dir, paste0(prefix, "upsetr", suffix, ".pdf")), w = 16, h = 6)
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
  
  pdf(file.path(pd1Dir, paste0(prefix, "ConsensusClusterPlus", suffix, ".pdf")), width = 7, height = 7)
  
  results <- ConsensusClusterPlus::ConsensusClusterPlus(t(data),
    maxK = k, reps = 100, pItem = 0.9, pFeature = 1, title = tempdir(),
    plot = NULL, verbose = FALSE, clusterAlg = "hc", distance = "euclidean", seed = rand_seed)
  
  dev.off()
  
  # get cluster ids
  fsom_mc <- results[[k]]$consensusClass
  clust <- fsom_mc[fsom$mapping[,1]]

  ### Save clustering results
  clust_out <- data.frame(cluster = clust, bm, stringsAsFactors = FALSE)
  
  write.table(clust_out, file = file.path(pd1Dir, paste0(prefix, "clustering", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  
  
  
  ### Save cluster frequencies and the median expression
  
  a <- aggregate(bm, by=list(clust), FUN = mean)
  
  # get cluster frequencies
  freq_clust <- table(clust)
  
  clusters_out <- data.frame(cluster = names(freq_clust), label = names(freq_clust), counts = as.numeric(freq_clust), frequencies = as.numeric(freq_clust)/sum(freq_clust) * 100, a[, -1])
  
  write.table(clusters_out, file.path(pd1Dir, paste0(prefix, "clusters", suffix, ".xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  
  # ----------------------------
  # Plot heatmaps
  # ----------------------------
  
  expr <- a[, -1]
  rownames(expr) <- a[, 1]
  
  cluster_rows <- hclust(dist(expr), method = "average")
  cluster_cols <- hclust(dist(t(expr)), method = "average")
  
  labels_row <- paste0(rownames(expr), " (", round(clusters_out$frequencies, 2), "%)")
  labels_col <- colnames(expr)
  
  pheatmap(expr, color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(100), cluster_cols = cluster_cols, cluster_rows = cluster_rows, labels_col = labels_col, labels_row = labels_row, breaks = seq(from = 0, to = 1, length.out = 101), legend_breaks = seq(from = 0, to = 1, by = 0.2), display_numbers = TRUE, number_color = "black", fontsize_number = 7,  fontsize_row = 10, fontsize_col = 10, fontsize = 7, filename = file.path(pd1Dir, paste0(prefix, "pheatmap_row_clust", suffix, ".pdf")), width = 10, height = 7)
  
  
  pheatmap(expr, color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(100), cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, breaks = seq(from = 0, to = 1, length.out = 101), legend_breaks = seq(from = 0, to = 1, by = 0.2), display_numbers = TRUE, number_color = "black", fontsize_number = 7, fontsize_row = 10, fontsize_col = 10, fontsize = 7, filename = file.path(pd1Dir, paste0(prefix, "pheatmap", suffix, ".pdf")), width = 10, height = 7)
  
  
}

























################################
### 07_pd1 done!
################################