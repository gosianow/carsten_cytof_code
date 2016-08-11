##############################################################################
## <<06_cytokines.R>>

# BioC 3.3
# Created 5 Aug 2016
# Updated 10 Aug 2016

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

##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_02_CD4_merging2'
cytokines_prefix='pnlCD4_pca1_merging_CD4_cytCM_'
path_panel='panel_CD4.xlsx'
path_cytokines_cutoffs='panel_CD4_cytokines_CM_RAW.xlsx'
path_clustering='pnlCD4_pca1_merging_CD4_clustering.xls'
path_clustering_labels='pnlCD4_pca1_merging_CD4_clustering_labels.xls'
clusters2analyse='CM'
cutoff_colname='positive_cutoff_raw'
data2analyse='raw'
cytokines_suffix='_raw'

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
suffix <- cytokines_suffix

if(!data2analyse %in% c("raw", "norm"))
  stop("data2analyse must be 'raw' or 'norm'!")


# ------------------------------------------------------------
# define directories
# ------------------------------------------------------------

fcsDir <- "010_cleanfcs"; if( !file.exists(fcsDir) ) dir.create(fcsDir)
pcaDir <- "020_pcascores"; if( !file.exists(pcaDir) ) dir.create(pcaDir)
hmDir <- "030_heatmaps"; if( !file.exists(hmDir) ) dir.create(hmDir)
cyDir <- "060_cytokines"; if( !file.exists(cyDir) ) dir.create(cyDir)

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

# read metadata
md <- read.xls("metadata.xlsx",stringsAsFactors=FALSE)

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


if(clusters2analyse == "all"){
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
# Plot expression of markers
# ------------------------------------------------------------


distros_wrapper <- function(e, suffix){
  
  df <- data.frame(samp = samp, e)
  dfm <- melt(df, id.var = "samp")
  
  # add group info
  mm <- match(dfm$samp, md$shortname)
  dfm$group <- factor(md$condition[mm])
  
  ggp <- ggplot(dfm, aes(x=value)) + 
    geom_density(adjust = 1, fill = "black", alpha = 0.3) + 
    facet_wrap(~ variable, nrow = 3, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  pdf(file.path(cyDir, paste0(prefix, "distrosmer", suffix,".pdf")), w = ncol(epn), h = 10)
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
    guides(color = guide_legend(nrow = 2)) +
    scale_color_manual(values = color_values)
  
  pdf(file.path(cyDir, paste0(prefix, "distrosgrp", suffix,".pdf")), w = ncol(epn), h = 10)
  print(ggp)
  dev.off()
  
}



distros_wrapper(e = e, suffix = paste0(suffix, "_cyt"))



# ------------------------------------------------------------
# Create the bimatrix - TRUE when cells are positive expressed for a given marker
# ------------------------------------------------------------

## get the corresponding cutoffs
mm <- match(colnames(e), cytokines_cutoffs$Antigen)

e_cut <- cytokines_cutoffs[mm, cutoff_colname]
names(e_cut) <- colnames(e)


## create the bimatrix
bimatrix <- t(t(e) > e_cut)


# ------------------------------------------------------------
# Upsetr plots
# ------------------------------------------------------------

bidf <- data.frame(apply(bimatrix, 2, as.numeric), row.names = 1:nrow(bimatrix))

pdf(file.path(cyDir, paste0(prefix, "upsetr", suffix, ".pdf")), w = 16, h = 6)
upset(bidf, sets = colnames(bidf), nintersects = 50, order.by = "freq")
dev.off()



# ------------------------------------------------------------
# Frequencies of single cytokines
# ------------------------------------------------------------

cs <- colSums(bimatrix)
names(cs) <- paste0("any_", names(cs))
cs <- c("all" = nrow(bimatrix), "none_positive" = sum(rowSums(bimatrix) == 0), "any_positive" = sum(rowSums(bimatrix) > 0), cs)

fc <- data.frame(set = names(cs), counts = cs)

write.table(fc, file=file.path(cyDir, paste0(prefix, "setsize", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")


# ------------------------------------------------------------
# Save frequencies for all combinations
# ------------------------------------------------------------


saving_wrapper <- function(bimatrix, suffix){
  
  ## freqs based on combination names
  bivec <- apply(bimatrix, 1, function(x){
    paste0(colnames(bimatrix)[x], collapse = ".")
  })
  
  ### Freqs per sample
  
  tabl2 <- table(bivec, samp)
  ## skipp the all negative combination ("") for calculating proportions
  tabl2 <- tabl2[rownames(tabl2) != "", ]
  
  ## save counts 
  
  freq2 <- data.frame(combination = rownames(tabl2), as.data.frame.matrix(tabl2), stringsAsFactors = FALSE)
  
  write.table(freq2, file=file.path(cyDir,paste0(prefix, "cytokines_counts", suffix ,".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  
  ## save proportions
  
  prop2 <- data.frame(combination = rownames(tabl2), as.data.frame.matrix(t(t(tabl2) / colSums(tabl2))) * 100, stringsAsFactors = FALSE)
  
  write.table(prop2, file=file.path(cyDir, paste0(prefix, "cytokines_frequencies", suffix ,".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  
  
  ### Freqs for all samples
  
  tabl1 <- table(bivec)
  freq1 <- as.numeric(tabl1)
  names(freq1) <- names(tabl1)
  names(freq1)[names(freq1) == ""] <- "none_positive"
  
  rs <- c("all" = sum(freq1), freq1)
  fc <- data.frame(set = names(rs), counts = rs, frequencies = rs/rs[1]*100)
  
  write.table(fc, file=file.path(cyDir, paste0(prefix, "setsize2", suffix ,".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  
  
  return(NULL)
}



saving_wrapper(bimatrix = bimatrix, suffix = suffix)



# ------------------------------------------------------------
# Test for frequency differences between groups
# ------------------------------------------------------------


## freqs based on combination names
bivec <- apply(bimatrix, 1, function(x){
  paste0(colnames(bimatrix)[x], collapse = ".")
})

### Freqs per sample

tabl2 <- table(bivec, samp)
## skipp the all negative combination ("") for calculating proportions
tabl2 <- tabl2[rownames(tabl2) != "", ]


freq2 <- data.frame(combination = rownames(tabl2), as.data.frame.matrix(tabl2), stringsAsFactors = FALSE)

prop2 <- data.frame(combination = rownames(tabl2), as.data.frame.matrix(t(t(tabl2) / colSums(tabl2))) * 100, stringsAsFactors = FALSE)


# add more info about samples
cond_split <- strsplit2(md$condition, "_")
colnames(cond_split) <- c("day", "response")

md[, c("day", "response")] <- cond_split
md$response <- factor(md$response, levels = c("NR", "R", "HD"))


# -----------------------------
### Run two-way ANOVA

pvs_anova <- t(apply(prop2[, md$shortname], 1, function(y){
  # y <- prop2[1, md$shortname]
  
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


### Plot histograms of p-values
ggp <- ggplot(pvs_anova, aes(x = pval_day)) +
  geom_histogram(breaks = seq(0,1,0.02)) +
  ggtitle("Histogram of p-values") +
  xlab("P-value") +
  ylab("Frequency")

pdf(file.path(cyDir, paste0(prefix, "pvs_anova_day", suffix ,".pdf")), w = 7, h = 7)
print(ggp)
dev.off()

ggp <- ggplot(pvs_anova, aes(x = pval_response)) +
  geom_histogram(breaks = seq(0,1,0.02)) +
  ggtitle("Histogram of p-values") +
  xlab("P-value") +
  ylab("Frequency")

pdf(file.path(cyDir, paste0(prefix, "pvs_anova_response", suffix ,".pdf")), w = 7, h = 7)
print(ggp)
dev.off()



## get adjusted p-values

adjp_anova <- data.frame(apply(pvs_anova, 2, p.adjust, method = "BH"))
colnames(adjp_anova) <- paste0("adjp_", movars)

## save the results
pvs_anova_out <- data.frame(combination = rownames(pvs_anova), pvs_anova, adjp_anova)

write.table(pvs_anova_out, file=file.path(cyDir, paste0(prefix, "pvs_anova", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")


table(adjp_anova$adjp_response < 0.05)
table(adjp_anova$adjp_day < 0.05)



# -----------------------------
### Fit a GLM


pvs_glm <- t(apply(prop2[, md$shortname], 1, function(y){
  # y <- prop2[1, md$shortname]
  
  ## there must be at least 10 proportions greater than 0
  if(sum(y > 0) < 10)
    return(rep(NA, 4))
  
  data_tmp <- data.frame(y = as.numeric(y), md[, c("day", "response")])
  
  res_tmp <- glm(y ~ response + day, data = data_tmp)
  
  sum_tmp <- summary(res_tmp)
  
  out <- as.numeric(sum_tmp$coefficients[, "Pr(>|t|)"])
  
  return(out)
  
}))

data_tmp <- data.frame(y = as.numeric(prop2[1, md$shortname]), md[, c("day", "response")])
momat <- model.matrix(y ~ response + day, data = data_tmp)
movars <- colnames(momat)

colnames(pvs_glm) <- paste0("pval_", movars)
pvs_glm <- data.frame(pvs_glm)


### Plot histograms of p-values
ggp <- ggplot(pvs_glm, aes(x = pval_daytx)) +
  geom_histogram(breaks = seq(0,1,0.02)) +
  ggtitle("Histogram of p-values") +
  xlab("P-value") +
  ylab("Frequency")

pdf(file.path(cyDir, paste0(prefix, "pvs_glm_daytx", suffix ,".pdf")), w = 7, h = 7)
print(ggp)
dev.off()

ggp <- ggplot(pvs_glm, aes(x = pval_responseR)) +
  geom_histogram(breaks = seq(0,1,0.02)) +
  ggtitle("Histogram of p-values") +
  xlab("P-value") +
  ylab("Frequency")

pdf(file.path(cyDir, paste0(prefix, "pvs_glm_responseR", suffix ,".pdf")), w = 7, h = 7)
print(ggp)
dev.off()



## get adjusted p-values

adjp_glm <- data.frame(apply(pvs_glm, 2, p.adjust, method = "BH"))
colnames(adjp_glm) <- paste0("adjp_", movars)

## save the results
pvs_glm_out <- data.frame(combination = rownames(pvs_glm), pvs_glm, adjp_glm)

write.table(pvs_glm_out, file=file.path(cyDir, paste0(prefix, "pvs_glm", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")


table(adjp_glm$adjp_responseR < 0.05)
table(adjp_glm$adjp_daytx < 0.05)





# ------------------------------------------------------------
# Plots of frequencies for combinations of interest
# ------------------------------------------------------------


plotting_wrapper <- function(bimatrix, plot_comb, suffix){
  # suffix = "_teeest"
  
  ## freqs based on combination names
  bivec <- apply(bimatrix, 1, function(x){
    paste0(colnames(bimatrix)[x], collapse = ".")
  })
  
  tabl2 <- table(bivec, samp)
  ## skipp the all negative combination ("") for calculating proportions
  tabl2 <- tabl2[rownames(tabl2) != "", ]

  prop2 <- data.frame(combination = rownames(tabl2), as.data.frame.matrix(t(t(tabl2) / colSums(tabl2))) * 100, stringsAsFactors = FALSE)
  
  ### Prepare data for plotting
  
  ggdf <- melt(prop2[plot_comb,], id.vars = "combination", value.name = "prop", variable.name = "samp")
  
  # order like in top freq
  ggdf$combination <- factor(ggdf$combination, levels = plot_comb)
  
  # add group info
  mm <- match(ggdf$samp, md$shortname)
  ggdf$group <- factor(md$condition[mm])
  
  ## merge base_HD and tx_HD into one level - HD
  # new_levels <- levels(ggdf$group)
  # new_levels[grep("HD", new_levels)] <- "HD"
  # levels(ggdf$group) <- new_levels
  
  ## replace _ with \n
  levels(ggdf$group) <- gsub("_", "\n", levels(ggdf$group))
  
  ## calculate mean and sd for the error bars on the plot
  ggds <- ddply(ggdf, .(group, combination), summarise, mean = mean(prop), sd = sd(prop))
  
  
  ## plot mean frequencies as a barplot
  ggp <- ggplot(data = ggds, aes(x = group, y = mean, fill = combination)) +
    geom_bar(stat="identity") +
    theme_bw() +
    ylab("Mean frequency") +
    xlab("") +
    theme(axis.text.x = element_text(size=14, face="bold"), 
      axis.title.y = element_text(size=14, face="bold"),
      panel.grid.major=element_blank(),
      legend.key = element_blank(),
      legend.title=element_blank())
  
  
  pdf(file.path(cyDir, paste0(prefix, "pbar", suffix ,".pdf")), w = 9, h = 7)
  print(ggp)
  dev.off()
  
  
  # ## plot freqs per sample as points
  # ggp <- ggplot(data = ggdf, aes(x = group, y = prop, fill = combination)) +
  #   geom_point(position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.8), pch = 21, size = 2) +
  #   geom_errorbar(data=ggds, aes(x=group, y=mean, ymin=mean, ymax=mean), position = position_dodge(width = 0.8), colour='black', width=0.4) +
  #   geom_errorbar(data=ggds, aes(x=group, y=mean, ymin=mean-sd, ymax=mean+sd), position = position_dodge(width = 0.8), colour='black', width=0.25) +
  #   geom_vline(xintercept = c(1:(nlevels(ggdf$group)-1) + 0.5), color="grey90") +
  #   theme_bw() +
  #   ylab("Frequency") +
  #   xlab("") +
  #   theme(axis.text.x = element_text(size=14, face="bold"), 
  #     axis.title.y = element_text(size=14, face="bold"),
  #     panel.grid.major=element_blank(),
  #     legend.key = element_blank(),
  #     legend.title=element_blank()) +
  #   guides(fill = guide_legend(override.aes = list(size = 4)))
  # 
  # 
  # pdf(file.path(cyDir, paste0(prefix, "ppoints", suffix ,".pdf")), w = 14, h = 5)
  # print(ggp)
  # dev.off()
  
  
  ## plot freqs per sample as points and facet per combination
  ggp <- ggplot(data = ggdf, aes(x = group, y = prop, fill = combination)) +
    geom_point(position = position_jitter(width = 0.05), pch = 21, size = 2) +
    geom_errorbar(data=ggds, aes(x=group, y=mean, ymin=mean, ymax=mean), colour='black', width=0.4) +
    geom_errorbar(data=ggds, aes(x=group, y=mean, ymin=mean-sd, ymax=mean+sd), colour='black', width=0.25) +
    theme_bw() +
    ylab("Frequency") +
    xlab("") +
    facet_wrap(~ combination, scales = "free_y") +
    theme(axis.text.x = element_text(size=14, face="bold"), 
      axis.title.y = element_text(size=14, face="bold"),
      panel.grid.major=element_blank(),
      legend.position="none",
      legend.key = element_blank(),
      legend.title=element_blank()) +
    guides(fill = guide_legend(override.aes = list(size = 4)))
  
  
  pdf(file.path(cyDir, paste0(prefix, "pfacet", suffix ,".pdf")), w = 14, h = 7)
  print(ggp)
  dev.off()
  
  
  return(NULL)
}





## plot top 5 combinations with diff freqs

## ANOVA

plot_comb <- pvs_anova_out[order(pvs_anova_out$adjp_response, decreasing = FALSE)[1:5], "combination"]

plotting_wrapper(bimatrix = bimatrix, plot_comb = plot_comb, suffix = paste0("_top5_anova_response", suffix))


plot_comb <- pvs_anova_out[order(pvs_anova_out$adjp_day, decreasing = FALSE)[1:5], "combination"]

plotting_wrapper(bimatrix = bimatrix, plot_comb = plot_comb, suffix = paste0("_top5_anova_day", suffix))


## GLM

plot_comb <- pvs_glm_out[order(pvs_glm_out$adjp_responseR, decreasing = FALSE)[1:5], "combination"]

plotting_wrapper(bimatrix = bimatrix, plot_comb = plot_comb, suffix = paste0("_top5_glm_responseR", suffix))


plot_comb <- pvs_glm_out[order(pvs_glm_out$adjp_daytx, decreasing = FALSE)[1:5], "combination"]

plotting_wrapper(bimatrix = bimatrix, plot_comb = plot_comb, suffix = paste0("_top5_glm_daytx", suffix))





























################################
### 06_cytokines done!
################################