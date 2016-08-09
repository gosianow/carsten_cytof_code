##############################################################################
## <<06_cytokines.R>>

# BioC 3.3
# Created 5 Aug 2016

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

##############################################################################
# Test arguments
##############################################################################

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_02_CD4_merging2'
# cytokines_prefix='pnlCD4_pca1_'
# path_panel='panel_CD4.xlsx'
# path_cytokines_cutoffs='panel_CD4_cytokines.xlsx'


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


# ------------------------------------------------------------
# define directories
# ------------------------------------------------------------

fcsDir <- "010_cleanfcs"; if( !file.exists(fcsDir) ) dir.create(fcsDir)
pcaDir <- "020_pcascores"; if( !file.exists(pcaDir) ) dir.create(pcaDir)
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

cytokines_cutoffs <- merge(panel, cytokines_cutoffs[, c("Isotope", "positive_cutoff_raw", "positive_cutoff_norm")], by = "Isotope", all.x = TRUE, sort = FALSE)


# ------------------------------------------------------------

### Indeces of observables used for positive-negative analysis

pncols <- which(fcs_panel$Isotope %in% cytokines_cutoffs[!is.na(cytokines_cutoffs$positive_cutoff_raw), "Isotope"])

## Order them in the way as in panel
mm <- match(cytokines_cutoffs[!is.na(cytokines_cutoffs$positive_cutoff_raw), "Antigen"], fcs_panel$Antigen[pncols])
pncols <- pncols[mm]


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



# ------------------------------------------------------------
# Plot expression of markers
# ------------------------------------------------------------


samp <- rep(names(fcs), sapply(fcs, nrow))

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




# Raw expression, included observables

distros_wrapper(e = epn, suffix = "_raw_cyt")

# Normalized expression, included observables

distros_wrapper(e = epnl, suffix = "_norm_cyt")




# ------------------------------------------------------------
# Create the bimatrix - TRUE when cells are positive expressed for a given marker
# ------------------------------------------------------------

## get the corresponding cutoffs
mm <- match(colnames(epnl), cytokines_cutoffs$Antigen)

epn_cut <- cytokines_cutoffs[mm, "positive_cutoff_raw"]
names(epn_cut) <- colnames(epn)

epnl_cut <- cytokines_cutoffs[mm, "positive_cutoff_norm"]
names(epnl_cut) <- colnames(epnl)


## create the bimatrix

bimatrix <- t(t(epn) > epn_cut)
bimatrixl <- t(t(epnl) > epnl_cut)


# ------------------------------------------------------------
# Upsetr plots
# ------------------------------------------------------------

bidf <- data.frame(apply(bimatrix, 2, as.numeric), row.names = 1:nrow(bimatrix))
bidfl <- data.frame(apply(bimatrixl, 2, as.numeric), row.names = 1:nrow(bimatrixl))


pdf(file.path(cyDir, paste0(prefix, "upsetr_raw.pdf")), w = 16, h = 6)
upset(bidf, sets = colnames(bidf), nintersects = 50, order.by = "freq")
dev.off()


pdf(file.path(cyDir, paste0(prefix, "upsetr_norm.pdf")), w = 16, h = 6)
upset(bidfl, sets = colnames(bidf), nintersects = 50, order.by = "freq")
dev.off()


# ------------------------------------------------------------
# Frequencies of single cytokines
# ------------------------------------------------------------

## raw data
cs <- colSums(bimatrix)
names(cs) <- paste0("any_", names(cs))
cs <- c("all" = nrow(bimatrix), "none_positive" = sum(rowSums(bimatrix) == 0), "any_positive" = sum(rowSums(bimatrix) > 0), cs)

fc <- data.frame(set = names(cs), counts = cs)

write.table(fc, file=file.path(cyDir, paste0(prefix, "setsize_raw.xls")), row.names=FALSE, quote=FALSE, sep="\t")


## norm data
cs <- colSums(bimatrixl)
names(cs) <- paste0("any_", names(cs))
cs <- c("all" = nrow(bimatrixl), "none_positive" = sum(rowSums(bimatrixl) == 0), "any_positive" = sum(rowSums(bimatrixl) > 0), cs)

fc <- data.frame(set = names(cs), counts = cs)

write.table(fc, file=file.path(cyDir, paste0(prefix, "setsize_norm.xls")), row.names=FALSE, quote=FALSE, sep="\t")


# ------------------------------------------------------------
# Save frequencies for all combinations
# ------------------------------------------------------------

samp <- rep( names(fcs), sapply(fcs, nrow) )


saving_wrapper <- function(bimatrix, suffix){
  # suffix = "_raw"
  
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



saving_wrapper(bimatrix = bimatrix, suffix = "_raw")


saving_wrapper(bimatrix = bimatrixl, suffix = "_norm")




# ------------------------------------------------------------
# Test for frequency differences between groups
# ------------------------------------------------------------









# ------------------------------------------------------------
# Plots of frequencies for combinations of interest
# ------------------------------------------------------------

samp <- rep( names(fcs), sapply(fcs, nrow) )


plotting_wrapper <- function(bimatrix, suffix){
  # suffix = "_raw"

  ## freqs based on combination names
  bivec <- apply(bimatrix, 1, function(x){
    paste0(colnames(bimatrix)[x], collapse = ".")
  })
  

  tabl2 <- table(bivec, samp)
  ## skipp the all negative combination ("") for calculating proportions
  tabl2 <- tabl2[rownames(tabl2) != "", ]
  

  prop2 <- data.frame(combination = rownames(tabl2), as.data.frame.matrix(t(t(tabl2) / colSums(tabl2))) * 100, stringsAsFactors = FALSE)
  
  
  ## take 10 top frequent positive combinations to plot
  tabl1 <- table(bivec)
  freq1 <- as.numeric(tabl1)
  names(freq1) <- names(tabl1)
  
  topnr <- 11
  top_freq <- names(freq1)[order(freq1, decreasing = TRUE)[1:topnr]]
  
  ## skipp the "" all negative combination
  plot_comb <- top_freq[top_freq != ""]
  
  
  
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
  
  
  ## plot mean frequencies
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
  
  
  ## plot freqs per sample
  ggp <- ggplot(data = ggdf, aes(x = group, y = prop, fill = combination)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.8), pch = 21, size = 2) +
    geom_errorbar(data=ggds, aes(x=group, y=mean, ymin=mean, ymax=mean), position = position_dodge(width = 0.8), colour='black', width=0.4) +
    geom_errorbar(data=ggds, aes(x=group, y=mean, ymin=mean-sd, ymax=mean+sd), position = position_dodge(width = 0.8), colour='black', width=0.25) +
    geom_vline(xintercept = c(1:(nlevels(ggdf$group)-1) + 0.5), color="grey90") +
    theme_bw() +
    ylab("Frequency") +
    xlab("") +
    theme(axis.text.x = element_text(size=14, face="bold"), 
      axis.title.y = element_text(size=14, face="bold"),
      panel.grid.major=element_blank(),
      legend.key = element_blank(),
      legend.title=element_blank()) +
    guides(fill = guide_legend(override.aes = list(size = 4)))
  
  
  pdf(file.path(cyDir, paste0(prefix, "ppoints", suffix ,".pdf")), w = 14, h = 5)
  print(ggp)
  dev.off()
  
  return(NULL)
}


plotting_wrapper(bimatrix = bimatrix, suffix = "_raw")


plotting_wrapper(bimatrix = bimatrixl, suffix = "_norm")
































################################
### 06_cytokines done!
################################