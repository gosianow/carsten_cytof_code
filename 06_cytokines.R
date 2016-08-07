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
# cytokines_prefix='pnlCD4_pca1_merging_CD4_'
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

# Raw expression, included observables

df <- data.frame(epn)
dfm <- melt(df)

ggp <- ggplot(dfm, aes(x=value)) + 
  geom_density(adjust=3, fill = "black", alpha = 0.3) + 
  facet_wrap(~ variable, nrow = 3, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf(file.path(cyDir, paste0(prefix, "distrosmer_raw_cyt.pdf")), w = ncol(epn)/2, h = 6)
print(ggp)
dev.off()


# Normalized expression, included observables

df <- data.frame(epnl)
dfm <- melt(df)

ggp <- ggplot(dfm, aes(x=value)) + 
  geom_density(adjust=3, fill = "black", alpha = 0.3) + 
  facet_wrap(~ variable, nrow = 3, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf(file.path(cyDir, paste0(prefix, "distrosmer_norm_cyt.pdf")), w = ncol(epnl)/2, h = 6)
print(ggp)
dev.off()



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
# Plots of frequencies
# ------------------------------------------------------------

samp <- rep( names(fcs), sapply(fcs, nrow) )


plotting_wrapper <- function(bimatrix, suffix){
  
  
  bivec01 <- apply(apply(bimatrix, 2, as.numeric), 1, paste0, collapse = "")
  
  freq <- table(bivec01)
  table(freq)
  
  
  bivec <- apply(bimatrix, 1, function(x){
    paste0(colnames(bimatrix)[x], collapse = ".")
  })
  
  tabl1 <- table(bivec)
  freq1 <- as.numeric(tabl1)
  names(freq1) <- names(tabl1)
  
  tabl2 <- table(bivec, samp)
  freq2 <- data.frame(combination = rownames(tabl2), as.data.frame.matrix(tabl2), stringsAsFactors = FALSE)
  prop2 <- data.frame(combination = rownames(tabl2), as.data.frame.matrix(t(t(tabl2) / colSums(tabl2))), stringsAsFactors = FALSE)
  
  
  ## save frequencies
  write.table(prop2, file=file.path(cyDir, paste0(prefix, "cytokines_frequencies", suffix ,".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  write.table(freq2, file=file.path(cyDir,paste0(prefix, "cytokines_counts", suffix ,".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  
  
  
  ## take top frequent positive combinations
  topnr <- 11
  
  top_freq <- names(freq1)[order(freq1, decreasing = TRUE)[1:topnr]]
  
  ## skipp the "" all negative combination
  top_freq <- top_freq[top_freq != ""]
  
  
  
  ### Plot frequencies
  
  ggdf <- melt(prop2[top_freq,], id.vars = "combination", value.name = "prop", variable.name = "samp")
  
  # order like in top freq
  ggdf$combination <- factor(ggdf$combination, levels = top_freq)
  
  # add group info
  mm <- match(ggdf$samp, md$shortname)
  ggdf$group <- factor(md$condition[mm])
  
  ## merge base_HD and tx_HD into one level - HD
  new_levels <- levels(ggdf$group)
  new_levels[grep("HD", new_levels)] <- "HD"
  levels(ggdf$group) <- new_levels
  
  ## replace _ with \n
  levels(ggdf$group) <- gsub("_", "\n", levels(ggdf$group))
  
  ## calculate mean and sd for the error bars on the plot
  ggds <- ddply(ggdf, .(group, combination), summarise, mean = mean(prop), sd = sd(prop))
  
  
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