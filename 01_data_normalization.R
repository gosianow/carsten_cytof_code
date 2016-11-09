##############################################################################
## <<01_data_normalization.R>>

# BioC 3.3
# Created 24 Aug 2016


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


##############################################################################
# Test arguments
##############################################################################

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
# data_prefix='23_01_'
# data_outdir='010_data'
# path_panel='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_panels/panel1.xlsx'
# path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.xlsx'

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-29_03all_myeloid_merging3'
data_prefix='29mye_03_'
data_outdir='010_data'
path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_29_03all.xlsx'
path_panel='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_panels/panel3.xlsx'

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

setwd(rwd)

prefix <- data_prefix

# ------------------------------------------------------------
# define directories
# ------------------------------------------------------------

fcsDir <- "010_cleanfcs"
outdir <- data_outdir

if(!file.exists(outdir)) 
  dir.create(outdir, recursive = TRUE)


# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

# read metadata
md <- read.xls(path_metadata, stringsAsFactors=FALSE)

## colors per sample
color_values <- md$color
names(color_values) <- md$shortname


# ------------------------------------------------------------
# Load fcs files
# ------------------------------------------------------------

# define FCS file names
f <- file.path(fcsDir, md$filename)
names(f) <- md$shortname

# read raw FCS files in
fcs <- lapply(f, read.FCS)

fcs_colnames <- colnames(fcs[[1]])
fcs_colnames

## Create sample info
samp <- rep(names(fcs), sapply(fcs, nrow))

# ------------------------------------------------------------
# Load panel
# ------------------------------------------------------------

# read panel, pick which columns to use
panel <- read.xls(path_panel, stringsAsFactors=FALSE)

if(!all(c("fcs_colname", "Isotope", "Antigen", "transform") %in% colnames(panel)))
  stop("Wrong columns in panel!!!")


# cols - get fcs columns that are in the panel with transform = 1
cols <- which(fcs_colnames %in% panel$fcs_colname[panel$transform==1])

# get the isotope and antigen for fcs markers
m <- match(fcs_colnames, panel$fcs_colname)

fcs_panel <- data.frame(fcs_colname = fcs_colnames, Isotope = panel$Isotope[m], Antigen = panel$Antigen[m], stringsAsFactors = FALSE)



# --------------------------------------------------------------------------
# Get the marker expression
# --------------------------------------------------------------------------


# arc-sin-h the columns specific
fcsT <- lapply(fcs, function(u) {
  e <- exprs(u)
  e[,cols] <- asinh( e[,cols] / 5 )
  exprs(u) <- e
  u
})


es <- lapply(fcsT, function(u) {
  exprs(u)[, cols]
})

e <- do.call("rbind",es)

# normalize to 0-1
el <- e
rng <- apply(el,2,quantile,p=c(.01,.99))
for(i in 1:ncol(el)) {
  el[,i] <- (el[,i]-rng[1,i])/(rng[2,i]-rng[1,i])
}
el[el<0] <- 0
el[el>1] <- 1


### Save the data

e_out <- data.frame(cell_id = 1:nrow(e), sample_id = samp, e, check.names = FALSE, stringsAsFactors = FALSE)

saveRDS(e_out, file.path(outdir, paste0(prefix, "expr_raw.rds")))



el_out <- data.frame(cell_id = 1:nrow(el), sample_id = samp, el, check.names = FALSE, stringsAsFactors = FALSE)

saveRDS(el_out, file.path(outdir, paste0(prefix, "expr_norm.rds")))





# ------------------------------------------------------------
# Plot expression of markers for (i) pooled/merged data and (ii) strat. per sample
# ------------------------------------------------------------




plotting_wrapper <- function(e, suffix){

  df <- data.frame(samp = samp, e, check.names = FALSE)
  dfm <- melt(df, id.var = "samp")

  ggp <- ggplot(dfm, aes(x=value)) +
    geom_density(adjust = 1, fill = "black", alpha = 0.3) +
    facet_wrap(~ variable, nrow = 4, scales = "free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  pdf(file.path(outdir, paste0(prefix, "distrosmer", suffix,".pdf")), w = ncol(e)*2/3, h = 10)
  print(ggp)
  dev.off()


  ggp <- ggplot(dfm, aes(x=value, color = samp)) +
    geom_density(adjust = 1) +
    facet_wrap(~ variable, nrow = 4, scales = "free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title = element_blank(), legend.position = "bottom") +
    guides(color = guide_legend(nrow = 2)) +
    scale_color_manual(values = color_values)

  pdf(file.path(outdir, paste0(prefix, "distrosgrp", suffix,".pdf")), w = ncol(e)*2/3, h = 11)
  print(ggp)
  dev.off()


  return(NULL)

}




## Expression of raw data

colnames(e) <- fcs_panel$Antigen[cols]

plotting_wrapper(e = e, suffix = "_raw")


## Expression of normalized data

colnames(el) <- fcs_panel$Antigen[cols]

plotting_wrapper(e = el, suffix = "_norm")




# ------------------------------------------------------------
# Plot number of cells per sample
# ------------------------------------------------------------

samp <- factor(samp, levels = md$shortname)

cell_table <- table(samp)

cell_counts <- rep(0, nrow(md))
names(cell_counts) <- md$shortname

cell_counts[names(cell_table)] <- as.numeric(cell_table)


ggdf <- data.frame(sample_id = factor(names(cell_counts), levels = md$shortname), cell_counts = cell_counts)


ggp <- ggplot(ggdf, aes(x = sample_id, y = cell_counts, fill = sample_id)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = cell_counts), hjust=0.5, vjust=-0.5, size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none") +
  scale_fill_manual(values = color_values, drop=FALSE) +
  scale_x_discrete(drop=FALSE)


pdf(file.path(outdir, paste0(prefix, "cell_counter.pdf")), w = nlevels(ggdf$sample_id)/3 + 2, h = 5)
print(ggp)
dev.off()







sessionInfo()













################################
### 01_data_normalization.R done!
################################