##############################################################################
## <<01_pcascores.R>>

# BioC 3.3
# Created 27 July 2016
# Updated 24 Aug 2016

##############################################################################
Sys.time()
##############################################################################

# Load packages
library(gdata)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2)
library(coop)
library(RColorBrewer)

##############################################################################
# Test arguments
##############################################################################

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
# pcas_prefix='23_01_'
# pcas_outdir='020_pcascores'
# path_data='010_data/23_01_expr_raw.rds'
# path_panel='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_panels/panel1.xlsx'
# path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.xlsx'

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

prefix <- pcas_prefix

# ------------------------------------------------------------
# define directories
# ------------------------------------------------------------

pcaDir <- pcas_outdir
if( !file.exists(pcaDir) ) dir.create(pcaDir)

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

# read expression data
if(grepl(".txt", path_data)){
  expr <- read.table(path_data, header = TRUE, sep = "\t", as.is = TRUE)
}
if(grepl(".rds", path_data)){
  expr <- readRDS(path_data)
}

samp <- expr[, "sample_id"]
fcs_colnames <- colnames(expr)[!grepl("cell_id|sample_id", colnames(expr))]
e <- expr[, fcs_colnames]


# read metadata
md <- read.xls(path_metadata, stringsAsFactors=FALSE)


# ------------------------------------------------------------
# Load panel
# ------------------------------------------------------------

# read panel, pick which columns to use
panel <- read.xls(path_panel, stringsAsFactors=FALSE)

if(!all(c("fcs_colname", "Isotope", "Antigen", "transform") %in% colnames(panel)))
  stop("Wrong columns in panel!!!")


# get the isotope and antigen for fcs markers
m <- match(fcs_colnames, panel$fcs_colname)

fcs_panel <- data.frame(fcs_colname = fcs_colnames, Isotope = panel$Isotope[m], Antigen = panel$Antigen[m], stringsAsFactors = FALSE)


# --------------------------------------------------------------------------
# run Levine et al. 2015 marker scoring,
# plot marginal distributions
# --------------------------------------------------------------------------

doPRINCOMP <- function(z, ncomp=3) {
  # z = es[[1]]; ncomp=3
  
  ## Score by Levine
  pr <- prcomp(z, center = TRUE, scale. = FALSE)  # default is to center but not scale
  
  score <- rowSums (outer( rep(1, ncol(z)), pr$sdev[1:ncomp]^2 ) * abs(pr$rotation[,1:ncomp]) )
  
  return(score)
  
}


# split expression per sample
es <- split(e, samp)

prs <- sapply(es, doPRINCOMP)
rmprs <- rowMeans(prs)

prs <- data.frame(mass = fcs_panel$fcs_colname, marker = fcs_panel$Antigen, avg_score = rmprs, round(prs, 4))

### Save ordered PCA scores
o <- order(rmprs, decreasing=TRUE)
prs <- prs[o,]

write.table(prs, file = file.path(pcaDir, paste0(prefix, "princompscore_by_sample.xls")), sep="\t", row.names=FALSE, quote=FALSE)




## plot PCA scores

ggp1 <- ggplot(prs, aes(x = marker, y = avg_score)) +
  geom_bar(stat="identity") +
  theme_bw() +
  ylab("Average PCA score") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12))


## plot PCA scores for ordered markers

prs$marker <- factor(prs$marker, levels = prs$marker)

ggp2 <- ggplot(prs, aes(x = marker, y = avg_score)) +
  geom_bar(stat="identity") +
  theme_bw() +
  ylab("Average PCA score") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12))


pdf(file.path(pcaDir, paste0(prefix, "princompscore_average.pdf")), width = 10, height = 7)
print(ggp1)
print(ggp2)
dev.off()








### Plot chanel distributions

pdf(file.path(pcaDir, paste0(prefix, "channel_distributions.pdf")))

ttl <- fcs_panel$Antigen

colors <- as.numeric(as.factor(md$condition))

# for every channel
for(i in 1:nrow(fcs_panel)){
  
  ii <- o[i]
  dens <- lapply(es, function(u) density( u[,ii] ))
  my <- max(sapply(dens, function(u) max(u$y)))
  
  maxx <- max(sapply(dens, function(u) max(u$x)))
  minx <- min(sapply(dens, function(u) min(u$x)))
  
  # for every sample
  k <- prs$marker==ttl[ii]
  for(s in 1:length(es)) {
    if(s==1)
      plot( dens[[s]], xlim = c(minx, maxx), ylim = c(0, my), xlab="", lwd=3, col=colors[s], 
        main=paste(fcs_panel$Isotope[ii],"/", as.character(ttl[ii]),"/", round(prs$avg_score[k],2)))
    else
      lines( dens[[s]], lwd=3, col=colors[s])
    nm <- paste0(names(es)," (", round(prs[k, names(es)], 2) ,")")
    legend("topright", nm, lwd=1, col = colors, cex = 0.7)
  }
  
}

dev.off()











sessionInfo()













################################
### 01_pcascores done!
################################