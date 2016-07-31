##############################################################################
## <<01_pcascores.R>>

# BioC 3.3
# Created 27 July 2016


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

##############################################################################
# Test arguments
##############################################################################

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-29_02'
# pcas_prefix=''
# path_panel='panel.xlsx'

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


# ------------------------------------------------------------
# define directories
# ------------------------------------------------------------

fcsDir <- "010_cleanfcs"; if( !file.exists(fcsDir) ) dir.create(fcsDir)
pcaDir <- "020_pcascores"; if( !file.exists(pcaDir) ) dir.create(pcaDir)


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
panel <- read.xls(path_panel, stringsAsFactors=FALSE)


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



# --------------------------------------------------------------------------
# run Levine et al. 2015 marker scoring,
# plot marginal distributions
# --------------------------------------------------------------------------

doPRINCOMP <- function(z, ncomp=3) {
  cat(".")
  pr <- prcomp(z)  # default is to center but not scale
  score <- outer( rep(1,ncol(z)), pr$sdev[1:ncomp]^2 ) * abs(pr$rotation[,1:ncomp])
  rowSums(score)
}


# extract columns of interest to list of tables
es <- lapply(fcsT, function(u) {
  exprs(u)[,cols]
})


prs <- sapply(es, doPRINCOMP)
rmprs <- rowMeans(prs)
o <- order(rmprs, decreasing=TRUE)

prs <- data.frame(mass = rownames(prs), marker = fcs_panel$Antigen[cols], avg_score=rmprs, round(prs,3))[o,]

write.table(prs, file = file.path(pcaDir, paste0(pcas_prefix, "princompscore_by_sample.xls")), sep="\t", row.names=FALSE, quote=FALSE)




### Plot chanel distributions

pdf(file.path(pcaDir, paste0(pcas_prefix, "channel_distributions.pdf")))

m <- match(panel_mass[cols], panel$Isotope)
ttl <- panel$Antigen[m]

colors <- as.numeric(as.factor(md$condition))

# for every channel
for(i in 1:length(cols)) {
  
  ii <- o[i]
  dens <- lapply(es, function(u) density( u[,ii] ))
  mx <- max(sapply(dens, function(u) max(u$y)))
  
  # for every sample
  k <- prs$marker==ttl[ii]
  for(s in 1:length(es)) {
    if(s==1)
      plot( dens[[s]], ylim=c(0,mx), xlab="", lwd=3, col=colors[s], 
        main=paste(panel_mass[cols][ii],"/",
          as.character(ttl[ii]),"/",
          round(prs$avg_score[k],2)))
    else
      lines( dens[[s]], lwd=3, col=colors[s])
    nm <- paste0(names(es)," (",round(prs[k,names(es)],2),")")
    legend("topright", nm, lwd=1, col = colors, cex = 0.7)
  }
  
}
dev.off()





### Plot chanel distributions with ggplot - not finished

# samp <- rep( names(fcsT), sapply(fcsT, nrow) )
# 
# e <- do.call("rbind",es)
# colnames(e) <- fcs_panel$Antigen[cols]
# 
# 
# ggdf <- data.frame(samp = samp, e)
# 
# 
# pdf(file.path(pcaDir, paste0(pcas_prefix, "channel_distributions.pdf")))
# for (i in seq(length(ggp)))
#   print(ggp[[i]])
# dev.off()





sessionInfo()













################################
### 01_pcascores done!
################################