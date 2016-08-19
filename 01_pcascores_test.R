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
library(RColorBrewer)

##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_03'
pcas_prefix='update_'
path_panel='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_panels/panel3.xlsx'
path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_03.xlsx'

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

fcsDir <- "010_cleanfcs"; if( !file.exists(fcsDir) ) dir.create(fcsDir)
pcaDir <- "020_pcascores"; if( !file.exists(pcaDir) ) dir.create(pcaDir)


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
panel <- read.xls(path_panel, stringsAsFactors=FALSE)

fcs_colnames <- colnames(fcs[[1]])

# get isotope mass of columns in fcs files.. to match against the panel
panel_mass <- as.numeric(gsub("[[:alpha:]]", "", fcs_colnames))

mm <- match(panel$Isotope, panel_mass)

panel$fcs_colname <- fcs_colnames[mm]


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


### Create sample info
samp <- rep(names(fcs), sapply(fcs, nrow))


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


# extract columns of interest to list of tables
es <- lapply(fcsT, function(u) {
  exprs(u)[,cols]
})


prs <- sapply(es, doPRINCOMP)
rmprs <- rowMeans(prs)


prs <- data.frame(mass = rownames(prs), marker = fcs_panel$Antigen[cols], avg_score=rmprs, round(prs,3))


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



# --------------------------------------------------------------------------
# My exploration of PCA
# --------------------------------------------------------------------------




z = es[[1]]; ncomp=3

## Score by Levine
pr <- prcomp(z, center = TRUE, scale. = FALSE)  # default is to center but not scale

score <- rowSums (outer( rep(1, ncol(z)), pr$sdev[1:ncomp]^2 ) * abs(pr$rotation[,1:ncomp]) )




## My testing
pr <- prcomp(z, center = TRUE, scale. = FALSE)  # default is to center but not scale
colSums(pr$rotation^2) # are all 1 because eigenvectors have length 1

sum(pr$sdev^2) # is equal to the number of variables when scale. = TRUE because all are scaled to have variance of 1


rowSums( (outer( rep(1, ncol(z)), pr$sdev) * pr$rotation)^2 ) # all 1 when scale. = TRUE (portion of the variables' variance being explained by components)


rowSums( (outer( rep(1, ncol(z)), pr$sdev[1:ncomp]) * pr$rotation[,1:ncomp])^2 )


rowSums( outer( rep(1, ncol(z)), pr$sdev^2) * pr$rotation^2 )


rowSums( outer( rep(1, ncol(z)), pr$sdev^2) * abs(pr$rotation) )

rowSums (outer( rep(1, ncol(z)), pr$sdev[1:ncomp]^2 ) * abs(pr$rotation[,1:ncomp]) )











sessionInfo()













################################
### 01_pcascores done!
################################