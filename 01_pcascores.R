##############################################################################
## <<01_pcascores.R>>

# BioC 3.3
# Created 27 July 2016
# Updated 31 Aug 2016

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

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-29_03all2_myeloid_merging3'
pcas_prefix='29mye_03_'
pcas_outdir='020_pcascores'
path_data='010_data/29mye_03_expr_raw.rds'
path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_29_03all2.xlsx'
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

prefix <- pcas_prefix

outdir <- pcas_outdir

if( !file.exists(outdir) ) 
  dir.create(outdir, recursive = TRUE)

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


# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

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

min_samp <- nrow(md)
  
  

doPRINCOMP <- function(z, ncomp=3) {
  # z = es[[1]]; ncomp=3
  
  if(nrow(z) < min_samp)
    return(rep(NA, ncol(z)))
  
  ## Score by Levine
  pr <- prcomp(z, center = TRUE, scale. = FALSE)  # default is to center but not scale

  score <- rowSums (outer( rep(1, ncol(z)), pr$sdev[1:ncomp]^2 ) * abs(pr$rotation[,1:ncomp]) )

  return(score)

}


# split expression per sample
es <- split(e, samp)

prs <- sapply(es, doPRINCOMP)
rmprs <- rowMeans(prs, na.rm = TRUE)

prs <- data.frame(mass = fcs_panel$fcs_colname, marker = fcs_panel$Antigen, avg_score = rmprs, round(prs, 4))

### Save ordered PCA scores
o <- order(rmprs, decreasing=TRUE)
prs <- prs[o,]

write.table(prs, file = file.path(outdir, paste0(prefix, "princompscore_by_sample.xls")), sep="\t", row.names=FALSE, quote=FALSE)


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


pdf(file.path(outdir, paste0(prefix, "princompscore_average.pdf")), width = 10, height = 7)
print(ggp1)
print(ggp2)
dev.off()


### Plot chanel distributions

pdf(file.path(outdir, paste0(prefix, "channel_distributions.pdf")))

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





# --------------------------------------------------------------------------
# My exploration of PCA
# --------------------------------------------------------------------------

# # split expression per sample
# es <- split(e, samp)
# 
# z = es[[1]]; ncomp=3
# 
# 
# ## Score by Levine
# pr <- prcomp(z, center = TRUE, scale. = FALSE)  # default is to center but not scale
# 
# rowSums (outer( rep(1, ncol(z)), pr$sdev[1:ncomp]^2 ) * abs(pr$rotation[,1:ncomp]) )
# 
# # ----------------------------
# ## My testing - scale = TRUE
# 
# pr <- prcomp(z, center = TRUE, scale. = TRUE)  # default is to center but not scale
# 
# colSums(pr$rotation^2) # are all 1 because eigenvectors have length 1
# 
# sum(pr$sdev^2) # is equal to the number of variables when scale. = TRUE because all are scaled to have variance of 1
# 
# 
# st_all <- rowSums( (outer( rep(1, ncol(z)), pr$sdev) * pr$rotation)^2 ) # all 1 when scale. = TRUE (portion of the variables' variance being explained by components)
# st_all
# st_comp <- rowSums( (outer( rep(1, ncol(z)), pr$sdev[1:ncomp]) * pr$rotation[,1:ncomp])^2 ) # (portion of the variables' variance being explained by components)
# st_comp
# 
# ### using squares
# 
# rowSums( outer( rep(1, ncol(z)), pr$sdev^2) * pr$rotation^2 )
# 
# score_ss_scaleT <- rowSums( outer( rep(1, ncol(z)), pr$sdev[1:ncomp]^2) * pr$rotation[,1:ncomp]^2 )
# 
# 
# ### using abs
# 
# rowSums( outer( rep(1, ncol(z)), pr$sdev^2) * abs(pr$rotation) )
# 
# score_levine_scaleT <- rowSums (outer( rep(1, ncol(z)), pr$sdev[1:ncomp]^2 ) * abs(pr$rotation[,1:ncomp]) )
# 
# 
# # ----------------------------
# ## My testing - scale = FALSE
# 
# pr <- prcomp(z, center = TRUE, scale. = FALSE)  # default is to center but not scale
# 
# colSums(pr$rotation^2) # are all 1 because eigenvectors have length 1
# rowSums(pr$rotation^2) # are also all 1 (but I do not know why...)
# 
# sum(pr$sdev^2) # is equal to the number of variables when scale. = TRUE because all are scaled to have variance of 1
# 
# 
# sf_all <- rowSums( (outer( rep(1, ncol(z)), pr$sdev) * pr$rotation)^2 ) # all 1 when scale. = TRUE (portion of the variables' variance being explained by components)
# sf_all
# 
# sf_comp <- rowSums( (outer( rep(1, ncol(z)), pr$sdev[1:ncomp]) * pr$rotation[,1:ncomp])^2 ) # (portion of the variables' variance being explained by components)
# sf_comp
# 
# 
# sf_comp/sf_all
# 
# 
# ### using squares
# 
# rowSums( outer( rep(1, ncol(z)), pr$sdev^2) * pr$rotation^2 )
# 
# score_ss_scaleF <- rowSums( outer( rep(1, ncol(z)), pr$sdev[1:ncomp]^2) * pr$rotation[,1:ncomp]^2 )
# 
# 
# ### using abs
# 
# rowSums( outer( rep(1, ncol(z)), pr$sdev^2) * abs(pr$rotation) )
# 
# score_levine_scaleF <- rowSums (outer( rep(1, ncol(z)), pr$sdev[1:ncomp]^2 ) * abs(pr$rotation[,1:ncomp]) )


# --------------------------------------------------------------------------
# My exploration of PCA - run for all the samples
# --------------------------------------------------------------------------

es <- split(e, samp)

# ------------------------------
### Plot the fraction of variance that is explained by PCs

pc_frac_expl <- sapply(es, function(z){
  # z = es[[1]]
  
  if(nrow(z) < min_samp)
    return(rep(NA, ncol(z)))
  
  pr <- prcomp(z, center = TRUE, scale. = FALSE)
  
  cumsum(pr$sdev^2) / sum(pr$sdev^2, na.rm = TRUE)
  
})


pc_frac_expl <- data.frame(NoPC = 1:nrow(pc_frac_expl), pc_frac_expl)

ggdf <- melt(pc_frac_expl, id.vars = "NoPC", variable.name = "sample", value.name = "fraction")


ggp <- ggplot(ggdf, aes(x = NoPC, y = fraction, group = sample)) +
  geom_line() +
  geom_vline(xintercept = 3, color = "green") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  xlab("Number of PCs") +
  ylab("Fraction of variance explained") 


pdf(file.path(outdir, paste0(prefix, "expl_frac_var_explained.pdf")), width = 10, height = 7)
print(ggp)
dev.off()




# ------------------------------
### Plot the Levine PCA score and the max that it can get

ncomp=3

score_levine <- sapply(es, function(z){
  
  if(nrow(z) < min_samp)
    return(rep(NA, ncol(z)))
  
  pr <- prcomp(z, center = TRUE, scale. = FALSE)
  rowSums (outer( rep(1, ncol(z)), pr$sdev[1:ncomp]^2 ) * abs(pr$rotation[,1:ncomp]) )
  
})

score_levine <- data.frame(fcs_panel[, c("fcs_colname", "Antigen")], score_levine)

ggdf <- melt(score_levine, id.vars = c("fcs_colname", "Antigen"), variable.name = "sample", value.name = "score")

ggp <- ggplot(ggdf, aes(x = Antigen, y = score)) +
  geom_jitter(width = 0.5, height = 0) +
  stat_summary(fun.y = "mean", colour = "orange", size = 4, geom = "point", shape = 18) +
  theme_bw() +
  xlab("Marker") +
  ylab("PCA score") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12))


pdf(file.path(outdir, paste0(prefix, "expl_score_levine.pdf")), width = 10, height = 7)
print(ggp)
dev.off()


### max
score_levine_max <- sapply(es, function(z){
  
  if(nrow(z) < min_samp)
    return(rep(NA, ncol(z)))
  
  pr <- prcomp(z, center = TRUE, scale. = FALSE)
  rowSums (outer( rep(1, ncol(z)), pr$sdev^2 ) * abs(pr$rotation) )
  
})

score_levine_max <- data.frame(fcs_panel[, c("fcs_colname", "Antigen")], score_levine_max)

ggdf <- melt(score_levine_max, id.vars = c("fcs_colname", "Antigen"), variable.name = "sample", value.name = "score")

ggp <- ggplot(ggdf, aes(x = Antigen, y = score)) +
  geom_jitter(width = 0.5, height = 0) +
  stat_summary(fun.y = "mean", colour = "orange", size = 4, geom = "point", shape = 18) +
  theme_bw() +
  xlab("Marker") +
  ylab("PCA score") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12))


pdf(file.path(outdir, paste0(prefix, "expl_score_levine_max.pdf")), width = 10, height = 7)
print(ggp)
dev.off()


### ratio

score_levine_ratio <- data.frame(fcs_panel[, c("fcs_colname", "Antigen")], score_levine[, md$shortname] / score_levine_max[, md$shortname] * 100)

ggdf <- melt(score_levine_ratio, id.vars = c("fcs_colname", "Antigen"), variable.name = "sample", value.name = "score")

ggp <- ggplot(ggdf, aes(x = Antigen, y = score)) +
  geom_jitter(width = 0.5, height = 0) +
  stat_summary(fun.y = "mean", colour = "orange", size = 4, geom = "point", shape = 18) +
  coord_cartesian(ylim = c(0, 100)) +
  theme_bw() +
  xlab("Marker") +
  ylab("PCA score ratio") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12))


pdf(file.path(outdir, paste0(prefix, "expl_score_levine_ratio.pdf")), width = 10, height = 7)
print(ggp)
dev.off()


# ------------------------------
### Plot the SS PCA score and the max that it can get

ncomp=3

score_ss <- sapply(es, function(z){
  
  if(nrow(z) < min_samp)
    return(rep(NA, ncol(z)))
  
  pr <- prcomp(z, center = TRUE, scale. = FALSE)
  rowSums (outer( rep(1, ncol(z)), pr$sdev[1:ncomp]^2 ) * pr$rotation[,1:ncomp]^2 )
  
})

score_ss <- data.frame(fcs_panel[, c("fcs_colname", "Antigen")], score_ss)

ggdf <- melt(score_ss, id.vars = c("fcs_colname", "Antigen"), variable.name = "sample", value.name = "score")


ggp <- ggplot(ggdf, aes(x = Antigen, y = score)) +
  geom_jitter(width = 0.5, height = 0) +
  stat_summary(fun.y = "mean", colour = "orange", size = 4, geom = "point", shape = 18) +
  theme_bw() +
  xlab("Marker") +
  ylab("PCA score") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12))


pdf(file.path(outdir, paste0(prefix, "expl_score_ss_scalef.pdf")), width = 10, height = 7)
print(ggp)
dev.off()


### max
score_ss_max <- sapply(es, function(z){
  
  if(nrow(z) < min_samp)
    return(rep(NA, ncol(z)))
  
  pr <- prcomp(z, center = TRUE, scale. = FALSE)
  rowSums (outer( rep(1, ncol(z)), pr$sdev^2 ) * pr$rotation^2 )
  
})

score_ss_max <- data.frame(fcs_panel[, c("fcs_colname", "Antigen")], score_ss_max)

ggdf <- melt(score_ss_max, id.vars = c("fcs_colname", "Antigen"), variable.name = "sample", value.name = "score")


ggp <- ggplot(ggdf, aes(x = Antigen, y = score)) +
  geom_jitter(width = 0.5, height = 0) +
  stat_summary(fun.y = "mean", colour = "orange", size = 4, geom = "point", shape = 18) +
  theme_bw() +
  xlab("Marker") +
  ylab("PCA score") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12))


pdf(file.path(outdir, paste0(prefix, "expl_score_ss_scalef_max.pdf")), width = 10, height = 7)
print(ggp)
dev.off()


### ratio
score_ss_ratio <- data.frame(fcs_panel[, c("fcs_colname", "Antigen")], score_ss[, md$shortname] / score_ss_max[, md$shortname] * 100)

ggdf <- melt(score_ss_ratio, id.vars = c("fcs_colname", "Antigen"), variable.name = "sample", value.name = "score")

ggp <- ggplot(ggdf, aes(x = Antigen, y = score)) +
  geom_jitter(width = 0.5, height = 0) +
  stat_summary(fun.y = "mean", colour = "orange", size = 4, geom = "point", shape = 18) +
  coord_cartesian(ylim = c(0, 100)) +
  theme_bw() +
  xlab("Marker") +
  ylab("PCA score ratio") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12))


pdf(file.path(outdir, paste0(prefix, "expl_score_ss_scalef_ratio.pdf")), width = 10, height = 7)
print(ggp)
dev.off()


## scale = TRUE
score_ss_scalet <- sapply(es, function(z){
  
  if(nrow(z) < min_samp)
    return(rep(NA, ncol(z)))
  
  pr <- prcomp(z, center = TRUE, scale. = TRUE)
  rowSums (outer( rep(1, ncol(z)), pr$sdev[1:ncomp]^2 ) * pr$rotation[,1:ncomp]^2 )
  
})

score_ss_scalet <- data.frame(fcs_panel[, c("fcs_colname", "Antigen")], score_ss_scalet)

ggdf <- melt(score_ss_scalet, id.vars = c("fcs_colname", "Antigen"), variable.name = "sample", value.name = "score")


ggp <- ggplot(ggdf, aes(x = Antigen, y = score)) +
  geom_jitter(width = 0.5, height = 0) +
  stat_summary(fun.y = "mean", colour = "orange", size = 4, geom = "point", shape = 18) +
  theme_bw() +
  xlab("Marker") +
  ylab("PCA score") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12))


pdf(file.path(outdir, paste0(prefix, "expl_score_ss_scalet.pdf")), width = 10, height = 7)
print(ggp)
dev.off()


# ------------------------------
### Plot correlation between the score and the score ratio 

## score ss
ggdfs <- melt(score_ss, id.vars = c("fcs_colname", "Antigen"), variable.name = "sample", value.name = "score")
ggdfsr <- melt(score_ss_ratio, id.vars = c("fcs_colname", "Antigen"), variable.name = "sample", value.name = "score_ratio")

ggdf <- merge(ggdfs, ggdfsr, by = c("fcs_colname", "Antigen", "sample"))

ggp <- ggplot(ggdf, aes(x = score_ratio, y = score)) +
  geom_point() +
  coord_cartesian(xlim = c(0, 100)) +
  theme_bw() +
  xlab("Score ratio") +
  ylab("Score")


pdf(file.path(outdir, paste0(prefix, "expl_score_ss_scalef_correlation.pdf")), width = 7, height = 7)
print(ggp)
dev.off()



## score Levine
ggdfs <- melt(score_levine, id.vars = c("fcs_colname", "Antigen"), variable.name = "sample", value.name = "score")
ggdfsr <- melt(score_levine_ratio, id.vars = c("fcs_colname", "Antigen"), variable.name = "sample", value.name = "score_ratio")

ggdf <- merge(ggdfs, ggdfsr, by = c("fcs_colname", "Antigen", "sample"))

ggp <- ggplot(ggdf, aes(x = score_ratio, y = score)) +
  geom_point() +
  coord_cartesian(xlim = c(0, 100)) +
  theme_bw() +
  xlab("Score ratio") +
  ylab("Score")


pdf(file.path(outdir, paste0(prefix, "expl_score_levine_correlation.pdf")), width = 7, height = 7)
print(ggp)
dev.off()












sessionInfo()













################################
### 01_pcascores done!
################################