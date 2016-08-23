##############################################################################
## <<04_expression.R>>

# BioC 3.3
# Created 22 Aug 2016

##############################################################################
Sys.time()
##############################################################################

# Load packages
library(flowCore)
library(gdata)
library(FlowSOM)
library(Repitools)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2)
library(RColorBrewer)


##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_03'
path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_03.xlsx'
expr_prefix='23_03_pca1_merging3_'
path_clustering_observables='23_03_pca1_clustering_observables.xls'
path_clustering='23_03_pca1_merging3_clustering.xls'
path_clustering_labels='23_03_pca1_merging3_clustering_labels.xls'


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

prefix <- expr_prefix


# ------------------------------------------------------------
# define directories
# ------------------------------------------------------------

fcsDir <- "010_cleanfcs"; if( !file.exists(fcsDir) ) dir.create(fcsDir)
hmDir <- "030_heatmaps"; if( !file.exists(hmDir) ) dir.create(hmDir)
exprDir <- "080_expression"; if( !file.exists(exprDir) ) dir.create(exprDir)

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

samp <- rep( names(fcs), sapply(fcs, nrow) )

# ------------------------------------------------------------
# Load more data
# ------------------------------------------------------------

# clustering
clust <- read.table(file.path(hmDir, path_clustering), header = TRUE, sep = "\t", as.is = TRUE)
clust <- clust[, 1]

# clustering labels
labels <- read.table(file.path(hmDir, path_clustering_labels), header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))
rownames(labels) <- labels$cluster


# clustering_observables
if(!grepl("/", path_clustering_observables)){
  clustering_observables <- read.table(file.path(hmDir, path_clustering_observables), header = TRUE, sep = "\t", as.is = TRUE)
}else{
  clustering_observables <- read.table(path_clustering_observables, header = TRUE, sep = "\t", as.is = TRUE)
}

rownames(clustering_observables) <- clustering_observables$mass

clust_observ <- clustering_observables[clustering_observables$clustering_observable, "mass"]



# -------------------------------------
# get the isotope and antigen for fcs markers

m <- match(fcs_colnames, clustering_observables$mass)

fcs_panel <- data.frame(colnames = fcs_colnames, Isotope = clustering_observables$mass[m], Antigen = clustering_observables$marker[m], stringsAsFactors = FALSE)


# -------------------------------------

cols <- which(fcs_colnames %in% clustering_observables$mass)

# arc-sin-h the columns specific 
fcsT <- lapply(fcs, function(u) {
  e <- exprs(u)
  e[ , cols] <- asinh( e[ , cols] / 5 )
  exprs(u) <- e
  u
})



# ------------------------------------------------------------

### Indeces of observables used for clustering 

scols <- which(fcs_colnames %in% clust_observ)

# ordered by decreasing pca score
scols <- scols[order(clustering_observables[fcs_colnames[scols], "avg_score"], decreasing = TRUE)]


### Indeces of other observables

xcols <- setdiff(cols, scols)

# ordered by decreasing pca score
xcols <- xcols[order(clustering_observables[fcs_colnames[xcols], "avg_score"], decreasing = TRUE)]




# ------------------------------------------------------------
# Get expression
# ------------------------------------------------------------

### Get marker expression

es <- lapply(fcsT, function(u) {
  exprs(u)[, cols]
})

e <- do.call("rbind",es)
colnames(e) <- fcs_panel$Antigen[cols]

a <- aggregate( e, by=list(clust, samp), FUN=median)

mlab <- match(a$Group.1, labels$cluster)
a$label <- labels$label[mlab]

colnames(a)[1:2] <- c("cluster", "sample")

# normalize to 0-1
el <- e
rng <- apply(el,2,quantile,p=c(.01,.99))
for(i in 1:ncol(el)) {
  el[,i] <- (el[,i]-rng[1,i])/(rng[2,i]-rng[1,i])
}
el[el<0] <- 0
el[el>1] <- 1

al <- aggregate( el, by=list(clust, samp), FUN=median)

mlab <- match(al$Group.1, labels$cluster)
al$label <- labels$label[mlab]

colnames(al)[1:2] <- c("cluster", "sample")


a <- a[, c("cluster", "label", "sample", fcs_panel$Antigen[c(scols, xcols)])]
al <- al[, c("cluster", "label", "sample", fcs_panel$Antigen[c(scols, xcols)])]

### Save the median expression per cluster and sample

# raw
write.table(a, file.path(exprDir, paste0(prefix, "expr_cluster_sample_raw.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# normalized 
write.table(al, file.path(exprDir, paste0(prefix, "expr_cluster_sample_norm.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



# ---------------------------------------
### Plot expression
# ---------------------------------------

expr <- a

ggdf <- melt(expr, id.vars = c("cluster", "label", "sample"), value.name = "expr", variable.name = "marker")

## use labels as clusters
ggdf$cluster <- factor(ggdf$cluster, levels = labels$cluster, labels = labels$label)

## add group info
mm <- match(ggdf$samp, md$shortname)
ggdf$group <- factor(md$condition[mm])

## replace _ with \n
levels(ggdf$group) <- gsub("_", "\n", levels(ggdf$group))

## calculate mean and sd for the error bars on the plot
ggds <- ddply(ggdf, .(group, cluster, marker), summarise, mean = mean(expr), sd = sd(expr))


## plot each cluster as a separate page in the pdf file
ggp <- list()

markers <- levels(ggdf$marker)

for(i in 1:nlevels(ggdf$marker)){
  # i = 1
  
  df <- ggdf[ggdf$marker == markers[i], , drop = FALSE]
  ds <- ggds[ggds$marker == markers[i], , drop = FALSE]
  
  ggp[[i]] <- ggplot(df, aes(x = group, y = expr)) +
    geom_jitter(size=2, shape = 16, width = 0.5, height = 0) +
    geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean, ymax=mean), colour='black', width=0.4) +
    geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean-sd, ymax=mean+sd), colour='black', width=0.25) +
    facet_wrap(~ cluster, scales = "fixed") +
    ggtitle(markers[i]) +
    theme_bw() +
    ylab("Expression") +
    xlab("") +
    theme(axis.text.x = element_text(size=12, face="bold"), 
      axis.title.y = element_text(size=12, face="bold"))
  
}


pdf(file.path(exprDir, paste0(prefix, "expr_cluster_sample_raw.pdf")), w=12, h=10, onefile=TRUE)
for(i in seq(length(ggp)))
  print(ggp[[i]])
dev.off()









sessionInfo()













################################
### 04_expression done!
################################