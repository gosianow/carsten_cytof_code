##############################################################################
## <<carsten_cytof_pipeline.R>>

# BioC 3.3
# Updated 25 July 2016 

##############################################################################
Sys.time()
##############################################################################

# Load packages
library(flowCore)
library(gdata)
library(FlowSOM)
library(Repitools)
library(gplots)
library(Rtsne)
library(ggplot2)
library(plyr)
library(reshape2)
library(coop)
library(pheatmap)
library(RColorBrewer)

##############################################################################
# Test arguments
##############################################################################

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
# outdir_cd4='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01_CD4/010_cleanfcs'
# outdir_cd8='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01_CD8/010_cleanfcs'
# save_cd4_cd8=TRUE
# keep_cd4='CD4'
# keep_cd8='CD8'
# pca_score_cutoff=3
# tsne_pmin=1500

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

if(save_cd4_cd8){
  
  if( !file.exists(outdir_cd4) ) dir.create(outdir_cd4, recursive = TRUE)
  if( !file.exists(outdir_cd8) ) dir.create(outdir_cd8, recursive = TRUE)
  
}

# define directories
fcsDir <- "010_cleanfcs"; if( !file.exists(fcsDir) ) dir.create(fcsDir)
pcaDir <- "020_pcascores"; if( !file.exists(pcaDir) ) dir.create(pcaDir)
hmDir <- "030_heatmaps"; if( !file.exists(hmDir) ) dir.create(hmDir)
sneDir <- "040_tsnemaps"; if( !file.exists(sneDir) ) dir.create(sneDir)
frqDir <- "050_frequencies"; if( !file.exists(frqDir) ) dir.create(frqDir)
crsDir <- "060_crossplots"; if( !file.exists(crsDir) ) dir.create(crsDir)
dmpDir <- "070_dumpfcs"; if( !file.exists(dmpDir) ) dir.create(dmpDir)


### Set seed to be reproducible with FlowSOM clustering
options(myseed = 1234)
set.seed(1234)


### Colors for 20 clusters 
# ggplot palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=60 , c=100)[1:n]
}

ramp <- c(colorRampPalette(brewer.pal(12,"Paired"))(12)[-c(11)],  gg_color_hue(9) )



# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

md <- read.xls("metadata.xlsx",stringsAsFactors=FALSE)

# define FCS file names
f <- file.path(fcsDir,md$filename)
names(f) <- md$shortname

# associate group
grp <- md$condition

# read panel, pick which columns to use
panel <- read.xls("panel.xlsx",stringsAsFactors=FALSE)


# read raw FCS files in
fcs <- lapply(f, read.FCS)




### ugly way to get the mass of each column .. to match against panel
ss1 <- strsplit( colnames(fcs[[1]]), "Di" )
ss2 <- strsplit( sapply(ss1,.subset,1), "[a-zY]" )
panel_mass <- as.integer(sapply(ss2, tail, 1))

cols <- which( panel_mass %in% panel$Isotope[panel$transform==1] )

m <- match(panel_mass[cols], panel$Isotope)
panel$Antigen[m]

# arc-sin-h the columns specific 
fcsT <- lapply(fcs, function(u) {
  e <- exprs(u)
  e[,cols] <- asinh( e[,cols] / 5 )
  exprs(u) <- e
  u
})

# extract columns of interest to list of tables
es <- lapply(fcsT, function(u) {
  exprs(u)[,cols]
})







# -------------------------------------
# run Levine et al. 2015 marker scoring,
# plot marginal distributions
# -------------------------------------


doPRINCOMP <- function(z, ncomp=3) {
  cat(".")
  pr <- prcomp(z)  # default is to center but not scale
  score <- outer( rep(1,ncol(z)), pr$sdev[1:ncomp]^2 ) * abs(pr$rotation[,1:ncomp])
  rowSums(score)
}


prs <- sapply(es, doPRINCOMP)
rmprs <- rowMeans(prs)
o <- order(rmprs,decreasing=TRUE)
m <- match(panel_mass[cols], panel$Isotope)
prs <- data.frame(mass=rownames(prs), marker=panel$Antigen[m], round(prs,3), avg_score=rmprs)[o,]

write.table(prs, file=file.path(pcaDir,"princompscore_by_sample.xls"), sep="\t", row.names=FALSE, quote=FALSE)


### Plot chanel distributions
pdf(file.path(pcaDir,"channel_distributions.pdf"))

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





# -------------------------------------
# selected columns


scols <- which(colnames(fcs[[1]]) %in% rownames(prs)[prs$avg_score > pca_score_cutoff])

xcols <- setdiff(cols, scols)



# -------------------------------------
# run FlowSOM and make heatmaps
# -------------------------------------

# Number of clusters
nmetaclusts <- 20

set.seed(1234)
fs <- as(fcsT,"flowSet")
system.time( fsom <- FlowSOM::ReadInput(fs, transform = FALSE, scale = FALSE) )
system.time( fsom <- FlowSOM::BuildSOM(fsom, colsToUse = scols) )
system.time( fsom_mc <- FlowSOM::metaClustering_consensus(fsom$map$codes, k=nmetaclusts) )


save(fsom, file = file.path(hmDir, "fsom.rda"))
save(fsom_mc, file = file.path(hmDir, "fsom_mc.rda"))


# get cluster ids
clust <- fsom_mc[fsom$map$mapping[,1]]
(freq_clust <- table(clust))


### Save clustering results

freq_clust_out <- data.frame(cluster = names(freq_clust), freq = as.numeric(freq_clust))
write.table(freq_clust_out, file=file.path(hmDir,"cluster_counts.xls"), row.names=FALSE, quote=FALSE, sep="\t")

clust_out <- data.frame(cluster = clust, stringsAsFactors = FALSE)
write.table(clust_out, file = file.path(hmDir, "clustering.xls"), row.names=FALSE, quote=FALSE, sep="\t")




heatmap_wrapper <- function(clust, labels, fcsT, scols, xcols, suffix = "", pdf_width = 9, pdf_height = 6, xspace = 2){
  # clust - vector with clusters
  # labels - data frame with cluster and label columns
  
 
  ### Calculate cluster frequencies
  samp <- rep( names(fcsT), sapply(fcsT,nrow) )
  freq <- table( cluster = clust, samp )
  # Use labels as new cluster names
  mlab <- match(rownames(freq), labels$cluster)
  rownames(freq) <- labels$label[mlab]
  
  prop <- t(t(freq) / colSums(freq))
  # normalize 
  rp <- sweep(prop, 1, STATS=rowMeans(prop), FUN="-")
  

  ### Plot expression of diff markers in each cluster

  # re-extract the data as list of tables
  es <- lapply(fcsT, function(u) {
    exprs(u)[,scols]
  })
  
  esX <- lapply(fcsT, function(u) {
    exprs(u)[,xcols]
  })
  
  e <- do.call("rbind",es)
  eX <- do.call("rbind",esX)
  
  
  ss1 <- strsplit( colnames(e), "Di" )
  ss2 <- strsplit( sapply(ss1,.subset,1), "[a-zY]" )
  panel_mass <- as.integer(sapply(ss2, tail, 1))
  
  ss1X <- strsplit( colnames(eX), "Di" )
  ss2X <- strsplit( sapply(ss1X,.subset,1), "[a-zY]" )
  panel_massX <- as.integer(sapply(ss2X, tail, 1))
  
  
  m <- match(panel_mass, panel$Isotope)
  el <- e
  colnames(el) <- panel$Antigen[m]
  
  
  mlab <- match(clust, labels$cluster)
  
  df <- data.frame(el, clust = labels$label[mlab], samp)
  df2 <- melt(df,id=c("clust","samp"))
  
  ggp <- ggplot(df2, aes(x=value)) + geom_density(size=1,adjust=5) + 
    facet_grid(clust~variable,scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  pdf(file.path(hmDir,paste0("clust_distros",suffix,".pdf")),w=ncol(el),h=nmetaclusts)
  print(ggp)
  dev.off()
  
  
  mX <- match(panel_massX, panel$Isotope)
  elX <- eX
  colnames(elX) <- panel$Antigen[mX]
  
  a <- aggregate( el, by=list(clust), FUN=median)
  aX <- aggregate( elX, by=list(clust), FUN=median)
  
  
  ### Make heatmaps
  
  ### Prepare expression data
  
  # normalize to 0-1
  rng <- apply(el,2,quantile,p=c(.01,.99))
  for(i in 1:ncol(el)) {
    el[,i] <- (el[,i]-rng[1,i])/(rng[2,i]-rng[1,i])
  }
  el[el<0] <- 0
  el[el>1] <- 1
  
  
  # normalize to 0-1
  rng <- apply(elX,2,quantile,p=c(.01,.99))
  for(i in 1:ncol(elX)) {
    elX[,i] <- (elX[,i]-rng[1,i])/(rng[2,i]-rng[1,i])
  }
  elX[elX<0] <- 0
  elX[elX>1] <- 1
  
  a <- aggregate( el, by=list(clust), FUN=median)
  aX <- aggregate( elX, by=list(clust), FUN=median)
  
  mlab <- match(a$Group.1, labels$cluster)
  rownames(a) <- rownames(aX) <- labels$label[mlab]
  
  
  ### Prepare colors and breaks for the heatmaps
  cL <- NULL
  ncol <- 20
  
  palette <- colorRampPalette(c("green", "black", "red"),space = "Lab")
  br <- seq(-20,20,length=ncol)
  col <- palette(ncol-1)
  cL[[1]] <- list(breaks=br,colors=col)
  
  palette <- colorRampPalette(c("blue", "orange", "red"),space = "Lab")
  col <- palette(ncol-1)
  br <- seq(0,40,length=ncol)
  cL[[2]] <- list(breaks=br,colors=col)
  
  palette <- colorRampPalette(c("darkblue", "yellow", "red"),space = "Lab")
  col <- palette(ncol-1)
  br <- seq(0,1,length=ncol)
  cL[[4]] <- cL[[3]] <- list(breaks=br,colors=col)
  
  
  testD <- list(rp*100, prop*100, a[,-1], aX[,-1])
  
  testD[[1]][testD[[1]]>max(cL[[1]]$breaks)] <- max(cL[[1]]$breaks)
  testD[[1]][testD[[1]]<min(cL[[1]]$breaks)] <- min(cL[[1]]$breaks)
  testD[[2]][testD[[2]]>max(cL[[2]]$breaks)] <- max(cL[[2]]$breaks)
  testD[[2]][testD[[2]]<min(cL[[2]]$breaks)] <- min(cL[[2]]$breaks)
  
  
  ### Heatmaps will all the clusters
  pdf(file.path(hmDir,paste0("multiheatmap",suffix,".pdf")), width = pdf_width, height = pdf_height)
  multiHeatmap(testD,cL, xspace = xspace, ystarts=c(.25,.9,.925,.95,.98), clabelcex=.5)
  dev.off()
  
  
  ### Heatmaps with larger clusters
  w <- which(rowSums(prop>=.05)>=5)
  testDs <- lapply(testD, function(u) u[w,,drop=FALSE])
  
  pdf(file.path(hmDir,paste0("multiheatmap_subset",suffix,".pdf")), width = pdf_width, height = pdf_height)
  multiHeatmap(testDs,cL, xspace = xspace, ystarts=c(.25,.9,.925,.95,.98), clabelcex=.5)
  dev.off()
  
  return(testD)
  
}



labels <- data.frame(cluster = 1:nmetaclusts, label = factor(1:nmetaclusts))

testD <- heatmap_wrapper(clust = clust, labels = labels, fcsT = fcsT, scols = scols, xcols = xcols, suffix = "", pdf_width = 10, pdf_height = 6)



# -------------------------------------
### Cluster merging

if(file.exists("cluster_merging.xlsx")){
  
  cm <- read.xls("cluster_merging.xlsx")
  
  # remove spaces in labels bcs they are problematic...
  cm$label <- factor(cm$label, labels = gsub(" ", "_", levels(cm$label))) 
  
  # Get new merged clustering
  clustm <- factor(clust, levels = cm$old_cluster)
  levels(clustm) <- cm$new_cluster
  clustm <- as.numeric(as.character(clustm))
  
  labels <- unique(cm[,c("new_cluster","label")])
  colnames(labels) <- c("cluster", "label")
  
  labels <- labels[order(labels$cluster, decreasing = FALSE), ]
  labels$label <- factor(labels$label, levels = unique(labels$label))
  
  
  testDm <- heatmap_wrapper(clust = clustm, labels = labels, fcsT = fcsT, scols = scols, xcols = xcols, suffix = "_merged", pdf_width = 15, pdf_height = 3, xspace = 7)
  
  
  ### Save fcs files for clusters CD4 and CD8
  if(save_cd4_cd8){
    
    nrow_fcs <- sapply(fcs,nrow)
    mlab <- match(clustm, labels$cluster)
    
    clustmList <- split(labels$label[mlab], rep(names(nrow_fcs), nrow_fcs))
    
    writeOutCluster <- function(u,v,z,keep="CD4", outdir) {
      # u - cluster; v - flowFrame; z - original filename
      fn <- file.path(outdir, basename(z))
      cat(fn,"\n")
      # out <- v[u==keep, ]
      out <- v[grep(keep, u), ]
      write.FCS(out, fn)
    }
    
    m <- match(names(fcs), names(clustmList))
    clustmList <- clustmList[m]
    
    dummy <- mapply(writeOutCluster, u = clustmList, v = fcs, z = f, keep=keep_cd4, outdir = outdir_cd4)
    dummy <- mapply(writeOutCluster, u = clustmList, v = fcs, z = f, keep=keep_cd8, outdir = outdir_cd8)
    
    
  }
  
  
  
  
}


if(file.exists("cluster_merging2.xlsx")){
  
  cm <- read.xls("cluster_merging2.xlsx")
  
  # remove spaces in labels bcs they are problematic...
  cm$label <- factor(cm$label, labels = gsub(" ", "_", levels(cm$label))) 
  
  # Get new merged clustering
  clustm <- factor(clust, levels = cm$old_cluster)
  levels(clustm) <- cm$new_cluster
  clustm <- as.numeric(as.character(clustm))
  
  labels <- unique(cm[,c("new_cluster","label")])
  colnames(labels) <- c("cluster", "label")
  
  labels <- labels[order(labels$cluster, decreasing = FALSE), ]
  labels$label <- factor(labels$label, levels = unique(labels$label))
  labels
  
  testDm2 <- heatmap_wrapper(clust = clustm, labels = labels, fcsT = fcsT, scols = scols, xcols = xcols, suffix = "_merged2", pdf_width = 15, pdf_height = 3, xspace = 7)
  
  
}






# ---------------------------------------
# tSNE analyses
# ---------------------------------------


# re-extract the data as list of tables
es <- lapply(fcsT, function(u) {
  exprs(u)[,scols]
})


e <- do.call("rbind",es)

ss1 <- strsplit( colnames(e), "Di" )
ss2 <- strsplit( sapply(ss1,.subset,1), "[a-zY]" )
panel_mass <- as.integer(sapply(ss2, tail, 1))

m <- match(panel_mass, panel$Isotope)
el <- e
colnames(el) <- panel$Antigen[m]


# normalize to 0-1
rng <- apply(el,2,quantile,p=c(.01,.99))
for(i in 1:ncol(el)) {
  el[,i] <- (el[,i]-rng[1,i])/(rng[2,i]-rng[1,i])
}
el[el<0] <- 0
el[el>1] <- 1


# find duplicates
dups <- duplicated(el)  
w <- which(!dups)


### Data subsampling
samp <- rep( names(fcsT), sapply(fcsT,nrow) )

inds <- split(1:nrow(el), samp)  # create indices by sample

ts <- table(samp)
(ns <- pmin(ts, tsne_pmin))  # per-sample, how many cells to downsample

# get subsampled indices
subs <- mapply(function(u,v) {
  set.seed(1234)
  s <- sample(u,ns[v],replace = FALSE)
  intersect(s,w)
},inds,names(inds))

s <- unlist(subs)


### Run tSNE
set.seed(1234)
rtsne_out <- Rtsne(el[s,], pca = FALSE, verbose = TRUE)

save(rtsne_out, file = file.path(sneDir, "rtsne_out.rda"))


plot_tsne_wrapper <- function(rtsne_out, clust, group, el, prefix = "tSNE_all_", cmin, colors){
  
  
  clust <- factor(clust)
  group <- factor(group)
  
  ### Plot of tsne - all cells, all clusters
  cat("Plot 1\n")
  df <- data.frame(X = rtsne_out$Y[,1], Y = rtsne_out$Y[,2], cluster = clust, group = group)
  
  ggp <- ggplot(df,  aes(x = X, y = Y, color=cluster)) +
    geom_point(size=1) +
    facet_wrap(~ group) +
    labs(x = "tSNE 1", y="tSNE 2")+ 
    theme_bw() +
    theme(strip.text.x = element_text(size=15, face="bold"),
      axis.title.x  = element_text(size=15, face="bold"),
      axis.title.y  = element_text(size=15, face="bold")) +
    scale_color_manual(values = colors) +
    guides(colour = guide_legend(override.aes = list(size = 5)))
  
  pdf(file.path(sneDir,paste0(prefix, "clusters.pdf")), w=15, h=10)                 
  print(ggp)
  dev.off()
  
  
  
  ### Plot of tsne - cells that are close to the centers of clusters (no outliers), all clusters
  cat("Plot 2\n")
  a <- aggregate(el, by = list(clust), FUN=median)
  
  dists <- distse <- rep(NA, nrow(el))
  
  clust_lev <- levels(clust)
  
  for(i in 1:nlevels(clust)) {
    cent <- a[i,-1,drop=FALSE]
    data <- el[clust == clust_lev[i], , drop=FALSE]
    d <- 1 - cosine( t(rbind(cent,data)) )[1,-1]
    dists[clust == clust_lev[i]] <- d
    cent <- matrix(rep( as.numeric(cent), nrow(data)), byrow=TRUE, ncol=ncol(cent))
    #distse[clust[s]==i] <- sqrt(rowSums((cent-data)^2))
    distse[clust == clust_lev[i]] <- rowMeans(abs(cent-data))
  }
  
  df_sub <- df[distse<.25,]
  
  ggp <- ggplot(df_sub,  aes(x = X, y = Y, color = cluster)) +
    geom_point(size=1) + 
    facet_wrap( ~ group) +
    labs(x = "tSNE 1", y="tSNE 2")+ 
    theme_bw() +
    theme(strip.text.x = element_text(size=15, face="bold"),
      axis.title.x  = element_text(size=15, face="bold"),
      axis.title.y  = element_text(size=15, face="bold"))+
    scale_color_manual(values = colors) +
    guides(colour = guide_legend(override.aes = list(size = 5)))
  
  
  pdf(file.path(sneDir,paste0(prefix, "clusters_filtered.pdf")), w=15, h=10)                 
  print(ggp)
  dev.off()
  
  
  ### Plot of tsne - large clusters
  cat("Plot 3\n")
  tc <- table(df$cluster)
  kk <- names(tc[tc > cmin])
  
  df_sub <- df[df$cluster %in% kk,]
  
  ggp <- ggplot(df_sub,  aes(x = X, y = Y, color=cluster )) +
    geom_point(size=1) + 
    facet_wrap( ~ group)+ 
    labs(x = "tSNE 1", y="tSNE 2")+ 
    theme_bw() +
    theme(strip.text.x = element_text(size=15, face="bold"),
      axis.title.x  = element_text(size=15, face="bold"),
      axis.title.y  = element_text(size=15, face="bold")) +
    scale_color_manual(values = colors) + 
    guides(colour = guide_legend(override.aes = list(size = 5)))
  
  pdf(file.path(sneDir, paste0(prefix, "clusters_subset.pdf")), w=15,h=10)                 
  print(ggp)
  dev.off()
  
  return(":-)")
  
}



clust_sub <- sprintf("%02d", clust[s]) 
mm <- match(samp, md$shortname)
group <- md$condition[mm]
group_sub <- group[s]
el_sub <- el[s,]


plot_tsne_wrapper(rtsne_out, clust = clust_sub, group = group_sub, el = el_sub, prefix = "tSNE_all_", cmin = tsne_pmin*2/3, colors = ramp)




### Plot of tsne - merged clusters
if(file.exists("cluster_merging.xlsx")){
  
  cm <- read.xls("cluster_merging.xlsx")
  
  # remove spaces in labels bcs they are problematic...
  cm$label <- factor(cm$label, labels = gsub(" ", "_", levels(cm$label))) 
  
  ### Get new merged clustering
  clustm <- factor(clust, levels = cm$old_cluster)
  levels(clustm) <- cm$label
  
  # reorder the levels
  labels <- unique(cm[, c("label", "new_cluster")])
  labels <- labels[order(labels$new_cluster, decreasing = FALSE), ]
  
  clustm <- factor(clustm, levels = labels$label)
  
  clust_sub <- clustm[s]
  mm <- match(samp, md$shortname)
  group <- md$condition[mm]
  group_sub <- group[s]
  el_sub <- el[s,]
  
  plot_tsne_wrapper(rtsne_out, clust = clust_sub, group = group_sub, el = el_sub, prefix = "tSNE_merged_", cmin = tsne_pmin*2/3, colors = ramp[1:nrow(labels)])
  
}





### Plot of tsne - merged clusters
if(file.exists("cluster_merging2.xlsx")){
  
  cm <- read.xls("cluster_merging2.xlsx")
  
  # remove spaces in labels bcs they are problematic...
  cm$label <- factor(cm$label, labels = gsub(" ", "_", levels(cm$label))) 
  
  ### Get new merged clustering
  clustm <- factor(clust, levels = cm$old_cluster)
  levels(clustm) <- cm$label
  
  # reorder the levels
  labels <- unique(cm[, c("label", "new_cluster")])
  labels <- labels[order(labels$new_cluster, decreasing = FALSE), ]
  
  clustm <- factor(clustm, levels = labels$label)
  
  clust_sub <- clustm[s]
  mm <- match(samp, md$shortname)
  group <- md$condition[mm]
  group_sub <- group[s]
  el_sub <- el[s,]
  
  plot_tsne_wrapper(rtsne_out, clust = clust_sub, group = group_sub, el = el_sub, prefix = "tSNE_merged2_", cmin = tsne_pmin*2/3, colors = ramp[1:nrow(labels)])
  
}










# ---------------------------------------
# plots and analysis of cluster frequencies
# ---------------------------------------




plot_frequencies_wrapper <- function(clust, labels, fcsT, md, suffix = ""){
  # clust - vector with clusters
  # labels - data frame with cluster and label columns
  # md - metadata
  
  # calculate frequencies
  samp <- rep( names(fcsT), sapply(es,nrow) )
  freq <- table( cluster = clust, samp )
  
  # use labels as names of clusters
  mlab <- match(rownames(freq), labels$cluster)
  rownames(freq) <- labels$label[mlab]
  
  prop <- t(t(freq) / colSums(freq))
  
  ### Save the frequencies and proportions
  prop_out <- data.frame(cluster = rownames(prop), as.data.frame.matrix(prop * 100))
  freq_out <- data.frame(cluster = rownames(freq), as.data.frame.matrix(freq))

  write.table(prop_out, file=file.path(frqDir, paste0("frequencies", suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  write.table(freq_out, file=file.path(frqDir,paste0("counts",suffix,".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  
  
  ### Data frame preparation for ggplot
  
  ggdf <- melt(prop, value.name = "prop")
  
  # Add group info
  mm <- match(ggdf$samp, md$shortname)
  ggdf$group <- md$condition[mm]
  
  ggds <- ddply(ggdf, .(group, cluster), summarise, mean = mean(prop), sd = sd(prop))
  
  ggp <- list()
  
  clusters <- levels(ggdf$cluster)
  
  for(i in 1:nlevels(ggdf$cluster)){
    
    df <- ggdf[ggdf$cluster == clusters[i], , drop = FALSE]
    ds <- ggds[ggds$cluster == clusters[i], , drop = FALSE]
    
    ggp[[i]] <- ggplot(df, aes(x = group, y = prop)) + 
      geom_point(size=3, shape = 1, stroke = 1) + 
      geom_point(data=ds, aes(x=group, y=mean), colour='blue', size=2) +
      geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean-sd, ymax=mean+sd), colour='blue', width=0.05) +
      ggtitle(clusters[i]) + 
      theme_bw() +
      ylab("proportion") + 
      xlab("") +
      theme(axis.text.x  = element_text(size=15, face="bold"),
        axis.title.y  = element_text(size=15, face="bold"))

  }
  
  
  pdf(file.path(frqDir, paste0("cluster_frequencies", suffix, ".pdf")), w=7, h=5, onefile=TRUE)
  for (i in seq(length(ggp)))
    print(ggp[[i]])
  dev.off()
  
  return("Done!")
  
}



labels <- data.frame(cluster = 1:nmetaclusts, label = factor(1:nmetaclusts, labels = paste0("Cluster ", sprintf("%02d", 1:nmetaclusts))))


plot_frequencies_wrapper(clust = clust, labels = labels, fcsT = fcsT, md = md, suffix = "")



if(file.exists("cluster_merging.xlsx")){
  
  cm <- read.xls("cluster_merging.xlsx")
  
  # remove spaces in labels bcs they are problematic...
  cm$label <- factor(cm$label, labels = gsub(" ", "_", levels(cm$label))) 
  
  # Get new merged clustering
  clustm <- factor(clust, levels = cm$old_cluster)
  levels(clustm) <- cm$new_cluster
  clustm <- as.numeric(as.character(clustm))
  
  # Create labels data frame
  labels <- unique(cm[,c("new_cluster","label")])
  colnames(labels) <- c("cluster", "label")
  
  labels <- labels[order(labels$cluster, decreasing = FALSE), ]
  labels$label <- factor(labels$label, levels = unique(labels$label))
  
  
  plot_frequencies_wrapper(clust = clustm, labels = labels, fcsT = fcsT, md = md, suffix = "_merged")
  
  
  
}



if(file.exists("cluster_merging2.xlsx")){
  
  cm <- read.xls("cluster_merging2.xlsx")
  
  # remove spaces in labels bcs they are problematic...
  cm$label <- factor(cm$label, labels = gsub(" ", "_", levels(cm$label))) 
  
  # Get new merged clustering
  clustm <- factor(clust, levels = cm$old_cluster)
  levels(clustm) <- cm$new_cluster
  clustm <- as.numeric(as.character(clustm))
  
  # Create labels data frame
  labels <- unique(cm[,c("new_cluster","label")])
  colnames(labels) <- c("cluster", "label")
  
  labels <- labels[order(labels$cluster, decreasing = FALSE), ]
  labels$label <- factor(labels$label, levels = unique(labels$label))
  labels
  
  plot_frequencies_wrapper(clust = clustm, labels = labels, fcsT = fcsT, md = md, suffix = "_merged2")
  

}





sessionInfo()
