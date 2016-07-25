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

##############################################################################
# Test arguments
##############################################################################

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
# outdir_cd4='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01_CD4/010_cleanfcs'
# outdir_cd8='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01_CD8/010_cleanfcs'
# save_cd4_cd8=TRUE
# keep_cd4='CD4'
# keep_cd8='CD8'

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


if( !file.exists(outdir_cd4) ) dir.create(outdir_cd4, recursive = TRUE)
if( !file.exists(outdir_cd8) ) dir.create(outdir_cd8, recursive = TRUE)

# define directories
fcsDir <- "010_cleanfcs"; if( !file.exists(fcsDir) ) dir.create(fcsDir)
pcaDir <- "020_pcascores"; if( !file.exists(pcaDir) ) dir.create(pcaDir)
hmDir <- "030_heatmaps"; if( !file.exists(hmDir) ) dir.create(hmDir)
sneDir <- "040_tsnemaps"; if( !file.exists(sneDir) ) dir.create(sneDir)
frqDir <- "050_frequencies"; if( !file.exists(frqDir) ) dir.create(frqDir)
crsDir <- "060_crossplots"; if( !file.exists(crsDir) ) dir.create(crsDir)
dmpDir <- "070_dumpfcs"; if( !file.exists(dmpDir) ) dir.create(dmpDir)


# settings for FlowSOM
pca_score_cutoff <- 1.5
nmetaclusts <- 20
set.seed(1234)


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
# ss1 <- strsplit( colnames(fcs[[1]]), "Di" )
# ss2 <- strsplit( sapply(ss1,.subset,1), "[a-zY]" )
# panel_mass <- as.integer(sapply(ss2, tail, 1))

### way to get the mass of each column .. to match against panel
panel_mass <- as.numeric(gsub("[[:alpha:]]", "", colnames(fcs[[1]])))

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

write.table(prs, file=file.path(pcaDir,"princompscore_by_sample.xls"), 
  sep="\t", row.names=FALSE, quote=FALSE)



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
    legend("topright",nm,lwd=3,col=colors)
  }
  
}
dev.off()

# limma::plotMDS(prs[,-c(1,2,18)])

# -------------------------------------
# selected columns
scols <- which(colnames(fcs[[1]]) %in% rownames(prs)[prs$avg_score > pca_score_cutoff])

xcols <- setdiff(cols,scols)

# re-extract the data as list of tables
es <- lapply(fcsT, function(u) {
  exprs(u)[,scols]
})

esX <- lapply(fcsT, function(u) {
  exprs(u)[,xcols]
})




# -------------------------------------
# run FlowSOM and make heatmaps
# -------------------------------------
fs <- as(fcsT,"flowSet")
system.time( fsom <- FlowSOM::ReadInput(fs, transform = FALSE, scale = FALSE) )
system.time( fsom <- FlowSOM::BuildSOM(fsom, colsToUse = scols) )
system.time( mc_fsom <- FlowSOM::metaClustering_consensus(fsom$map$codes, k=nmetaclusts) )


# get cluster ids
clust <- mc_fsom[fsom$map$mapping[,1]]
table(clust)


# calculate frequencies
samp <- rep( names(f), sapply(es,nrow) )
freq <- table( cluster=clust, samp )
prop <- t(t(freq) / colSums(freq))


round( prop*100, 2)
k <- rowSums( prop > .005 ) >= 3
table(k)


# plot expression of diff markers in each cluster


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

df <- data.frame(el,clust,samp)
df2 <- melt(df,id=c("clust","samp"))

ggp <- ggplot(df2, aes(x=value)) + geom_density(size=1,adjust=5) + 
  facet_grid(clust~variable,scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf(file.path(hmDir,"clust_distros.pdf"),w=ncol(el),h=nmetaclusts)
print(ggp)
dev.off()


mX <- match(panel_massX, panel$Isotope)
elX <- eX
colnames(elX) <- panel$Antigen[mX]

a <- aggregate( el, by=list(clust), FUN=median)
aX <- aggregate( elX, by=list(clust), FUN=median)

#dists <- rep(NA,nrow(el))
#for(i in 1:nmetaclusts) {
#  cat(i,"\n")
#  cent <- a[i,-1,drop=FALSE]
#  data <- el[clust==i,,drop=FALSE]
#  d <- 1 - cosine( t(rbind(cent,data)) )[1,-1]
#  dists[clust==i] <- d
#}

#dfdist <- data.frame(clust,dists,stringsAsFactors=FALSE)

#pdf(file.path(hmDir,"cluster_cosinedistances.pdf"),width=9,h=6)
#p <- ggplot(dfdist, aes(dists, colour = factor(clust))) +
#       geom_density(size=2,alpha = 0.1)
#print(p)
#print(p + xlim(c(0,.4)))
#dev.off()

#ggplot(dfdist, aes(x=rm, y=dists, colour = factor(clust))) + 
#  geom_density_2d() + geom_point() + facet_wrap(~clust)

#ggplot(dfdist, aes(x=rm, y=dists, colour = factor(clust))) +
#  geom_point(size=.5) + geom_density_2d(colour="black") + facet_wrap(~clust) + ylim(c(0,.5))


#ggplot(dfdist, aes(x=rm, y=dists, colour = factor(clust))) + 
#  stat_density_2d(geom = "raster", aes(fill = ..density..^.25), contour = FALSE) + facet_wrap(~clust)

#p + scale_y_log10() + ylim(c(.01,50))


# make heatmaps

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

#dists <- rep(NA,nrow(el))
#for(i in 1:nmetaclusts) {
#  cat(i,"\n")
#  cent <- a[i,-1,drop=FALSE]
#  data <- el[clust==i,,drop=FALSE]
#  d <- 1 - cosine( t(rbind(cent,data)) )[1,-1]
#  dists[clust==i] <- d
#}

#dfdist <- data.frame(clust,dists,stringsAsFactors=FALSE)

#pdf(file.path(hmDir,"cluster_cosinedistances_afternorm.pdf"),width=9,h=6)
#p <- ggplot(dfdist, aes(dists, colour = factor(clust))) +
#  geom_density(size=2,alpha = 0.1)
#print(p)
#print(p + xlim(c(0,.3)))
#dev.off()


### Cluster ordering

# d <- 1-abs(cor(t(prop)))
# hc <- hclust(as.dist(d), method="average")
# o <- hc$order
# 
# dc <- 1-abs(cor(a[,-1]))
# co <- hclust(as.dist(dc), method="average")$order

# override the hierarchical clustering
o <- 1:nmetaclusts
co <- 1:length(scols)


### Prepare the cluster proportion data

rp <- sweep(prop, 1, STATS=rowMeans(prop), FUN="-")


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


testD <- list(rp[o,,drop=FALSE]*100, prop[o,,drop=FALSE]*100, a[o,-1][,co], aX[o,-1])

testD[[1]][testD[[1]]>max(cL[[1]]$breaks)] <- max(cL[[1]]$breaks)
testD[[1]][testD[[1]]<min(cL[[1]]$breaks)] <- min(cL[[1]]$breaks)
testD[[2]][testD[[2]]>max(cL[[2]]$breaks)] <- max(cL[[2]]$breaks)
testD[[2]][testD[[2]]<min(cL[[2]]$breaks)] <- min(cL[[2]]$breaks)


### Heatmaps will all the clusters
pdf(file.path(hmDir,"multiheatmap.pdf"),width=9,h=6)
multiHeatmap(testD,cL,xspace=2,ystarts=c(.25,.9,.925,.95,.98), clabelcex=.5)
dev.off()


### Heatmaps with larger clusters
w <- which(rowSums(prop>=.05)>=5)
oo <- o %in% w
testDs <- lapply(testD, function(u) u[oo,,drop=FALSE])

pdf(file.path(hmDir,"multiheatmap_subset.pdf"),width=9,h=3)
multiHeatmap(testDs,cL,xspace=2,ystarts=c(.25,.9,.925,.95,.98), clabelcex=.5)
dev.off()


### Read in the cluster merging file
if(file.exists("cluster_merging.xlsx")){
  
  cm <- read.xls("cluster_merging.xlsx")
  
  # remove spaces in labels bcs they are problematic...
  cm$label <- factor(cm$label, labels = gsub(" ", "_", levels(cm$label))) 
  
  clu <- cm$new_cluster
  names(clu) <- cm$old_cluster
  
  lablu <- unique(cm[,c("new_cluster","label")])
  (labels <- lablu$label)
  
  names(labels) <- lablu$new_cluster
  
  
  clustm <- clu[as.character(clust)] ### new merged cluster labels 
  
  
  am <- aggregate( el, by=list(clustm), FUN=median)
  aXm <- aggregate( elX, by=list(clustm), FUN=median)
  
  df <- data.frame(el,clust=labels[as.character(clustm)],samp)
  df2 <- melt(df,id=c("clust","samp"))
  
  ggp <- ggplot(df2, aes(x=value)) + geom_density(size=1,adjust=5) + 
    facet_grid(clust~variable,scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  pdf(file.path(hmDir,"clust_distros_merged.pdf"),w=ncol(el),h=length(unique(clustm)))
  print(ggp)
  dev.off()
  
  
  # manually remove cluster=0
  am <- am[am$Group.1>0,]
  aXm <- aXm[aXm$Group.1>0,]
  
  rownames(am) <- am$Group.1
  rownames(aXm) <- aXm$Group.1
  
  freqm <- table( cluster=clustm, samp )
  freqm <- freqm[rownames(freqm)>0,,drop=FALSE]
  propm <- t(t(freqm) / colSums(freqm))
  
  
  rpm <- sweep(propm, 1, STATS=rowMeans(propm), FUN="-")
  
  testDm <- list(rpm*100,propm*100,am[,-1],aXm[,-1])
  testDm <- lapply(testDm, function(u) {
    rownames(u) <- as.character(labels[rownames(u)])
    u
  })
  
  pdf(file.path(hmDir,"multiheatmap_merged.pdf"),width=15,h=3)
  multiHeatmap(testDm,cL,xspace=7,ystarts=c(.25,.9,.925,.95,.98), clabelcex=.5)
  dev.off()
  
  
  
  ### Save fcs files for clusters CD4 and CD8
  if(save_cd4_cd8){
    
    nrow_fcs <- sapply(fcs,nrow)
    
    clustmList <- split(labels[clustm], rep(names(nrow_fcs), nrow_fcs))
    
    writeOutCluster <- function(u,v,z,keep="CD4", outdir) {
      # u - cluster; v - flowFrame; z - original filename
      fn <- file.path(outdir, basename(z))
      cat(fn,"\n")
      out <- v[u==keep, ]
      write.FCS(out, fn)
    }
    
    m <- match(names(fcs), names(clustmList))
    clustmList <- clustmList[m]
    
    dummy <- mapply(writeOutCluster, u = clustmList, v = fcs, z = f, keep=keep_cd4, outdir = outdir_cd4)
    dummy <- mapply(writeOutCluster, u = clustmList, v = fcs, z = f, keep=keep_cd8, outdir = outdir_cd8)
    
    
  }
  
  
}

# ---------------------------------------
# tSNE analyses
# ---------------------------------------

dups <- duplicated(el)  # find duplicates
w <- which(!dups)

inds <- split(1:nrow(el), samp)  # create indices by sample

ts <- table(samp)
#ns <- pmin(ts,5000)  # per-sample, how many cells to downsample
(ns <- pmin(ts,1500))  # per-sample, how many cells to downsample

# get subsampled indices
subs <- mapply(function(u,v) {
  s <- sample(u,ns[v],replace = FALSE)
  intersect(s,w)
},inds,names(inds))

# spot check
# sapply(subs,length)

s <- unlist(subs)
set.seed(1234)
rtsne_out <- Rtsne(el[s,], pca = FALSE, verbose = TRUE)


### Plot of tsne - all cells, all clusters
df <- data.frame(X=rtsne_out$Y[,1],Y=rtsne_out$Y[,2],
  cluster=sprintf("%02d",clust[s]), sample=samp[s],
  group=gsub("[1-5]","",samp[s]))

ggp <- ggplot(df,  aes(x = X, y = Y, color=cluster)) +
  geom_point(size=1) +
  facet_wrap( ~ group) #+  scale_fill_brewer(palette="Spectral")

pdf(file.path(sneDir,"tSNE_all_clusters.pdf"),w=15,h=10)                 
print(ggp)
dev.off()




### Plot of tsne - cells that are close to the centers of clusters (no outliers), all clusters

a <- aggregate( el[s,], by=list(clust[s]), FUN=median)

dists <- distse <- rep(NA,length(s))
for(i in 1:nmetaclusts) {
  cat(i,"\n")
  cent <- a[i,-1,drop=FALSE]
  data <- el[s,][clust[s]==i,,drop=FALSE]
  d <- 1 - cosine( t(rbind(cent,data)) )[1,-1]
  dists[clust[s]==i] <- d
  cent <- matrix(rep( as.numeric(cent), nrow(data)), byrow=TRUE, ncol=ncol(cent))
  #distse[clust[s]==i] <- sqrt(rowSums((cent-data)^2))
  distse[clust[s]==i] <- rowMeans(abs(cent-data))
}



#dfdist <- data.frame(clust[s],dists,distse,stringsAsFactors=FALSE)

#p <- ggplot(dfdist, aes(distse, colour = factor(clust.s.))) +
#  geom_density(size=4,alpha = 0.1)
#print(p)

ggp <- ggplot(df[distse<.25,],  aes(x = X, y = Y, color=cluster)) +
  geom_point(size=1) + 
  facet_wrap( ~ group)

pdf(file.path(sneDir,"tSNE_all_clusters_filtered.pdf"),w=15,h=10)                 
print(ggp)
dev.off()



### Plot of tsne - large clusters

tc <- table(df$cluster)
kk <- names(tc[tc>1000])

df1 <- df[df$cluster %in% kk,]

ggp <- ggplot(df1,  aes(x = X, y = Y, color=cluster )) +
  geom_point(size=1) + scale_fill_brewer(palette="Spectral")+
  facet_wrap( ~ group)+ labs(x = "tSNE 1", y="tSNE 2")+ theme_bw()+
  theme(strip.text.x = element_text(size=15, face="bold"),
    axis.title.x  = element_text(size=15, face="bold"),
    axis.title.y  = element_text(size=15, face="bold"),
    panel.background = element_rect(fill = "white"))

pdf(file.path(sneDir,"tSNE_subset_clusters.pdf"),w=15,h=10)                 
print(ggp)
dev.off()



### Plot of tsne - merged clusters
if(file.exists("cluster_merging.xlsx")){
  
  df2 <- data.frame(X=rtsne_out$Y[,1],Y=rtsne_out$Y[,2],
    cluster=sprintf("%s",labels[as.character(clustm[s])]),sample=samp[s],
    clusterno=clustm[s],
    group=gsub("[1-5]","",samp[s]))
  
  k <- df2$cluster!="drop"
  df2 <- df2[k,]
  
  ggp <- ggplot(df2,  aes(x = X, y = Y, color=cluster)) +
    geom_point(size=1) + scale_fill_brewer(palette="Spectral")+
    facet_wrap( ~ group)+ labs(x = "tSNE 1", y="tSNE 2")+ theme_bw() +
    theme(strip.text.x = element_text(size=15, face="bold"),
      axis.title.x  = element_text(size=15, face="bold"),
      axis.title.y  = element_text(size=15, face="bold"),
      panel.background = element_rect(fill = "white"))
  
  pdf(file.path(sneDir,"tSNE_merged_clusters.pdf"),w=15,h=10)                 
  print(ggp)
  dev.off()
  
  
}




# ---------------------------------------
# plots and analysis of cluster frequencies
# ---------------------------------------


### Frequencies for 20 clusters 
p <- list()

condlu <- md$condition
names(condlu) <- md$shortname


for(i in 1:nrow(prop)) {
  
  df <- data.frame(prop=prop[i,], group=condlu[colnames(prop)])
  ds <- ddply(df, .(group), summarise, mean = mean(prop), sd = sd(prop))
  
  pd <- position_dodge(0.1)
  
  p[[i]] <- ggplot(df,aes(group,prop)) + 
    geom_point(size=3, shape = 1, stroke = 1) + 
    geom_point(data=ds, aes(x=group, y=mean), colour='blue', size=2) +
    geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean-sd, ymax=mean+sd),colour='blue', width=0.05) +
    ggtitle(paste0("cluster",rownames(prop)[i]))+ theme_bw()+
    ylab("proportion") + xlab("") +
    theme(axis.text.x  = element_text(size=15, face="bold"),
      axis.title.y  = element_text(size=15, face="bold"),
      panel.background = element_rect(fill = "white"))
  
}



pdf(file.path(frqDir,"cluster_frequences.pdf"),w=8,h=6,onefile=TRUE)
for (i in seq(length(p)))
  print(p[[i]])
dev.off()


#limma::plotMDS(prop,col=as.numeric(as.factor(gsub("[1-5]","",colnames(prop)))))


### Frequencies for merged clusters 

if(file.exists("cluster_merging.xlsx")){
  
  propmdf <- data.frame(cluster = rownames(testDm[[2]]), as.data.frame.matrix(testDm[[2]])) # testDm is used for heatmaps
  
  write.table(propmdf,file=file.path(frqDir,"frequencies_merged.xls"),
    row.names=FALSE,quote=FALSE,sep="\t")
  
  
  
  pm <- list()
  
  
  for(i in 1:nrow(propmdf)) {
    
    df <- data.frame(prop=unname(t(propmdf[i,-1])), group=condlu[colnames(prop)])
    ds <- ddply(df, .(group), summarise, mean = mean(prop), sd = sd(prop))
    
    pd <- position_dodge(0.1)
    
    pm[[i]] <- ggplot(df,aes(group,prop)) + 
      geom_point(size=3, shape = 1, stroke = 1) + 
      geom_point(data=ds, aes(x=group, y=mean), colour='blue', size=2) +
      geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean-sd, ymax=mean+sd),colour='blue', width=0.05) +
      ggtitle(rownames(propmdf)[i])+ theme_bw()+
      ylab("proportion") + xlab("") +
      theme(axis.text.x  = element_text(size=15, face="bold"),
        axis.title.y  = element_text(size=15, face="bold"),
        panel.background = element_rect(fill = "white"))
    
  }
  
  
  pdf(file.path(frqDir,"merged_cluster_frequences.pdf"),w=8,h=6,onefile=TRUE)
  for (i in seq(length(pm)))
    print(pm[[i]])
  dev.off()
  
  #limma::plotMDS(prop,col=as.numeric(as.factor(gsub("[1-5]","",colnames(prop)))))
  
}






# ---------------------------------------
# looking at changes in individual markers
# each cluster
# ---------------------------------------

mcols <- setdiff(cols,scols)

# re-extract the data as list of tables
es <- lapply(fcsT, function(u) {
  exprs(u)[,mcols]
})
e <- do.call("rbind",es)

ss1 <- gsub("Di","",colnames(e))
ss2 <- strsplit(ss1, "[a-zY]" )
panel_mass <- as.integer(sapply(ss2, tail, 1))

m <- match(panel_mass, panel$Isotope)
colnames(e) <- panel$Antigen[m]



a <- aggregate( e, by=list(clust,samp), FUN=median)
as <- split(a[,-1],a[,1])

pvals <- data.frame(cluster=numeric(),marker=numeric(),Pval=numeric(),stringsAsFactors=FALSE)

for(i in 1:length(as)) {
  g <- gsub("[1-5]","",as[[i]]$Group.2)
  thisp <- apply( as[[i]][,-1], 2, function(u) kruskal.test(split(u,g))$p.value)
  pvals <- rbind(pvals, data.frame(cluster=i,marker=names(thisp),Pval=thisp,stringsAsFactors = FALSE))
}

o <- order(pvals$Pval)
head(pvals[o,])

pvals <- subset(pvals,Pval<.1)

p <- list()

for(i in 1:nrow(pvals)) {
  k <- a$Group.1==pvals$cluster[i]
  g <- gsub("[1-5]","",a$Group.2[k])
  df <- data.frame(group=g,signal=a[k,pvals$marker[i]])
  p[[i]] <- ggplot(df,aes(group,signal)) + geom_jitter(width=.2) +
    ylab(pvals$marker[i]) + ggtitle(paste0("Cluster ",pvals$cluster[i]))+ theme_bw() + xlab("") +
    theme(axis.text.x  = element_text(size=15, face="bold"),
      axis.title.y  = element_text(size=15, face="bold"),
      panel.background = element_rect(fill = "white"))
}


pdf(file.path(crsDir,"cross_plots.pdf"),w=8,h=6,onefile=TRUE)
for (i in seq(length(p)))
  print(p[[i]])
dev.off()


f <- file.path(crsDir,"cross_pvals.xls")
write.table(pvals,f,row.names=FALSE,quote=FALSE,sep="\t")


sessionInfo()
