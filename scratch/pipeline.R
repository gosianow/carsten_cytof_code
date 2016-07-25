# Load packages
library(flowCore)
library(gdata)
library(FlowSOM)
library(Repitools)
library(gplots)
library(Rtsne)
library(ggplot2)
library(plyr)
library(gplots)


#setwd("/Users/mark/Dropbox/Work/DB_projects/Mark_Carsten_CYTOF/CK_2015-10-18")
setwd("~/projects/cytof/carsten/Mark_Carsten_CYTOF/CK_2015-10-18")


# define directories
fcsDir <- "010_cleanfcs"; if( !file.exists(fcsDir) ) dir.create(fcsDir)
pcaDir <- "020_pcascores"; if( !file.exists(pcaDir) ) dir.create(pcaDir)
hmDir <- "030_heatmaps"; if( !file.exists(hmDir) ) dir.create(hmDir)
sneDir <- "040_tsnemaps"; if( !file.exists(sneDir) ) dir.create(sneDir)
frqDir <- "050_frequencies"; if( !file.exists(frqDir) ) dir.create(frqDir)
crsDir <- "060_crossplots"; if( !file.exists(crsDir) ) dir.create(crsDir)


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

# ugly way to get the mass of each column .. to match against panel
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
# run FlowSOM, make heatmaps
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
# -------------------------------------


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

mX <- match(panel_massX, panel$Isotope)
elX <- eX
colnames(elX) <- panel$Antigen[mX]

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


d <- 1-abs(cor(t(prop)))
hc <- hclust(as.dist(d), method="average")
o <- hc$order

dc <- 1-abs(cor(a[,-1]))
co <- hclust(as.dist(dc), method="average")$order

rp <- sweep(prop, 1, STATS=rowMeans(prop), FUN="-")


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

testD <- list(rp[o,]*100,prop[o,]*100,a[o,-1][,co],aX[o,-1])
testD[[1]][testD[[1]]>max(cL[[1]]$breaks)] <- max(cL[[1]]$breaks)
testD[[1]][testD[[1]]<min(cL[[1]]$breaks)] <- min(cL[[1]]$breaks)
testD[[2]][testD[[2]]>max(cL[[2]]$breaks)] <- max(cL[[2]]$breaks)
testD[[2]][testD[[2]]<min(cL[[2]]$breaks)] <- min(cL[[2]]$breaks)


pdf(file.path(hmDir,"multiheatmap.pdf"),width=9,h=6)
multiHeatmap(testD,cL,xspace=2,ystarts=c(.25,.9,.925,.95,.98), clabelcex=.5)
dev.off()

w <- which(rowSums(prop>=.10)>=2)
oo <- o %in% w
testDs <- lapply(testD, function(u) u[oo,])

pdf(file.path(hmDir,"multiheatmap_subset.pdf"),width=9,h=3)
multiHeatmap(testDs,cL,xspace=2,ystarts=c(.25,.9,.925,.95,.98), clabelcex=.5)
dev.off()








# ---------------------------------------
# tSNE analyses
# ---------------------------------------

dups <- duplicated(el)  # find duplicates
w <- which(!dups)

inds <- split(1:nrow(el),samp)  # create indices by sample

ts <- table(samp)
#ns <- pmin(ts,5000)  # per-sample, how many cells to downsample
(ns <- pmin(ts,1000))  # per-sample, how many cells to downsample

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


df <- data.frame(X=rtsne_out$Y[,1],Y=rtsne_out$Y[,2],
                 cluster=sprintf("%02d",clust[s]),sample=samp[s],
                 group=gsub("[1-5]","",samp[s]))

pdf(file.path(sneDir,"tSNE_all_clusters.pdf"),w=15,h=5)                 
ggplot(df,  aes(x = X, y = Y, color=cluster)) +
       geom_point(size=2) +
       facet_grid( ~ group) #+  scale_fill_brewer(palette="Spectral")
dev.off()


tc <- table(df$cluster)
kk <- names(tc[tc>500])

df1 <- df[df$cluster %in% kk,]
pdf(file.path(sneDir,"tSNE_subset_clusters.pdf"),w=15,h=5)                 
ggplot(df1,  aes(x = X, y = Y, color=cluster )) +
  geom_point(size=1.5) + scale_fill_brewer(palette="Spectral")+
  facet_grid( ~ group)+ labs(x = "tSNE 1", y="tSNE 2")+ theme_bw()+
  theme(strip.text.x = element_text(size=15, face="bold"),
        panel.background = element_rect(fill = "white"))
dev.off()


#df2 <- df1[!(df1$cluster %in% c("09","15")) & (df1$group %in% c("hPBMC","pPBMC")),]
df2 <- df1[(df1$group %in% c("hPBMC","pPBMC")),]
pdf(file.path(sneDir,"tSNE_subset_clusters_PBMC.pdf"),w=10,h=5)
ggplot(df2,  aes(x = X, y = Y, color=cluster )) +
  geom_point(size=1.5) + scale_fill_brewer(palette="Spectral")+
  facet_grid( ~ group)+ labs(x = "tSNE 1", y="tSNE 2")+ theme_bw()+
  theme(strip.text.x = element_text(size=15, face="bold"),
        panel.background = element_rect(fill = "white"))
        #panel.grid.major = NULL, panel.grid.minor=NULL) #+ #
dev.off()


#df3 <- df1[!(df1$cluster %in% c("09","15")) & (df1$group %in% c("hPMN","pPMN")),]
df3 <- df1[(df1$group %in% c("hPMN","pPMN")),]
pdf(file.path(sneDir,"tSNE_subset_clusters_PMN.pdf"),w=10,h=5)
ggplot(df3,  aes(x = X, y = Y, color=cluster )) +
  geom_point(size=1.5) + scale_fill_brewer(palette="Spectral")+
  facet_grid( ~ group)+ labs(x = "tSNE 1", y="tSNE 2")+ theme_bw()+
  theme(strip.text.x = element_text(size=15, face="bold"),
        panel.background = element_rect(fill = "white"))
dev.off()







# ---------------------------------------
# plots and analysis of cluster frequencies
# ---------------------------------------


p <- list()

for(i in 1:nrow(prop)) {
 
  df <- data.frame(prop=prop[i,], group=md$condition)
  ds <- ddply(df, .(group), summarise, mean = mean(prop), sd = sd(prop))

  pd <- position_dodge(0.1)

  p[[i]] <- ggplot(df,aes(group,prop)) + 
    geom_point() + 
    geom_point(data=ds, aes(x=group, y=mean), colour='blue', size=4) +
    geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean-sd, ymax=mean+sd),colour='blue', width=0.2) +
    ggtitle(paste0("cluster",rownames(prop)[i]))
  
}


limma::plotMDS(prop,col=as.numeric(as.factor(gsub("[1-5]","",colnames(prop)))))



pdf(file.path(frqDir,"cluster_frequences.pdf"),w=8,h=6,onefile=TRUE)
for (i in seq(length(p)))
  print(p[[i]])
dev.off()


kt <- apply(prop,1, function(u) kruskal.test(split(u,md$condition))$p.value)

propm <- matrix(as.numeric(prop),nrow=nrow(prop),ncol=ncol(prop),dimnames=dimnames(prop))
propm <- data.frame(cluster=rownames(propm),round(propm,3),Pvalue=kt)

write.table(propm, file=file.path(frqDir,"cluster_freqs.xls"), 
            sep="\t", row.names=FALSE, quote=FALSE)






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
      ylab(pvals$marker[i]) + ggtitle(paste0("Cluster ",pvals$cluster[i]))
}


pdf(file.path(crsDir,"cross_plots.pdf"),w=8,h=6,onefile=TRUE)
for (i in seq(length(p)))
  print(p[[i]])
dev.off()


f <- file.path(crsDir,"cross_pvals.xls")
write.table(pvals,f,row.names=FALSE,quote=FALSE,sep="\t")



