##############################################################################
## <<02_flowsom_validation.R>>

# BioC 3.3
# Created 29 Sep 2016

##############################################################################
Sys.time()
##############################################################################

# Load packages
library(gdata)
library(FlowSOM)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2)
library(ConsensusClusterPlus)


##############################################################################
# Test arguments
##############################################################################

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
# flowsom_prefix='23_01_pca1_cl20_'
# flowsom_outdir='030_flowsom_validation'
# path_data='010_data/23_01_expr_raw.rds'
# path_clustering_observables='030_heatmaps/23_01_pca1_clustering_observables.xls'
# rand_seed_consensus=123
# nmetaclusts=20

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
rand_seed <- 1234

prefix <- flowsom_prefix
outdir <- flowsom_outdir

if(!file.exists(outdir)) 
  dir.create(outdir)

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

cell_id <- expr[, "cell_id"]
samp <- expr[, "sample_id"]
fcs_colnames <- colnames(expr)[!grepl("cell_id|sample_id", colnames(expr))]
e <- expr[, fcs_colnames]


# ------------------------------------------------------------
# Load clustering_observables
# ------------------------------------------------------------


clustering_observables <- read.table(path_clustering_observables, header = TRUE, sep = "\t", as.is = TRUE)

clust_observ <- clustering_observables[clustering_observables$clustering_observable, "mass"]

# selected columns for clustering 

scols <- which(fcs_colnames %in% clust_observ)

ef <- as.matrix(e[, scols])



# -------------------------------------
# Overwrite the function CDF from ConsensusClusterPlus so that it returns deltaK
# -------------------------------------


library(R.utils)

my_CDF <- function(ml,breaks=100){
  #plot CDF distribution
  plot(c(0),xlim=c(0,1),ylim=c(0,1),col="white",bg="white",xlab="consensus index",ylab="CDF",main="consensus CDF", las=2)
  k=length(ml)
  this_colors = rainbow(k-1)
  areaK = c()
  for (i in 2:length(ml)){
    v=triangle(ml[[i]],mode=1)
    
    #empirical CDF distribution. default number of breaks is 100    
    h = hist(v, plot=FALSE, breaks=seq(0,1,by=1/breaks))
    h$counts = cumsum(h$counts)/sum(h$counts)
    
    #calculate area under CDF curve, by histogram method.
    thisArea=0
    for (bi in 1:(length(h$breaks)-1)){
      thisArea = thisArea + h$counts[bi]*(h$breaks[bi+1]-h$breaks[bi]) #increment by height by width
      bi = bi + 1
    }
    areaK = c(areaK,thisArea)
    lines(h$mids,h$counts,col=this_colors[i-1],lwd=2,type='l')
  }
  legend(0.8,0.5,legend=paste(rep("",k-1),seq(2,k,by=1),sep=""),fill=this_colors)
  
  #plot area under CDF change.
  deltaK=areaK[1] #initial auc at k=2
  for(i in 2:(length(areaK))){
    #proportional increase relative to prior K.
    deltaK = c(deltaK,( areaK[i] - areaK[i-1])/areaK[i-1])
  }
  plot(1+(1:length(deltaK)),y=deltaK,xlab="k",ylab="relative change in area under CDF curve",main="Delta area",type="b")
  return(deltaK)
}


reassignInPackage("CDF", pkgName="ConsensusClusterPlus", my_CDF);




# -------------------------------------
# run FlowSOM with different rand seed 
# -------------------------------------

if(nmetaclusts == 20){
  
  rand_seed_list <- c(rand_seed, 1:10)
  rand_seed_consensus_list <- c(rand_seed_consensus, 1:10)
  
  results <- list()
  deltaK <- list()
  
  for(i in 1:length(rand_seed_list)){
    # i = 1
    
    rand_seed <- rand_seed_list[i]
    rand_seed_consensus <- rand_seed_consensus_list[i]
    
    # SOM
    set.seed(rand_seed)
    fsom <- FlowSOM::SOM(ef)
    
    # consensus clustering that is reproducible with seed
    data <- fsom$codes
    k <- nmetaclusts
    
    results[[i]] <- ConsensusClusterPlus::ConsensusClusterPlus(t(data),
      maxK = k, reps = 100, pItem = 0.9, pFeature = 1, title = tempdir(),
      plot = NULL, verbose = FALSE, clusterAlg = "hc", distance = "euclidean", seed = rand_seed_consensus)
    
    ### Generate the plot of relative change in area under CDF curve and get the values deltaK
    ml <- list()
    
    for(j in 2:k){
      ml[[j]] <- results[[i]][[j]]$ml
    }
    
    # pdf(file.path(outdir, paste0(prefix, "ConsensusClusterPlus_rs", rand_seed, "rsc", rand_seed_consensus, ".pdf")), width = 7, height = 7)
    
    deltaK[[i]] <- ConsensusClusterPlus:::CDF(ml)
    
    # dev.off()
    
    
  }
  
  
  
  ggdf <- data.frame(deltaK)
  colnames(ggdf) <- paste0("rsc", rand_seed_consensus_list)
  ggdf$no_clusters <- 2:k
  
  ggdfm <- melt(ggdf, id.vars = "no_clusters", variable.name = "rand_seed_consensus", value.name = "deltaK")
  
  ggdfm$rand_seed_consensus <- factor(ggdfm$rand_seed_consensus, levels = paste0("rsc", rev(rand_seed_consensus_list)))
  
  ggdfm$used_seed <- factor(ggdfm$rand_seed_consensus == paste0("rsc", rand_seed_consensus_list[1]))
  ggdfm$used_seed <- relevel(ggdfm$used_seed, ref = "TRUE")
  
  ggp <- ggplot(ggdfm, aes(x = no_clusters, y = deltaK, group = rand_seed_consensus, color = used_seed)) +
    geom_line(aes(linetype = used_seed)) +
    geom_point(shape = 1, size = 2)
  
  
  pdf(file.path(outdir, paste0(prefix, "ConsensusClusterPlus_CDFchange.pdf")), width = 8, height = 5)
  print(ggp)
  dev.off()
  
  
  
}


# -------------------------------------
# run FlowSOM to see if clusters are the same when running with nmetaclusts = 20 and taking k or running with nmetaclusts = k
# -------------------------------------

if(nmetaclusts < 20){
  
  # SOM
  set.seed(rand_seed)
  fsom <- FlowSOM::SOM(ef)
  
  # consensus clustering that is reproducible with seed
  data <- fsom$codes
  k <- nmetaclusts
  
  results <- ConsensusClusterPlus::ConsensusClusterPlus(t(data),
    maxK = k, reps = 100, pItem = 0.9, pFeature = 1, title = tempdir(),
    plot = NULL, verbose = FALSE, clusterAlg = "hc", distance = "euclidean", seed = rand_seed_consensus)
  
  # get cluster ids
  fsom_mc <- results[[k]]$consensusClass
  clust <- fsom_mc[fsom$mapping[,1]]
  

  k <- 20
  
  results20 <- ConsensusClusterPlus::ConsensusClusterPlus(t(data),
    maxK = k, reps = 100, pItem = 0.9, pFeature = 1, title = tempdir(),
    plot = NULL, verbose = FALSE, clusterAlg = "hc", distance = "euclidean", seed = rand_seed_consensus)
  
  # get cluster ids
  fsom_mc20 <- results20[[nmetaclusts]]$consensusClass
  clust20 <- fsom_mc20[fsom$mapping[,1]]
  
  out <- all(clust == clust20)
  
  write.table(out, file.path(outdir, paste0(prefix, "ConsensusClusterPlus_agreement.txt")), quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  
  
}











sessionInfo()













################################
### 02_flowsom_validation.R done!
################################