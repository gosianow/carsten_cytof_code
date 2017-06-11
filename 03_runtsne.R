

Sys.time()

# Load packages
library(Rtsne)


##############################################################################
# Test arguments
##############################################################################


prefix='23_01_pca1_'
outdir='CK_2016-06-23_01/040_tsnemaps'
path_data='CK_2016-06-23_01/010_data/23_01_expr_raw.rds'
path_clustering_observables='CK_2016-06-23_01/030_heatmaps/23_01_pca1_clustering_observables.xls'
tsne_pmin=1500


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



if(!file.exists(outdir)) 
  dir.create(outdir)

rand_seed <- 1234
perplexity <- 30

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------


expr <- readRDS(path_data)


cell_id <- expr[, "cell_id"]
samp <- expr[, "sample_id"]
fcs_colnames <- colnames(expr)[!grepl("cell_id|sample_id", colnames(expr))]
e <- expr[, fcs_colnames]


# ------------------------------------------------------------
# Load clustering_observables
# ------------------------------------------------------------

clustering_observables <- read.table(path_clustering_observables, header = TRUE, sep = "\t", as.is = TRUE)

clust_observ <- clustering_observables[clustering_observables$clustering_observable, "mass"]

# -------------------------------------

# Indeces of observables used for clustering 
scols <- which(fcs_colnames %in% clust_observ)


# ---------------------------------------
# tSNE analyses
# ---------------------------------------

et <- e[, scols]

### find duplicates
dups <- duplicated(et)  
w <- which(!dups)


### Data subsampling
# create indices by sample
inds <- split(1:length(samp), samp) 

# per-sample, how many cells to downsample
ts <- table(samp)
ns <- pmin(ts, tsne_pmin)  

# get subsampled indices
subs <- mapply(function(u,v) {
  set.seed(rand_seed)
  s <- sample(u, ns[v], replace = FALSE)
  intersect(s,w)
}, inds, names(inds))


cells2keep <- c(unlist(subs))


et_sub <- et[cells2keep, ]



### Run tSNE
set.seed(rand_seed)
rtsne_out <- Rtsne(et_sub, perplexity = perplexity, pca = FALSE, max_iter = 1000, verbose = TRUE, check_duplicates = FALSE)


# Save rtsne results

rtsne_data <- data.frame(cell_index = cells2keep, sample_name = samp[cells2keep], et_sub)

write.table(rtsne_data, file.path(outdir, paste0(prefix, "rtsne_data.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


saveRDS(rtsne_out, file = file.path(outdir, paste0(prefix, "rtsne_out.rds")))










sessionInfo()




