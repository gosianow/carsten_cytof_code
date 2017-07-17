

Sys.time()


##############################################################################
# Test arguments
##############################################################################


prefix='23_01_pca1_cl20_'
outdir='../carsten_cytof/PD1_project/CK_2016-06-23_01/080_expression_auto'
path_data='../carsten_cytof/PD1_project/CK_2016-06-23_01/010_data/23_01_expr_raw.rds'
path_clustering_observables='../carsten_cytof/PD1_project/CK_2016-06-23_01/030_heatmaps/23_01_pca1_clustering_observables.xls'
path_clustering='../carsten_cytof/PD1_project/CK_2016-06-23_01/030_heatmaps/23_01_pca1_cl20_clustering.xls'
path_clustering_labels='../carsten_cytof/PD1_project/CK_2016-06-23_01/030_heatmaps/23_01_pca1_cl20_clustering_labels.xls'

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

if( !file.exists(outdir) ) 
  dir.create(outdir, recursive = TRUE)


# ------------------------------------------------------------
# Load expression data
# ------------------------------------------------------------

expr <- readRDS(path_data)

fcs_colnames <- colnames(expr)[!grepl("cell_id|sample_id", colnames(expr))]
e <- expr[, fcs_colnames]


# ------------------------------------------------------------
# Load clustering data
# ------------------------------------------------------------

# clustering observables
clustering_observables <- read.table(path_clustering_observables, header = TRUE, sep = "\t", as.is = TRUE)
rownames(clustering_observables) <- clustering_observables$mass

clust_observ <- clustering_observables[clustering_observables$clustering_observable, "mass"]


# clustering labels
labels <- read.table(path_clustering_labels, header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))
rownames(labels) <- labels$cluster


# clustering
clustering <- read.table(path_clustering, header = TRUE, sep = "\t", as.is = TRUE)

## drop the "drop" cluster
if("drop" %in% labels$label){
  
  clust2drop <- labels$cluster[labels$label == "drop"]
  cells2drop <- clustering$cluster != clust2drop
  clustering <- clustering[cells2drop, , drop = FALSE]
  labels <- labels[labels$label != "drop", ,drop = FALSE]
  labels$label <- factor(labels$label)
  e <- e[cells2drop, , drop = FALSE]
  
}


clust <- clustering[, "cluster"]
samp <- clustering[, "sample_id"]


# ------------------------------------------------------------
# Get the isotope and antigen for fcs markers

m <- match(fcs_colnames, clustering_observables$mass)

fcs_panel <- data.frame(colnames = fcs_colnames, Isotope = clustering_observables$mass[m], Antigen = clustering_observables$marker[m], stringsAsFactors = FALSE)


# Indeces of observables used for clustering 
scols <- which(fcs_colnames %in% clust_observ)

# Indeces of other observables
xcols <- which(!fcs_colnames %in% clust_observ)


# Ordered by decreasing pca score
if("avg_score" %in% colnames(clustering_observables)){
  scols <- scols[order(clustering_observables[fcs_colnames[scols], "avg_score"], decreasing = TRUE)]
  xcols <- xcols[order(clustering_observables[fcs_colnames[xcols], "avg_score"], decreasing = TRUE)]
}


# ------------------------------------------------------------
# Get the median expression per cluster and overall
# If sample has not enough cells, set expression to NA
# ------------------------------------------------------------

colnames(e) <- fcs_panel$Antigen


### Median expression per cluster

min_cells <- 5
table_samp <- aggregate(e[, 1, drop = FALSE], by = list(cluster = clust, sample = samp), FUN = length, drop = FALSE)
keep_samps <- table_samp[, 3] > min_cells


a <- aggregate(e, by = list(cluster = clust, sample = samp), FUN = median, drop = FALSE)
head(a)

mm <- match(a$cluster, labels$cluster)
a$label <- labels$label[mm]

a <- a[, c("cluster", "label", "sample", fcs_panel$Antigen[c(scols, xcols)])]

a[!keep_samps, fcs_panel$Antigen[c(scols, xcols)]] <- NA


### Save the median expression per cluster and sample
write.table(a, file.path(outdir, paste0(prefix, "clust_expr.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)




### Overall median expression


min_cells <- 50
table_samp <- table(samp)
keep_samps <- names(table_samp)[which(table_samp > min_cells)]


a <- aggregate(e, by = list(sample = samp), FUN = median)
a$cluster <- -1
a$label <- "all"

a <- a[, c("cluster", "label", "sample", fcs_panel$Antigen[c(scols, xcols)])]

a[!a$sample %in% keep_samps, fcs_panel$Antigen[c(scols, xcols)]] <- NA


### Save the median expression per cluster and sample
write.table(a, file.path(outdir, paste0(prefix, "all_expr.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)














sessionInfo()






