

Sys.time()


##############################################################################
# Test arguments
##############################################################################


path_pca_score='CK_2016-06-23_01/020_pcascores/23_01_princompscore_by_sample.xls'
path_pca_score_cutoff='CK_2016-06-23_01/010_helpfiles/23_01_pca1_princompscore_cutoff.txt'
prefix='23_01_pca1_'
outdir='CK_2016-06-23_01/030_heatmaps'


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

# ------------------------------------------------------------
# define directories
# ------------------------------------------------------------

if(!file.exists(outdir)) 
  dir.create(outdir)

# ------------------------------------------------------------
# Load PCA scores
# ------------------------------------------------------------

prs <- read.table(path_pca_score, header = TRUE, sep = "\t", as.is = TRUE)

prs <- prs[order(prs$avg_score, decreasing = TRUE), ]

# ------------------------------------------------------------
# Load the cutoff
# ------------------------------------------------------------


pca_score_cutoff <- read.table(path_pca_score_cutoff)[1,1]


# -------------------------------------
# selected observables with PCA scores higher than a cutoff but skipp top X scores 
# -------------------------------------


clustering_observables <- prs[, c("mass", "marker", "avg_score")]

clustering_observables$clustering_observable <- prs$avg_score > pca_score_cutoff


### Save the observables

write.table(clustering_observables, file = file.path(outdir, paste0(prefix, "clustering_observables.xls")), row.names=FALSE, quote=FALSE, sep="\t")




sessionInfo()









