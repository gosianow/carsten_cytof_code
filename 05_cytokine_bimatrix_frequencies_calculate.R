
Sys.time()


library(gdata)
library(limma)
library(UpSetR)
library(plyr)

##############################################################################
# Test arguments
##############################################################################

prefix='23CD4TmemCD69_02CD4_'
outdir='../carsten_cytof/PD1_project/CK_2016-06-23_02_CD4_merging3_Tmem_merging1_CD69/090_cytokine_bimatrix'
path_bimatrix='../carsten_cytof/PD1_project/CK_2016-06-23_02_CD4_merging3_Tmem_merging1_CD69/090_cytokine_bimatrix/23CD4TmemCD69_02CD4_bimatrix.rds'

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
  dir.create(outdir, recursive = TRUE)


# ------------------------------------------------------------
# Load bimatrix data
# ------------------------------------------------------------

bm <- lapply(1:length(path_bimatrix), function(i){
  # i = 1
  bm_tmp <- readRDS(path_bimatrix[i])
  
})


bimatrixdf <- rbind.fill(bm)

table(complete.cases(t(bimatrixdf)))

bimatrixdf <- bimatrixdf[, complete.cases(t(bimatrixdf)), drop = FALSE]

bimatrix <- bimatrixdf[, !grepl("cell_id|sample_id", colnames(bimatrixdf)), drop = FALSE]


# ------------------------------------------------------------
# Upsetr plots
# ------------------------------------------------------------


pdf(file.path(outdir, paste0(prefix, "upsetr.pdf")), w = 16, h = 6)
upset(data.frame(bimatrix, check.names = FALSE), sets = colnames(bimatrix), nintersects = 50, order.by = "freq")
dev.off()


# ------------------------------------------------------------
# Create a table with clustering observables - needed to run the 02_flowsom.R script
# ------------------------------------------------------------

clustering_observables <- data.frame(mass = colnames(bimatrix), marker = colnames(bimatrix), clustering_observable = TRUE, stringsAsFactors = FALSE)
clustering_observables

write.table(clustering_observables, file.path(outdir, paste0(prefix, "clustering_observables.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# --------------------------------------------------------------------------
# Calculate the positive marker frequencies per marker
# --------------------------------------------------------------------------

samp <- bimatrixdf$sample_id

freq_list <- lapply(colnames(bimatrix), function(i){
  freq_tmp <- table(bimatrix[, i], samp)
  rownames(freq_tmp) <- paste0(i, c("-", "+"))
  freq_tmp
})


prop_list <- lapply(freq_list, function(x){
  prop_tmp <- t(t(x[2, , drop = FALSE]) / colSums(x, na.rm = TRUE)) * 100
})


prop <- do.call(rbind, prop_list)
prop_out <- data.frame(cluster = rownames(prop), label = rownames(prop), prop)


freq <- do.call(rbind, freq_list)
freq_out <- data.frame(cluster = rownames(freq), label = rownames(freq), freq)


write.table(prop_out, file=file.path(outdir, paste0(prefix, "frequencies.xls")), row.names=FALSE, quote=FALSE, sep="\t")
write.table(freq_out, file=file.path(outdir, paste0(prefix, "counts.xls")), row.names=FALSE, quote=FALSE, sep="\t")



if(length(path_bimatrix) > 1){
  write.table(bimatrixdf, file = file.path(outdir, paste0(prefix, "bimatrix.txt")), row.names = FALSE, quote = FALSE, sep = "\t")
  saveRDS(bimatrixdf, file = file.path(outdir, paste0(prefix, "bimatrix.rds")))
  
}






sessionInfo()




