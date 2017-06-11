
Sys.time()


library(gdata)
library(limma)
library(UpSetR)

##############################################################################
# Test arguments
##############################################################################

prefix='23CD4TmemCD69_02CD4_'
outdir='../carsten_cytof/PD1_project/CK_2016-06-23_02_CD4_merging2_Tmem_merging2_CD69/090_cytokine_bimatrix'
path_data='../carsten_cytof/PD1_project/CK_2016-06-23_02_CD4_merging2_Tmem_merging2_CD69/010_data/23CD4TmemCD69_02CD4_expr_raw.rds'
path_panel_cytokines='../carsten_cytof/PD1_project/CK_panels/panel2CD4_23_cytokines.xlsx'


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
# Load expression data
# ------------------------------------------------------------

expr <- readRDS(path_data)

# ------------------------------------------------------------
# Load cytokine panel 
# ------------------------------------------------------------

# read panel, pick which columns to use
panel_cytokines <- read.xls(path_panel_cytokines, stringsAsFactors=FALSE)

if(!all(c("fcs_colname", "Isotope", "Antigen", "transform") %in% colnames(panel_cytokines)))
  stop("Wrong columns in the panel!!!")


cutoff_cytokines <- panel_cytokines[complete.cases(panel_cytokines), , drop = FALSE]
cutoff_cytokines


if(nrow(cutoff_cytokines) == 0)
  stop("No thresholds for cytokines in the panel file!!!")


cutoff_cytokines <- cutoff_cytokines[, c("fcs_colname", "positive_cutoff_raw_base", "positive_cutoff_raw_tx"), drop = FALSE]
colnames(cutoff_cytokines) <- c("fcs_colname", "base", "tx")
cutoff_cytokines


stopifnot(all(cutoff_cytokines$fcs_colname %in% colnames(expr)))


# ------------------------------------------------------------
# Calculate the bimatrix
# ------------------------------------------------------------


e_split <- split(expr[, cutoff_cytokines$fcs_colname, drop = FALSE], factor(expr$sample_id, levels = unique(expr$sample_id)))


bimatrix <- lapply(names(e_split), function(i){
  # i = "base_HD1"
  
  e <- e_split[[i]]
  
  day <- strsplit2(i, "_")[, 1]
  
  bm <- t(t(e) > cutoff_cytokines[, day])
  
  return(bm)
  
})


bimatrix <- do.call("rbind", bimatrix)
bimatrix <- apply(bimatrix, 2, as.numeric)

mm <- match(colnames(bimatrix), panel_cytokines$fcs_colname)
colnames(bimatrix) <- panel_cytokines$Antigen[mm]


bimatrixdf <- data.frame(cell_id = expr$cell_id, sample_id = expr$sample_id, bimatrix, check.names = FALSE)
head(bimatrixdf)


write.table(bimatrixdf, file = file.path(outdir, paste0(prefix, "bimatrix.txt")), row.names = FALSE, quote = FALSE, sep = "\t")
saveRDS(bimatrixdf, file = file.path(outdir, paste0(prefix, "bimatrix.rds")))


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

samp <- expr$sample_id

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










sessionInfo()




