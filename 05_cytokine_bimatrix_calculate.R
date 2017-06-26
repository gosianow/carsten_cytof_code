
Sys.time()


library(gdata)
library(limma)

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

if(!all(c("fcs_colname", "Isotope", "Antigen", "transform", "positive_cutoff_raw_base", "positive_cutoff_raw_tx") %in% colnames(panel_cytokines)))
  stop("Wrong columns in the cytokine panel!!!")


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










sessionInfo()




