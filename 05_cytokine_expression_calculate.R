

Sys.time()

library(gdata)

##############################################################################
# Test arguments
##############################################################################


prefix='23CD4TmemCD69_02CD4_'
outdir='../carsten_cytof/PD1_project/CK_2016-06-23_02_CD4_merging2_Tmem_merging2_CD69/090_cytokine_expression'
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

if( !file.exists(outdir) ) 
  dir.create(outdir, recursive = TRUE)


# ------------------------------------------------------------
# Load expression data
# ------------------------------------------------------------

expr <- readRDS(path_data)

fcs_colnames <- colnames(expr)[!grepl("cell_id|sample_id", colnames(expr))]
e <- expr[, fcs_colnames]

samp <- expr$sample_id

# ------------------------------------------------------------
# Load cytokine panel 
# ------------------------------------------------------------


# read panel, pick which columns to use
panel <- read.xls(path_panel_cytokines, stringsAsFactors=FALSE)

if(!all(c("fcs_colname", "Isotope", "Antigen", "transform") %in% colnames(panel)))
  stop("Wrong columns in panel!!!")


# ------------------------------------------------------------
# Get the overall median expression
# If sample has not enough cells, set expression to NA
# ------------------------------------------------------------

## Use Antigen as a marker name
mm <- match(colnames(e), panel$fcs_colname)
colnames(e) <- panel$Antigen[mm]

markers <- colnames(e)

### Median expression per cluster
a <- aggregate(e, by = list(sample = samp), FUN = median, drop = FALSE)
head(a)


### If sample has not enough cells, set expression to NA
min_cells <- 20
table_samp <- aggregate(e[, 1, drop = FALSE], by = list(sample = samp), FUN = length, drop = FALSE)

keep_samps <- table_samp[, 2] > min_cells

a[!keep_samps, markers] <- NA


### Add cluster and label columns to be able to use the 04_expression_analysis.R script
a$cluster <- 1
a$label <- "all"

a <- a[, c("cluster", "label", "sample", markers)]


### Save the median expression
write.table(a, file.path(outdir, paste0(prefix, "all_expr.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)














sessionInfo()






