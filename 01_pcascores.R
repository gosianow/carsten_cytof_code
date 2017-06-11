
Sys.time()

# Load packages
library(gdata)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2)
library(coop)
library(RColorBrewer)

##############################################################################
# Test arguments
##############################################################################

path_data='CK_2016-06-23_01/010_data/23_01_expr_raw.rds'
path_panel='CK_panels/panel1.xlsx'
path_metadata='CK_metadata/metadata_23_01.xlsx'
prefix='23_01_'
outdir='CK_2016-06-23_01/020_pcascores'

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
# Load expression
# ------------------------------------------------------------

expr <- readRDS(path_data)

samp <- expr[, "sample_id"]
fcs_colnames <- colnames(expr)[!grepl("cell_id|sample_id", colnames(expr))]
e <- expr[, fcs_colnames]


# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

# read metadata
md <- read.xls(path_metadata, stringsAsFactors=FALSE)


# ------------------------------------------------------------
# Load panel
# ------------------------------------------------------------

# read panel, pick which columns to use
panel <- read.xls(path_panel, stringsAsFactors=FALSE)

if(!all(c("fcs_colname", "Isotope", "Antigen", "transform") %in% colnames(panel)))
  stop("Wrong columns in panel!!!")


# get the isotope and antigen for fcs markers
m <- match(fcs_colnames, panel$fcs_colname)

fcs_panel <- data.frame(fcs_colname = fcs_colnames, Isotope = panel$Isotope[m], Antigen = panel$Antigen[m], stringsAsFactors = FALSE)


# --------------------------------------------------------------------------
# run Levine et al. 2015 marker scoring,
# --------------------------------------------------------------------------

min_samp <- nrow(md)
  

doPRINCOMP <- function(z, ncomp=3) {
  # z = es[[1]]; ncomp=3
  
  if(nrow(z) < min_samp)
    return(rep(NA, ncol(z)))
  
  ## Score by Levine
  pr <- prcomp(z, center = TRUE, scale. = FALSE)  # default is to center but not scale

  score <- rowSums (outer( rep(1, ncol(z)), pr$sdev[1:ncomp]^2 ) * abs(pr$rotation[,1:ncomp]) )

  return(score)

}


# split expression per sample
es <- split(e, samp)

prs <- sapply(es, doPRINCOMP)
rmprs <- rowMeans(prs, na.rm = TRUE)

prs <- data.frame(mass = fcs_panel$fcs_colname, marker = fcs_panel$Antigen, avg_score = rmprs, round(prs, 4))

### Save ordered PCA scores
o <- order(rmprs, decreasing=TRUE)
prs <- prs[o,]

write.table(prs, file = file.path(outdir, paste0(prefix, "princompscore_by_sample.xls")), sep="\t", row.names=FALSE, quote=FALSE)



# --------------------------------------------------------------------------
# Plot 
# --------------------------------------------------------------------------


## plot PCA scores

ggp1 <- ggplot(prs, aes(x = marker, y = avg_score)) +
  geom_bar(stat="identity") +
  theme_bw() +
  ylab("Average PCA score") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12))


pdf(file.path(outdir, paste0(prefix, "princompscore_average.pdf")), width = 10, height = 7)
print(ggp1)
dev.off()


## plot PCA scores for ordered markers

prs$marker <- factor(prs$marker, levels = prs$marker)

ggp2 <- ggplot(prs, aes(x = marker, y = avg_score)) +
  geom_bar(stat="identity") +
  theme_bw() +
  ylab("Average PCA score") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12))


pdf(file.path(outdir, paste0(prefix, "princompscore_average_ord.pdf")), width = 10, height = 7)
print(ggp2)
dev.off()






sessionInfo()





