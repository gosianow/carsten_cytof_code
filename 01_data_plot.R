
Sys.time()

# Load packages
library(flowCore)
library(gdata)
library(ggplot2)
library(reshape2)
library(limma)

##############################################################################
# Test arguments
##############################################################################

args <- NULL

path_data='../carsten_cytof/PD1_project/CK_2016-06-23_01/010_data/23_01_expr_raw.rds'
path_panel='../carsten_cytof/PD1_project/CK_panels/panel1.xlsx'
path_metadata='../carsten_cytof/PD1_project/CK_metadata/metadata_23_01.xlsx'
prefix='23_01_'
suffix='_raw'
outdir='../carsten_cytof/PD1_project/CK_2016-06-23_01/010_data'

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
md <- gdata::read.xls(path_metadata, stringsAsFactors=FALSE)

# add more info about samples
cond_split <- limma::strsplit2(md$condition, "_")
colnames(cond_split) <- c("day", "response")

md[, c("day", "response")] <- cond_split
md$response <- factor(md$response, levels = c("NR", "R", "HD"))
md$response <- factor(md$response)
md$day <- factor(md$day, levels = c("base", "tx"))
md$day <- factor(md$day)
md$patient_id <- factor(md$patient_id)


## colors per sample
color_samples <- md$color
names(color_samples) <- md$shortname


# ------------------------------------------------------------
# Load panel
# ------------------------------------------------------------

# read panel, pick which columns to use
panel <- gdata::read.xls(path_panel, stringsAsFactors=FALSE)

if(!all(c("fcs_colname", "Isotope", "Antigen", "transform") %in% colnames(panel)))
  stop("Wrong columns in panel!!!")


# ------------------------------------------------------------
# Plot expression of markers stratified per sample
# ------------------------------------------------------------


df <- data.frame(samp = samp, e, check.names = FALSE)
dfm <- melt(df, id.var = "samp")

mm <- match(dfm$samp, md$shortname)
dfm$day <- md$day[mm]

mm <- match(dfm$variable, panel$fcs_colname)
dfm$Antigen <- panel$Antigen[mm]

dfm$marker <- paste0(dfm$variable, " / ", dfm$Antigen)



ggp <- ggplot(dfm, aes(x = value, color = samp, linetype = day)) +
  geom_density(adjust = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title = element_blank(), legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2)) +
  scale_color_manual(values = color_samples) +
  facet_wrap(~ marker, nrow = 4, scales = "free")


pdf(file.path(outdir, paste0(prefix, "distrosgrp", suffix,".pdf")), w = ncol(e)*2/3, h = 11)
print(ggp)
dev.off()



# ------------------------------------------------------------
# Plot number of cells per sample
# ------------------------------------------------------------

samp <- factor(samp, levels = md$shortname)

cell_table <- table(samp)

cell_counts <- rep(0, nrow(md))
names(cell_counts) <- md$shortname

cell_counts[names(cell_table)] <- as.numeric(cell_table)


ggdf <- data.frame(sample_id = factor(names(cell_counts), levels = md$shortname), cell_counts = cell_counts)


ggp <- ggplot(ggdf, aes(x = sample_id, y = cell_counts, fill = sample_id)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = cell_counts), hjust=0.5, vjust=-0.5, size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none") +
  scale_fill_manual(values = color_samples, drop=FALSE) +
  scale_x_discrete(drop=FALSE)


pdf(file.path(outdir, paste0(prefix, "cell_counter.pdf")), w = nlevels(ggdf$sample_id)/3 + 2, h = 5)
print(ggp)
dev.off()







sessionInfo()


