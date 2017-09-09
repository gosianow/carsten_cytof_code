
Sys.time()

# Load packages
library(flowCore)
library(gdata)
library(ggplot2)
library(reshape2)
library(limma)
library(dplyr)
library(ggrepel)
library(plyr)

##############################################################################
# Test arguments
##############################################################################

args <- NULL

path_data='../carsten_cytof/PD1_project/CK_2016-06-23_01/010_data/23_01_expr_raw.rds'
path_metadata='../carsten_cytof/PD1_project/CK_metadata/metadata_23_01.xlsx'
prefix='23_01_'
suffix='_raw'
outdir='../carsten_cytof/PD1_project/CK_2016-06-23_01/010_data'



path_data=c('../carsten_cytof/PD1_project/CK_2016-06-23_02/010_data/23_02v2_expr_raw.rds','../carsten_cytof/PD1_project/CK_2016-06-29_02/010_data/29_02v2_expr_raw.rds')
path_metadata=c('../carsten_cytof/PD1_project/CK_metadata/metadata_23_02.xlsx','../carsten_cytof/PD1_project/CK_metadata/metadata_29_02.xlsx')
prefix='02v2_23_29_'
suffix='_raw'
outdir='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/02v2/010_data'



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

expr <- lapply(path_data, readRDS)

expr <- rbind.fill(expr)

samp <- expr[, "sample_id"]
fcs_colnames <- colnames(expr)[!grepl("cell_id|sample_id", colnames(expr))]
e <- expr[, fcs_colnames]


# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

md <- lapply(1:length(path_metadata), function(i){
  md <- read.xls(path_metadata[i], stringsAsFactors = FALSE)
  md
})

md <- plyr::rbind.fill(md)
rownames(md) <- md$shortname

### Factor arrangment
md$response <- factor(md$response, levels = c("NR", "R", "HD"))
md$response <- factor(md$response)
md$day <- factor(md$day, levels = c("base", "tx"))
md$day <- factor(md$day)
md$patient_id <- factor(md$patient_id)
md$data <- factor(md$data)
md$data_day <- interaction(md$data, md$day, lex.order = TRUE, drop = TRUE)

### Colors 
colors <- unique(md[, c("condition", "color")])
colors$condition <- factor(colors$condition)

color_groups <- colors$color
names(color_groups) <- colors$condition

color_groupsb <- adjustcolor(color_groups, alpha = 0.3)
names(color_groupsb) <- colors$condition

color_samples <- md$color
names(color_samples) <- md$shortname

colors <- unique(md[, c("response", "color")])
color_response <- colors$color
names(color_response) <- colors$response



# ------------------------------------------------------------
# Plot the MDS plot
# ------------------------------------------------------------


df <- data.frame(samp = samp, e, check.names = FALSE)


expr_median_sample_tbl <- df %>%
  group_by(samp) %>% 
  summarize_all(funs(median))



expr_median_sample <- t(expr_median_sample_tbl[, -1])
colnames(expr_median_sample) <- expr_median_sample_tbl$samp


mds <- plotMDS(expr_median_sample)


ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y, 
  samp = colnames(expr_median_sample))


mm <- match(ggdf$samp, md$shortname)

ggdf$condition <- md$condition[mm]
ggdf$response <- md$response[mm]
ggdf$data <- md$data[mm]
ggdf$day <- md$day[mm]
ggdf$data_day <- md$data_day[mm]



ggp <- ggplot(ggdf, aes(x = MDS1, y = MDS2, color = response, shape = data)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_label_repel(aes(label = samp)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, face = "bold", color = "black"), 
        axis.title = element_text(size = 16, face = "bold", color = "black"), 
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  scale_color_manual(values = color_response) 


pdf(file.path(outdir, paste0(prefix, "mds", suffix,".pdf")), width = 10, height = 8)
print(ggp)
dev.off()








sessionInfo()


