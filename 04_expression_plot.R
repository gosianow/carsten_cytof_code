

Sys.time()

# Load packages
library(gdata)
library(ggplot2)
library(reshape2)
library(plyr)
library(limma)

##############################################################################
# Test arguments
##############################################################################


prefix='23_01_pca1_cl20_all_'
outdir='../carsten_cytof/PD1_project/CK_2016-06-23_01/080_expression_auto'
path_metadata='../carsten_cytof/PD1_project/CK_metadata/metadata_23_01.xlsx'
path_expression='../carsten_cytof/PD1_project/CK_2016-06-23_01/080_expression_auto/23_01_pca1_cl20_all_expr.xls'
path_fun_plot_expression='00_plot_expression.R'


prefix='01CD4_23CD4merging5_29CD4merging5_clust_'
outdir='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/01_CD4_noHD10/080_expression/3responses_both'
path_expression=c('../carsten_cytof/PD1_project/CK_2016-06-23_01_CD4_mergingNEW2/080_expression/23CD4_01CD4_pca1_merging5_clust_expr.xls','../carsten_cytof/PD1_project/CK_2016-06-29_01_CD4_merging2/080_expression/29CD4_01CD4_pca1_merging5_clust_expr.xls')
path_metadata=c('../carsten_cytof/PD1_project/CK_metadata/metadata_23_01.xlsx','../carsten_cytof/PD1_project/CK_metadata/metadata_29_01_noHD10.xlsx')
path_fun_plot_expression='./00_plot_expression.R'


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
## replace _ with \n
levels(colors$condition) <- gsub("_", "\n", levels(colors$condition ))

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
# Load median expression data
# ------------------------------------------------------------


a <- lapply(1:length(path_expression), function(i){
  # i = 1
  a <- read.table(path_expression[i], header = TRUE, sep = "\t", as.is = TRUE, check.names = FALSE)
  a <- a[, -which(colnames(a) == "cluster")]
  a
})


## keep only those markers that are common for all the merged datasets 
marker_list <- unlist(lapply(a, colnames))
marker_list <- marker_list[!grepl("label|sample", marker_list)]

marker_table <- table(marker_list)
marker_overlap <- names(marker_table)[marker_table == length(path_expression)]

a <- rbind.fill(a)

a <- a[, colnames(a) %in% c(marker_overlap, "label", "sample"), drop = FALSE]


## drop the "drop" cluster
a <- a[a$label != "drop", , drop = FALSE]

## keep only those samples that are also in the metadata file
a <- a[a$sample %in% md$shortname, , drop = FALSE]


## keep only those clusters that are present in all the merged datasets
labels_keep <- names(which(table(a$label) == nrow(md)))

labels <- unique(a$label)
labels <- labels[labels %in% labels_keep]

labels <- data.frame(cluster = 1:length(labels), label = labels)
labels$label <- factor(labels$label, levels = unique(labels$label))
labels

a <- a[a$label %in% labels_keep, , drop = FALSE]

mm <- match(a$label, labels$label)
a <- cbind(cluster = labels$cluster[mm], a)
head(a)

markers_ordered <- colnames(a)[!colnames(a) %in% c("cluster", "label", "sample")]
markers_ordered


### Skipp samples that have NA for all the clusters
### Clusters with at least one NA will be skipped in the expression analysis (see 00_models.R)

table_sample <- table(a[complete.cases(a), "sample"])

keep_samps <- names(table_sample)[table_sample > 0]

a <- a[a$sample %in% keep_samps, , drop = FALSE]

md <- md[md$shortname %in% keep_samps, , drop = FALSE]


# Drop unused levels
md$response <- factor(md$response)
md$day <- factor(md$day)
md$data <- factor(md$data)
md$patient_id <- factor(md$patient_id)


# -----------------------------------------------------------------------------
# Plot expression per cluster
# -----------------------------------------------------------------------------


ggdf <- melt(a, id.vars = c("cluster", "label", "sample"), value.name = "expr", variable.name = "marker")

## use labels as clusters
ggdf$cluster <- factor(ggdf$label, levels = labels$label)


# add more info about samples
mm <- match(ggdf$samp, md$shortname)
ggdf$group <- factor(md$condition[mm])
ggdf$day <- factor(md$day[mm])
ggdf$data <- factor(md$data[mm])

## replace _ with \n
levels(ggdf$group) <- gsub("_", "\n", levels(ggdf$group))


## order markers as for heatmaps
ggdf$marker <- factor(ggdf$marker, levels = markers_ordered)



source(path_fun_plot_expression)

plot_expression(ggdf = ggdf, color_groups = color_groups, color_groupsb = color_groupsb, outdir = outdir, prefix = prefix)

















sessionInfo()




