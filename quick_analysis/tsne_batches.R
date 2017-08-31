

Sys.time()

# Load packages
library(Rtsne)
library(plyr)

##############################################################################
# Test arguments
##############################################################################


prefix='2329_02_pca0_'
outdir='../tsne_batches'
path_data=c('../carsten_cytof/PD1_project/CK_2016-06-23_02/010_data/23_02_expr_raw.rds', '../carsten_cytof/PD1_project/CK_2016-06-29_02/010_data/29_02_expr_raw.rds')
path_clustering_observables=c('../carsten_cytof/PD1_project/CK_2016-06-23_02/030_heatmaps/23_02v2_pca0_clustering_observables.xls', '../carsten_cytof/PD1_project/CK_2016-06-29_02/030_heatmaps/29_02v2_pca0_clustering_observables.xls')
tsne_pmin=1500


if(!file.exists(outdir)) 
  dir.create(outdir)

rand_seed <- 1234
perplexity <- 30


# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

expr <- lapply(path_data, readRDS)

expr <- rbind.fill(expr)

cell_id <- expr[, "cell_id"]
samp <- expr[, "sample_id"]
fcs_colnames <- colnames(expr)[!grepl("cell_id|sample_id", colnames(expr))]
e <- expr[, fcs_colnames]


# ------------------------------------------------------------
# Load clustering_observables
# ------------------------------------------------------------

clustering_observables <- lapply(path_clustering_observables, read.table, header = TRUE, sep = "\t", as.is = TRUE)

clustering_observables <- Reduce(function(...) merge(..., by = c("marker", "mass", "clustering_observable"), all=TRUE, sort = FALSE), clustering_observables)


clust_observ <- clustering_observables[clustering_observables$clustering_observable, "mass"]


# -------------------------------------

# Indeces of observables used for clustering 
scols <- which(fcs_colnames %in% clust_observ)


# ---------------------------------------
# tSNE analyses
# ---------------------------------------

et <- e[, scols]

### find duplicates
dups <- duplicated(et)  
w <- which(!dups)


### Data subsampling
# create indices by sample
inds <- split(1:length(samp), samp) 

# per-sample, how many cells to downsample
ts <- table(samp)
ns <- pmin(ts, tsne_pmin)  

# get subsampled indices
subs <- mapply(function(u,v) {
  set.seed(rand_seed)
  s <- sample(u, ns[v], replace = FALSE)
  intersect(s,w)
}, inds, names(inds))


cells2keep <- c(unlist(subs))


et_sub <- et[cells2keep, ]



### Run tSNE
set.seed(rand_seed)
rtsne_out <- Rtsne(et_sub, perplexity = perplexity, pca = FALSE, max_iter = 1000, verbose = TRUE, check_duplicates = FALSE)


# Save rtsne results

rtsne_data <- data.frame(cell_index = cells2keep, sample_name = samp[cells2keep], et_sub)

write.table(rtsne_data, file.path(outdir, paste0(prefix, "rtsne_data.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


saveRDS(rtsne_out, file = file.path(outdir, paste0(prefix, "rtsne_out.rds")))




##############################################################################


Sys.time()

# Load packages
library(gdata)
library(Rtsne)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(limma) 


path_metadata=c('../carsten_cytof/PD1_project/CK_metadata/metadata_23_02.xlsx', '../carsten_cytof/PD1_project/CK_metadata/metadata_29_02.xlsx')

path_clustering=c('../carsten_cytof/PD1_project/CK_2016-06-23_02/030_heatmaps/23_02v2_pca0_merging3_clustering.xls', '../carsten_cytof/PD1_project/CK_2016-06-29_02/030_heatmaps/29_02v2_pca0_merging3_clustering.xls')

path_clustering_labels=c('../carsten_cytof/PD1_project/CK_2016-06-23_02/030_heatmaps/23_02v2_pca0_merging3_clustering_labels.xls', '../carsten_cytof/PD1_project/CK_2016-06-29_02/030_heatmaps/29_02v2_pca0_merging3_clustering_labels.xls')



suffix <- ""
pdf_width=22
pdf_height=7


# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

md <- lapply(path_metadata, read.xls, stringsAsFactors=FALSE)

md <- rbind.fill(md)

# add more info about samples
cond_split <- strsplit2(md$condition, "_")
colnames(cond_split) <- c("day", "response")

md[, c("day", "response")] <- cond_split
md$response <- factor(md$response, levels = c("NR", "R", "HD"))
md$response <- factor(md$response)
md$day <- factor(md$day, levels = c("base", "tx"))
md$day <- factor(md$day)
md$patient_id <- factor(md$patient_id)


# ------------------------------------------------------------
# Load cluster data
# ------------------------------------------------------------

## clustering labels

labels <- lapply(path_clustering_labels, read.table, header = TRUE, sep = "\t", as.is = TRUE)

labels <- Reduce(function(...) merge(..., by = c("cluster", "label"), all=TRUE, sort = FALSE), labels)

labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))
labels


## clustering

clustering <- lapply(path_clustering, read.table, header = TRUE, sep = "\t", as.is = TRUE)

clustering <- rbind.fill(clustering)

clust <- clustering[, "cluster"]


# ------------------------------------------------------------
# Colors for tSNE maps
# ------------------------------------------------------------


# ggplot palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=60 , c=100)[1:n]
}

# color blind palette

colors_muted <- c("#DC050C", "#E8601C", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6", "#F1932D", "#F6C141", "#F7EE55", "#4EB265", "#90C987", "#CAEDAB")
color_ramp <- c(colors_muted, gg_color_hue(max(1, nlevels(labels$label) - length(colors_muted))))

colors_tsne <- color_ramp[1:nlevels(labels$label)]
names(colors_tsne) <- levels(labels$label)


# ------------------------------------------------------------
# Load tSNE data
# ------------------------------------------------------------

rtsne_data$cell_index2 <- paste0(rtsne_data$cell_index, rtsne_data$sample_name) 
clustering$cell_id2 <- paste0(1:nrow(clustering), clustering$sample_id)


## keep the tSNE cells that were used in clustering
cells2keep_rtsne <- rtsne_data$cell_index2 %in% clustering[, "cell_id2"]
table(cells2keep_rtsne)

# ------------------------------------------------------------
# tSNE plots
# ------------------------------------------------------------

# get clustering for cells that were used in tSNE
names(clust) <- clustering[, "cell_id2"]
clust_tsne <- clust[as.character(rtsne_data$cell_index2[cells2keep_rtsne])]

ggdf <- data.frame(tSNE1 = rtsne_out$Y[cells2keep_rtsne, 1], tSNE2 = rtsne_out$Y[cells2keep_rtsne, 2], cluster = clust_tsne, sample = rtsne_data$sample_name[cells2keep_rtsne])

# add group info
mm <- match(ggdf$sample, md$shortname)
ggdf$group <- md$response[mm]


# use cluster labels instead of numbers
mm <- match(ggdf$cluster, labels$cluster)
ggdf$cluster <- labels$label[mm]
ggdf$cluster <- factor(ggdf$cluster, levels = levels(labels$label))


# skipp the "drop" cluster

skipp_drop <- ggdf$cluster != "drop"

ggdf <- ggdf[skipp_drop, ]
ggdf$cluster <- factor(ggdf$cluster)



# ----------------------------------------------------------------------
### Plot of tsne - all cells, all clusters

## one plot 
ggp <- ggplot(ggdf,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
  geom_point(size = 0.5) +
  labs(x = "t-SNE1", y="t-SNE2")+ 
  theme_bw() +
  theme(axis.text = element_text(size = 12), 
    axis.title  = element_text(size = 15),
    strip.text = element_text(size = 15, hjust = 0),
    strip.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(size = 0.3, linetype = "solid", color = "black")) +
  scale_color_manual(values = colors_tsne[levels(ggdf$cluster)]) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

pdf(file.path(outdir, paste0(prefix, "tSNEone", suffix, ".pdf")), width = 9, height = 7)                 
print(ggp)
dev.off()


## facet per group
ggp <- ggplot(ggdf,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
  geom_point(size=0.5) +
  facet_wrap(~ group) +
  labs(x = "t-SNE1", y="t-SNE2")+
  theme_bw() +
  theme(axis.text = element_text(size = 12), 
    axis.title  = element_text(size = 15),
    strip.text = element_text(size = 15, hjust = 0), 
    strip.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(size = 0.3, linetype = "solid", color = "black")) +
  scale_color_manual(values = colors_tsne[levels(ggdf$cluster)]) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

pdf(file.path(outdir, paste0(prefix, "tSNEgroup", suffix, ".pdf")), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()


## facet per sample
ggp <- ggplot(ggdf,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
  geom_point(size=0.5) +
  facet_wrap(~ sample) +
  labs(x = "t-SNE1", y="t-SNE2")+
  theme_bw() +
  theme(axis.text = element_text(size = 12), 
    axis.title  = element_text(size = 15),
    strip.text = element_text(size = 15, hjust = 0), 
    strip.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(size = 0.3, linetype = "solid", color = "black")) +
  scale_color_manual(values = colors_tsne[levels(ggdf$cluster)]) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

pdf(file.path(outdir, paste0(prefix, "tSNEsample", suffix, ".pdf")), width = 20, height = 20)
print(ggp)
dev.off()









sessionInfo()




