

Sys.time()

# Load packages
library(gdata)
library(ggplot2)
library(reshape2)
library(limma) # for strsplit2
library(RColorBrewer)
library(plyr) # for rbind.fill

##############################################################################
# Test arguments
##############################################################################

## CD4
prefix='23CD4allall_29CD4allall_02CD4v2_cl49_clustering_data23CD4_'
outdir='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/cytokine_profiles_summary'
path_pvs='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/02v2_CD4/090_cytokine_bimatrix_frequencies_clustering/3responses_both/23CD4allall_29CD4allall_02CD4v2_cl49_frequencies_pvs_glmer_binomial_interglht_top10.xls'
path_profiles_prefix='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/02v2_CD4/090_cytokine_bimatrix_frequencies_clustering/cytokine_profiles/23CD4allall_29CD4allall_02CD4v2_cl49_clustering_data23CD4_'
path_clustering_labels='../carsten_cytof/PD1_project/CK_2016-06-23_02_CD4_merging3/030_heatmaps/23CD4_02CD4v2_pca0_merging1_clustering_labels.xls'

## CD8
# prefix='23CD8allall_29CD8allall_02CD8v2_cl49_clustering_data23CD8_'
# outdir='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/cytokine_profiles_summary'
# path_pvs='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/02v2_CD8/090_cytokine_bimatrix_frequencies_clustering/3responses_both/23CD8allall_29CD8allall_02CD8v2_cl49_frequencies_pvs_glmer_binomial_interglht_top10.xls'
# path_profiles_prefix='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/02v2_CD8/090_cytokine_bimatrix_frequencies_clustering/cytokine_profiles/23CD8allall_29CD8allall_02CD8v2_cl49_clustering_data23CD8_'
# path_clustering_labels='../carsten_cytof/PD1_project/CK_2016-06-23_02_CD8_merging3/030_heatmaps/23CD8_02CD8v2_pca0_merging1_clustering_labels.xls'


FDR_cutoff='10'

##############################################################################
# Read in the arguments
##############################################################################

# rm(list = ls())
# 
# args <- (commandArgs(trailingOnly = TRUE))
# for (i in 1:length(args)) {
#   eval(parse(text = args[[i]]))
# }
# 
# cat(paste0(args, collapse = "\n"), fill = TRUE)


##############################################################################

if(!file.exists(outdir)) 
  dir.create(outdir, recursive = TRUE)

suffix <- paste0("_top", FDR_cutoff)
FDR_cutoff <- as.numeric(paste0("0.", FDR_cutoff))
FDR_cutoff


# -----------------------------------------------------------------------------
# Prepare the ggadf and gglabels objects
# -----------------------------------------------------------------------------

comparison <- "adjp_NRvsR_tx"


pvs <- read.table(path_pvs, header = TRUE, sep = "\t", as.is = TRUE)

pvs <- pvs[order(pvs[, gsub("adjp_", "pval_", comparison)]), , drop = FALSE]


which_ccg <- pvs[pvs[, comparison] < FDR_cutoff & !is.na(pvs[, comparison]), "label"]



# -----------------------------------------------------------------------------
# Load cluster profiles for the significant CCGs
# -----------------------------------------------------------------------------


mexpr <- lapply(1:length(which_ccg), function(i){
  # i = 1
  mexpr <- read.table(paste0(path_profiles_prefix, "cl", which_ccg[i], "_cluster_median_expression_all_raw.xls"), header = TRUE, sep = "\t", as.is = TRUE)
  
  mexpr$ccg <- which_ccg[i]
  
  return(mexpr)
})

mexpr <- rbind.fill(mexpr)

mexpr$ccg <- factor(mexpr$ccg, levels = which_ccg)
mexpr$ccg2 <- factor(paste0("CCG ", mexpr$ccg), levels = paste0("CCG ", which_ccg))
mexpr$frequencies <- mexpr$frequencies * 100

# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------

nr_ccgs <- nlevels(mexpr$ccg)
nr_in_one_row <- 7
nrow <- ceiling(nr_ccgs/nr_in_one_row)
hh <- nrow * 2.5
ncol <- ceiling(nr_ccgs/nrow)
ww <- ncol * 2 + 1

ggp <- ggplot(mexpr, aes(x = label, y = frequencies)) +
  geom_col() +
  ylab("Frequency (%)") + 
  xlab("") + 
  theme_bw() +
  theme(axis.text = element_text(size = 12, face = "bold", color = "black"), 
    axis.title  = element_text(size = 15, face = "bold"),
    strip.text = element_text(size = 15, face = "bold", hjust = 0), 
    strip.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(size = 0.3, linetype = "solid", color = "black")) +
  facet_wrap(~ ccg2, scales = "free", ncol = nr_in_one_row)

pdf(file.path(outdir, paste0(prefix, "cyt_prof_freqs.pdf")), width = ww, height = hh)
print(ggp)
dev.off()



# ------------------------------------------------------------
# Load clustering data
# ------------------------------------------------------------

# clustering labels
labels <- read.table(path_clustering_labels, header = TRUE, sep = "\t", as.is = TRUE)

labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))
rownames(labels) <- labels$cluster
labels


# --------------------
# Colors for clusters
# --------------------

# ggplot palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=60 , c=100)[1:n]
}

# color blind palette

colors_muted <- c("#DC050C", "#E8601C", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6", "#F1932D", "#F6C141", "#F7EE55", "#4EB265", "#90C987", "#CAEDAB")
color_ramp <- c(colors_muted, gg_color_hue(max(1, nlevels(labels$label) - length(colors_muted))))

colors_clusters <- color_ramp[1:nlevels(labels$label)]
names(colors_clusters) <- levels(labels$label)
colors_clusters




ggp <- ggplot(mexpr, aes(x = ccg, y = frequencies, fill = label)) +
  geom_col() +
  ylab("Frequency (%)") + 
  xlab("CCG") + 
  theme_bw() +
  theme(axis.text = element_text(size = 15, face = "bold", color = "black"), 
    axis.title  = element_text(size = 15, face = "bold"),
    strip.text = element_text(size = 15, face = "bold", hjust = 0), 
    strip.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(size = 0.3, linetype = "solid", color = "black")) +
  scale_fill_manual(values = colors_clusters)

pdf(file.path(outdir, paste0(prefix, "cyt_prof_freqs2.pdf")), width = nr_ccgs/3 + 2, height = 4)
print(ggp)
dev.off()










sessionInfo()




