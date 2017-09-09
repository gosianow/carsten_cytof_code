

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

## CD8
# prefix='23CD8allall_29CD8allall_02CD8v2_cl49_clustering_data23CD8_'
# outdir='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/cytokine_profiles_summary'
# path_pvs='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/02v2_CD8/090_cytokine_bimatrix_frequencies_clustering/3responses_both/23CD8allall_29CD8allall_02CD8v2_cl49_frequencies_pvs_glmer_binomial_interglht_top10.xls'
# path_profiles_prefix='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/02v2_CD8/090_cytokine_bimatrix_frequencies_clustering/cytokine_profiles/23CD8allall_29CD8allall_02CD8v2_cl49_clustering_data23CD8_'


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


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------

nr_ccgs <- nlevels(mexpr$ccg)
nr_in_one_row <- 5
nrow <- ceiling(nr_ccgs/nr_in_one_row)
hh <- nrow * 2.5
if(nr_ccgs < nr_in_one_row){
  ww <- nr_ccgs * 2.5 + 1
}else{
  ww <- 2.5 * nr_in_one_row + 1
}

ggp <- ggplot(mexpr, aes(x = label, y = frequencies)) +
  geom_col() +
  ylab("Frequencies") + 
  xlab("") + 
  theme_bw() +
  theme(axis.text = element_text(size = 14, face = "bold", color = "black"), 
    axis.title  = element_text(size = 22, face = "bold"),
    strip.text = element_text(size = 20, hjust = 0), 
    strip.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(size = 0.3, linetype = "solid", color = "black")) +
  facet_wrap(~ ccg, scales = "free", ncol = nr_in_one_row)

pdf(file.path(outdir, paste0(prefix, "cyt_prof_freqs.pdf")), width = ww, height = hh)
print(ggp)
dev.off()























sessionInfo()




