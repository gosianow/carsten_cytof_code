##############################################################################
## <<10_facs_validation.R>>

# BioC 3.3
# Created 9 Dec 2016

##############################################################################
Sys.time()
##############################################################################

# Load packages
library(gdata)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2)
library(limma) # for strsplit2
library(RColorBrewer)
library(pheatmap)
library(gtools) # for logit


##############################################################################
# Test arguments
##############################################################################

prefix='facs_valid4_'
outdir='../carsten_cytof/PD1_project/FACSvalidation4/ck_analysis'

path_metadata='../carsten_cytof/PD1_project/FACSvalidation4/ck_orig_files/CK_2017-01-27_ValidationDataset2_metadata.xlsx'
path_freqs='../carsten_cytof/PD1_project/FACSvalidation4/ck_orig_files/CK_2017-01-27_ValidationDataset2_frequencies.xlsx'

path_fun_models='00_models.R'
path_fun_formulas='00_formulas_1dataset_3responses_base.R'
path_fun_plot_heatmaps='00_plot_heatmaps_for_sign_freqs.R'

pdf_hight=4

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

suffix <- ""

if(!file.exists(outdir)) 
  dir.create(outdir, recursive = TRUE)



# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

md <- read.xls(path_metadata, stringsAsFactors=FALSE)

# add more info about samples
cond_split <- strsplit2(md$condition, "_")
colnames(cond_split) <- c("day", "response")

md[, c("day", "response")] <- cond_split
md$response <- factor(md$response, levels = c("HD", "NR", "R"))
md$response <- factor(md$response)
md$day <- factor(md$day, levels = c("base", "tx"))
md$day <- factor(md$day)
md$patient_id <- factor(md$patient_id)


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

# ---------------------------------------
# Load frequencies
# ---------------------------------------


prop_out <- read.xls(path_freqs, stringsAsFactors=FALSE)

labels <- prop_out[, c("cluster", "label")]

# ------------------------------------------------------------
### Plot frequencies
# ------------------------------------------------------------

ggdf <- melt(prop_out, id.vars = c("cluster", "label"), value.name = "prop", variable.name = "samp")

## use labels as clusters
ggdf$cluster <- factor(ggdf$label, levels = labels$label)
ggdf <- ggdf[, c("cluster", "samp", "prop")]

## add group info
mm <- match(ggdf$samp, md$shortname)
ggdf$group <- factor(md$condition[mm])
ggdf$response <- factor(md$response[mm])
ggdf$data <- "FACS"

## replace _ with \n
levels(ggdf$group) <- gsub("_", "\n", levels(ggdf$group))

## calculate mean and sd for the error bars on the plot
ggds <- ddply(ggdf, .(group, cluster), summarise, mean = mean(prop), sd = sd(prop))

# add more info about samples
ggdf$day <- strsplit2(ggdf$group, "\n")[, 1]
ggdf$day <- factor(ggdf$day)



# ------------------------------------
# plot all clusters in one pdf; colors per group; boxplots + points

nr_cluster <- nlevels(ggdf$cluster)
nr_in_one_row <- 7
nrow <- ceiling(nr_cluster/nr_in_one_row)
hh <- nrow * 2.5
ncol <- ceiling(nr_cluster/nrow)
ww <- ncol * 2 + 1

ggp <- ggplot() +
  geom_boxplot(data = ggdf, aes(x = group, y = prop, color = response, shape = data, fill = response), width = 0.9, 
    position = position_dodge(width = 1), outlier.colour = NA, alpha = 0.4, show.legend = FALSE) +
  geom_point(data = ggdf, aes(x = group, y = prop, color = response, shape = data), size = 2.5, alpha = 0.9, 
    position = position_jitterdodge(jitter.width = 1, jitter.height = 0, dodge.width = 1)) +
  theme_bw() +
  ylab("Frequency (%)") +
  xlab("") +
  theme(axis.text.x = element_blank(), 
    axis.text.y = element_text(size = 15, face = "bold", color = "black"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 15, face = "bold"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    axis.line.x = element_line(size = 0.5, linetype = "solid", color = "black"), 
    axis.line.y = element_line(size = 0.5, linetype = "solid", color = "black"),
    legend.title = element_blank(), 
    legend.position = "right", 
    legend.key = element_blank(), 
    legend.text = element_text(size = 18, face = "bold"),
    strip.background = element_blank(), 
    strip.text = element_text(size = 18, hjust = 0, face = "bold")) +
  scale_color_manual(values = color_response) +
  scale_fill_manual(values = color_response) +
  facet_wrap(~ cluster, scales = "free", nrow = nrow) +
  guides(shape = guide_legend(override.aes = list(size = 4)), 
    color = guide_legend(override.aes = list(size = 4)))

pdf(file.path(outdir, paste0(prefix, "frequencies_plot.pdf")), w = ww, h = hh)
print(ggp)
dev.off()



# ------------------------------------------------------------
### normalize the expression
# ------------------------------------------------------------
# 
# ass_freq_out <- prop_out
# ass_freq_out[md$shortname] <- asin(sqrt(prop_out[md$shortname] / 100))
# 
# logit_freq_out <- prop_out
# logit_freq_out[md$shortname] <- logit(prop_out[md$shortname] / 100)
# 
# n01_freq_out <- prop_out
# n01_freq_out[md$shortname] <- prop_out[md$shortname] / 100
# 
# 
# expr_norm <- ass_freq_out[, c("cluster", "label", md[md$response != "HD", "shortname"])]
# th <- 2.5
# 
# days <- levels(md$day)
# 
# ### Normalized to mean = 0 and sd = 1 per day
# for(i in days){
#   # i = "base"
#   expr_norm[, md[md$response != "HD" & md$day == i, "shortname"]] <- t(apply(expr_norm[, md[md$response != "HD" & md$day == i, "shortname"], drop = FALSE], 1, function(x){ 
#     
#     if(sum(!is.na(x)) == 0)
#       return(x)
#     
#     if(sum(!is.na(x)) < 2)
#       return(x-mean(x, na.rm = TRUE))
#     
#     sdx <- sd(x, na.rm = TRUE)
#     if(sdx == 0)
#       x <- (x-mean(x, na.rm = TRUE))
#     else 
#       x <- (x-mean(x, na.rm = TRUE))/sdx
#     
#     x[x > th] <- th
#     x[x < -th] <- -th
#     
#     return(x)}))
# }
# 
# breaks = seq(from = -th, to = th, length.out = 101)
# legend_breaks = seq(from = -round(th), to = round(th), by = 1)
# 
# # ------------------------------------------------------------
# # Test for frequency differences between groups
# # ------------------------------------------------------------
# ## The model functions do not anlyse a cluster with NAs; 
# ## For merged data it means such cluster was not present in all the datasets
# ## For expression data clusters with no cells are skipped
# 
# 
# ### Load functions fitting models
# source(path_fun_models)
# ### Load formulas that are fit in the models - this function may change the md object!!!
# source(path_fun_formulas)
# 
# source(path_fun_plot_heatmaps)
# 
# 
# levels(md$day)
# levels(md$response)
# 
# 
# models2fit <- c("lm_arcsinesqrt_interglht", "lm_logit_interglht", "glmmadmb_fixed_beta_interglht")
# 
# for(k in models2fit){
#   # k = "lm_arcsinesqrt_interglht"
#   print(k)
#   
#   switch(k,
#     lm_arcsinesqrt_interglht = {
#       
#       fit_out <- fit_lm_interglht(data = ass_freq_out, md, formula = formula_lm, K = K)
#       
#     },
#     lm_logit_interglht = {
#       
#       ## Be carefull about Inf and -Inf for prop = 0, 1
#       fit_out <- fit_lm_interglht(data = logit_freq_out, md, formula = formula_lm, K = K)
#       
#     },
#     glmmadmb_fixed_beta_interglht = {
#       
#       fit_out <- fit_glm_interglht(data = n01_freq_out, md, family = "beta", formula = formula_glm_beta_01, K = K)
#       
#     }
#     
#     
#   )
#   
#   
#   # ----------------------------------------
#   # Extract p-values and coeffs
#   # ----------------------------------------
#   
#   pvs <- data.frame(ass_freq_out[, c("cluster", "label")], fit_out[["pvals"]])
#   coeffs <- data.frame(ass_freq_out[, c("cluster", "label")], fit_out[["coeffs"]])
#   
#   oo <- order(pvs[, pval_name1], decreasing = FALSE)
#   pvs <- pvs[oo, , drop = FALSE]
#   coeffs <- coeffs[oo, , drop = FALSE]
#   
#   ## save the results
#   write.table(pvs, file=file.path(outdir, paste0(prefix, "frequencies_pvs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
#   write.table(coeffs, file=file.path(outdir, paste0(prefix, "frequencies_coeffs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
#   
#   
#   # ----------------------------------------
#   # Plot a heatmap of significant cases - transform proportions with arcsin-sqrt so the dispersion is the same for low and high props.
#   # ----------------------------------------
#   
#   ### add p-value info
#   expr_all <- merge(pvs, expr_norm, by = c("cluster", "label"), all.x = TRUE, sort = FALSE)
#   
#   prefix2 <- paste0(k, "_")
#   
#   FDR_cutoff <- 0.05
#   
#   plot_heatmaps_for_sign_freqs(expr_all = expr_all, md = md, FDR_cutoff = FDR_cutoff, pval_name2 = pval_name2, adjpval_name2 = adjpval_name2, pval_name_list = pval_name_list, adjpval_name_list = adjpval_name_list, breaks = breaks, legend_breaks = legend_breaks, color_response = color_response, outdir = outdir, prefix = prefix, prefix2 = prefix2, suffix = suffix)
#   
#   
#   
# }
# 
# 




sessionInfo()













################################
### 04_frequencies done!
################################