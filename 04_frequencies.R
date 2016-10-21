##############################################################################
## <<04_frequencies.R>>

# BioC 3.3
# Created 27 July 2016
# Updated 25 Aug 2016

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

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
# freq_prefix='23_01_pca1_merging6_'
# freq_outdir='050_frequencies'
# path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.xlsx'
# path_clustering='030_heatmaps/23_01_pca1_merging6_clustering.xls'
# path_clustering_labels='030_heatmaps/23_01_pca1_merging6_clustering_labels.xls'
# path_fun_models='/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_models.R'

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

setwd(rwd)

prefix <- freq_prefix
suffix <- ""
outdir <- freq_outdir

if(!file.exists(outdir)) 
  dir.create(outdir)

source(path_fun_models)

# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

md <- read.xls(path_metadata, stringsAsFactors=FALSE)

# add more info about samples
cond_split <- strsplit2(md$condition, "_")
colnames(cond_split) <- c("day", "response")

md[, c("day", "response")] <- cond_split
md$response <- factor(md$response, levels = c("NR", "R", "HD"))
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

color_samples <- md$color
names(color_samples) <- md$shortname


# ------------------------------------------------------------
# Load cluster data
# ------------------------------------------------------------

## clustering labels
labels <- read.table(path_clustering_labels, header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))

## clustering
clustering <- read.table(path_clustering, header = TRUE, sep = "\t", as.is = TRUE)

## drop the "drop" cluster
if("drop" %in% labels$label){
  
  clust2drop <- labels$cluster[labels$label == "drop"]
  cells2drop <- clustering$cluster != clust2drop
  clustering <- clustering[cells2drop, , drop = FALSE]
  labels <- labels[labels$label != "drop", ,drop = FALSE]
  labels$label <- factor(labels$label)
  
}

clust <- clustering[, "cluster"]

# ---------------------------------------
# Calculate the cluster frequencies per sample
# ---------------------------------------

samp <- clustering[, "sample_id"]

# calculate frequencies
freq <- table(clust, samp)

prop <- t(t(freq) / colSums(freq)) * 100 # proportion of clusters in samples

# use labels as names of clusters
mlab <- match(rownames(freq), labels$cluster)


### Save the frequencies and proportions
prop_out <- data.frame(labels[mlab, c("cluster", "label")], as.data.frame.matrix(prop))

freq_out <- data.frame(labels[mlab, c("cluster", "label")], as.data.frame.matrix(freq))

write.table(prop_out, file=file.path(outdir, paste0(prefix, "frequencies.xls")), row.names=FALSE, quote=FALSE, sep="\t")
write.table(freq_out, file=file.path(outdir, paste0(prefix, "counts.xls")), row.names=FALSE, quote=FALSE, sep="\t")



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

## replace _ with \n
levels(ggdf$group) <- gsub("_", "\n", levels(ggdf$group))

## calculate mean and sd for the error bars on the plot
ggds <- ddply(ggdf, .(group, cluster), summarise, mean = mean(prop), sd = sd(prop))

# add more info about samples
ggdf$day <- strsplit2(ggdf$group, "\n")[, 1]
ggds$day <- strsplit2(ggds$group, "\n")[, 1]


ggdf$day <- factor(ggdf$day)
ggds$day <- factor(ggds$day)


# ------------------------------------
### plot each cluster as a separate page in the pdf file
ggp <- list()
clusters <- levels(ggdf$cluster)

for(i in 1:nlevels(ggdf$cluster)){
  # i = 1
  
  df <- ggdf[ggdf$cluster == clusters[i], , drop = FALSE]
  ds <- ggds[ggds$cluster == clusters[i], , drop = FALSE]
  
  ggp[[i]] <- ggplot(df, aes(x = group, y = prop, color = group)) +
    geom_jitter(size=2.5, shape = 17, width = 0.5, height = 0) +
    geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean, ymax=mean), color='black', width=0.4) +
    geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean-sd, ymax=mean+sd), color='black', width=0.25) +
    ggtitle(clusters[i]) +
    theme_bw() +
    ylab("Frequency") +
    xlab("") +
    theme(axis.text.x = element_text(size=12, face="bold"), 
      axis.title.y = element_text(size=12, face="bold"), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(), 
      axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"), 
      axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
      legend.position = "none") +
    scale_color_manual(values = color_groups)
  
}

pdf(file.path(outdir, paste0(prefix, "frequencies_plot.pdf")), w=5, h=4, onefile=TRUE)
for(i in seq(length(ggp)))
  print(ggp[[i]])
dev.off()


# ------------------------------------
# plot all clusters in one pdf; colors per group; points; separate pdf for base and tx

days <- levels(ggdf$day)

for(i in 1:nlevels(ggdf$day)){
  # i = 1
  
  df <- ggdf[ggdf$day == days[i], , drop = FALSE]
  ds <- ggds[ggds$day == days[i], , drop = FALSE]
  
  ggp <- ggplot(df, aes(x = cluster, y = prop, color = group)) +
    geom_point(size=2, shape = 16, alpha = 0.8, position = position_jitterdodge(jitter.width = 3, jitter.height = 0, dodge.width = 0.7)) +
    geom_errorbar(data=ds, aes(x=cluster, y=mean, ymin=mean, ymax=mean), width=0.4, position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.7), size = 1) +
    geom_errorbar(data=ds, aes(x=cluster, y=mean, ymin=mean-sd, ymax=mean+sd), width=0.25, position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.7), size = 1) +
    theme_bw() +
    ylab("Frequency") +
    xlab("") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12, face="bold"), 
      axis.title.y = element_text(size=12, face="bold"), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(), 
      axis.line.x = element_line(size = 0.5, linetype = "solid", color = "black"), 
      axis.line.y = element_line(size = 0.5, linetype = "solid", color = "black"),
      legend.title = element_blank(), legend.position = "right", legend.key = element_blank()) +
    guides(color = guide_legend(ncol = 1)) +
    scale_color_manual(values = color_groups)
  
  pdf(file.path(outdir, paste0(prefix, "frequencies_plot_", days[i] ,".pdf")), w=7, h=4)
  print(ggp)
  dev.off()
  
}


# ------------------------------------------------------------
# Test for frequency differences between groups
# ------------------------------------------------------------

source(path_fun_models)

levels(md$day)
levels(md$response)


if(identical(levels(md$day), c("base", "tx")) && identical(levels(md$response), c("NR", "R", "HD"))){
  ## create formulas
  formula_lm <- y ~ response + day + response:day
  formula_lmer <- y ~ response + day + response:day + (1|patient_id)
  
  formula_glm_binomial <- cbind(y, total-y) ~ response + day + response:day
  formula_glm_beta <- y/total ~ response + day + response:day
  formula_glmer_binomial <- y/total ~ response + day + response:day + (1|patient_id)
  
  ## create contrasts
  contrast_names <- c("NRvsR", "NRvsR_base", "NRvsR_tx", "NRvsR_basevstx")
  k0 <- c(0, 1, 0, 0, 1/2, 0) # NR vs R
  k1 <- c(0, 1, 0, 0, 0, 0) # NR vs R in base
  k2 <- c(0, 1, 0, 0, 1, 0) # NR vs R in tx
  k3 <- c(0, 0, 0, 0, 1, 0) # whether NR vs R is different in base and tx
  K <- matrix(c(k0, k1, k2, k3), nrow = 4, byrow = TRUE)
  rownames(K) <- contrast_names
  
  ### p-value for sorting the output
  pval_name <- "pval_NRvsR"
  ### p-value for plotting the pheatmap2
  adjpval_name <- "adjp_NRvsR"
  ### p-value for plotting the pheatmap3
  adjpval_name_list <- c("adjp_NRvsR", "adjp_NRvsR_base", "adjp_NRvsR_tx", "adjp_NRvsR_basevstx")
  
  
}else if(identical(levels(md$day), "base") && identical(levels(md$response), c("NR", "R", "HD"))){
  ## create formulas
  formula_lm <- y ~ response
  formula_lmer <- y ~ response + (1|patient_id)
  
  formula_glm_binomial <- cbind(y, total-y) ~ response
  formula_glm_beta <- y/total ~ response
  formula_glmer_binomial <- y/total ~ response + (1|patient_id)
  
  ## create contrasts
  contrast_names <- c("NRvsR_base")
  k1 <- c(0, 1, 0) # NR vs R in base
  K <- matrix(k1, nrow = 1, byrow = TRUE)
  rownames(K) <- contrast_names
  
  ### p-value for sorting the output
  pval_name <- "pval_NRvsR_base"
  ### p-value for plotting the pheatmap2
  adjpval_name <- NULL
  ### p-value for plotting the pheatmap3
  adjpval_name_list <- "adjp_NRvsR_base"
  
}else{
  stop("Metadata does not fit to any the models that are specified !!!")  
}


models2fit <- c("glm_binomial_interglht", "glm_quasibinomial_interglht", "glmer_binomial_interglht", "test_wilcoxon", "glmmadmb_fixed_betabinomial_interglht", "glmmadmb_fixed_beta_interglht", "lmer_arcsinesqrt_interglht", "lmer_logit_interglht")


for(k in models2fit){
  # k = "glmer_logit_interglht"
  print(k)
  
  switch(k,
    glm_binomial_interglht = {
      # Fit a GLM binomial with interactions + test contrasts with multcomp pckg
      fit_out <- fit_glm_interglht(data = freq_out, md, family = "binomial", formula = formula_glm_binomial, K = K)
      
    }, 
    glm_quasibinomial_interglht = {
      # Fit a GLM quasibinomial with interactions + test contrasts with multcomp pckg
      fit_out <- fit_glm_interglht(data = freq_out, md, family = "quasibinomial", formula = formula_glm_binomial, K = K)
      
    },
    glmer_binomial_interglht = {
      # Fit a GLMM binomial with interactions + test contrasts with multcomp pckg
      fit_out <- fit_glmer_interglht(data = freq_out, md, family = "binomial", formula = formula_glmer_binomial, K = K)
      
    },
    lmer_logit_interglht = {
      
      logit_freq_out <- freq_out
      logit_freq_out[md$shortname] <- logit(t(t(freq_out[md$shortname]) / colSums(freq_out[md$shortname])))
      ## Be carefull about Inf and -Inf for prop = 0, 1
      fit_out <- fit_lmer_interglht(data = logit_freq_out, md, formula = formula_lmer, K = K)
      
    },
    lmer_arcsinesqrt_interglht = {
      
      ass_freq_out <- freq_out
      ass_freq_out[md$shortname] <- asin(sqrt((t(t(freq_out[md$shortname]) / colSums(freq_out[md$shortname])))))
      
      fit_out <- fit_lmer_interglht(data = ass_freq_out, md, formula = formula_lmer, K = K)
      
    },
    glmmadmb_fixed_beta_interglht = {
      fit_out <- fit_glm_interglht(data = freq_out, md, family = "beta", formula = formula_glm_beta, K = K)
    },
    glmmadmb_fixed_betabinomial_interglht = {
      fit_out <- fit_glm_interglht(data = freq_out, md, family = "betabinomial", formula = formula_glm_binomial, K = K)
    },
    test_wilcoxon = {
      
      n01_freq_out <- freq_out
      n01_freq_out[md$shortname] <- t(t(freq_out[md$shortname]) / colSums(freq_out[md$shortname]))
      
      fit_out <- test_wilcoxon(data = n01_freq_out, md)
    }
    
  )
  
  # ----------------------------------------
  # Extract p-values and coeffs
  # ----------------------------------------
  
  pvs <- data.frame(freq_out[, c("cluster", "label")], fit_out[["pvals"]])
  coeffs <- data.frame(freq_out[, c("cluster", "label")], fit_out[["coeffs"]])
  
  oo <- order(pvs[, pval_name], decreasing = FALSE)
  pvs <- pvs[oo, , drop = FALSE]
  coeffs <- coeffs[oo, , drop = FALSE]
  
  ## save the results
  write.table(pvs, file=file.path(outdir, paste0(prefix, "frequencies_pvs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  write.table(coeffs, file=file.path(outdir, paste0(prefix, "frequencies_coeffs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  
  # table(pvs$adjp_NRvsR < 0.05, useNA = "always")
  # table(pvs$adjp_NRvsR_base < 0.05, useNA = "always")
  # table(pvs$adjp_NRvsR_tx < 0.05, useNA = "always")
  # table(pvs$adjp_NRvsR_basevstx < 0.05, useNA = "always")
  # 
  # table(pvs$pval_NRvsR < 0.05, useNA = "always")
  # table(pvs$pval_NRvsR_base < 0.05, useNA = "always")
  # table(pvs$pval_NRvsR_tx < 0.05, useNA = "always")
  # table(pvs$pval_NRvsR_basevstx < 0.05, useNA = "always")
  
  
  # ----------------------------------------
  # Plot a heatmap of significant cases - transform proportions with arcsin-sqrt so the dispersion is the same for low and high props.
  # ----------------------------------------
  
  ### normalize the expression
  ass_freq_out <- freq_out
  ass_freq_out[md$shortname] <- asin(sqrt((t(t(freq_out[md$shortname]) / colSums(freq_out[md$shortname])))))
  
  expr_norm <- ass_freq_out[, c("cluster", "label", md[md$response != "HD", "shortname"])]
  th <- 2.5
  
  days <- levels(md$day)
  
  ### Normalized to mean = 0 and sd = 1 per day
  for(i in days){
    # i = "base"
    expr_norm[, md[md$response != "HD" & md$day == i, "shortname"]] <- t(apply(expr_norm[, md[md$response != "HD" & md$day == i, "shortname"], drop = FALSE], 1, function(x){ x <- (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE); x[x > th] <- th; x[x < -th] <- -th; return(x)}))
  }
  
  breaks = seq(from = -th, to = th, length.out = 101)
  legend_breaks = seq(from = -round(th), to = round(th), by = 1)
  
  ### add p-value info
  expr_all <- merge(pvs, expr_norm, by = c("cluster", "label"), all.x = TRUE, sort = FALSE)
  
  
  # -----------------------------
  ### Plot one heatmap with R vs NR
  
  ## group the expression by cluster
  if(is.null(adjpval_name)){
    which_top_pvs <- FALSE
  }else{
    expr_all <- expr_all[order(expr_all[, adjpval_name]), , drop = FALSE]
    
    which_top_pvs <- expr_all[, adjpval_name] < 0.05 & !is.na(expr_all[, adjpval_name])
    which(which_top_pvs)
  }
  
  
  if(sum(which_top_pvs) > 0) {
    
    expr_heat <- expr_all[which_top_pvs, , drop = FALSE]
    
    # -----------------------------
    ## order the samples by base and tx
    samples2plot <- md[md$response %in% c("NR", "R"), ]
    samples2plot <- samples2plot$shortname[order(samples2plot$day, samples2plot$response)]
    
    ## gap in the heatmap 
    gaps_col <- c(max(grep("base_NR", samples2plot)), rep(max(grep("base", samples2plot)), 2), max(grep("tx_NR", samples2plot)))
    gaps_row <- NULL
    
    ## expression scaled by row
    expr <- expr_heat[ , samples2plot, drop = FALSE]
    
    labels_row <- paste0(expr_heat$label, " (", sprintf( "%.02e", expr_heat[, adjpval_name]), ")") 
    labels_col <- colnames(expr)
    
    pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "frequencies_", k, "_pheatmap2", suffix, ".pdf")))
    
    
  }
  
  
  # -----------------------------
  ### Plot two heatmaps with R vs NR for base and tx
  
  for(i in levels(md$day)){
    # i = "base"
    
    adjpval_name <- paste0("adjp_NRvsR_", i)
    
    ## group the expression by cluster
    expr_all <- expr_all[order(expr_all[, adjpval_name]), , drop = FALSE]
    
    which_top_pvs <- expr_all[, adjpval_name] < 0.05 & !is.na(expr_all[, adjpval_name])
    which(which_top_pvs)
    
    if(sum(which_top_pvs) > 0){
      
      expr_heat <- expr_all[which_top_pvs, , drop = FALSE]
      
      # -----------------------------
      ## order the samples by NR and R
      
      samples2plot <- md[md$response %in% c("NR", "R"), ]
      samples2plot <- samples2plot$shortname[order(samples2plot$response, samples2plot$day)]
      samples2plot <- samples2plot[grep(i, samples2plot)]
      
      ## gap in the heatmap 
      gaps_col <- sum(grepl("_NR", samples2plot))
      gaps_row <- NULL
      
      ## expression scaled by row
      expr <- expr_heat[ , samples2plot, drop = FALSE]
      
      labels_row <- paste0(expr_heat$label, " (", sprintf( "%.02e", expr_heat[, adjpval_name]), ")") 
      labels_col <- colnames(expr)
      
      pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "frequencies_", k, "_pheatmap_", i, suffix, ".pdf")))
      
      
    }
    
  }
  
  # -----------------------------
  ### Plot one heatmap with R vs NR + heatmap with p-values for NRvsR_base, NRvsR_tx and NRvsR_basevstx
  
  ## group the expression by cluster and order by adjpval
  for(i in length(adjpval_name_list):1){
    expr_all <- expr_all[order(expr_all[, adjpval_name_list[i]]), , drop = FALSE]
  }
  
  
  # which_top_pvs <- rowSums(expr_all[, adjpval_name_list] < 0.05, na.rm = TRUE) > 0 & rowSums(is.na(expr_all[, adjpval_name_list])) == 0
  # which_top_pvs <- rowSums(is.na(expr_all[, adjpval_name_list])) < length(adjpval_name_list)
  which_top_pvs <- rep(TRUE, nrow(expr_all))
  which(which_top_pvs)
  
  if(sum(which_top_pvs) > 0) {
    
    expr_heat <- expr_all[which_top_pvs, , drop = FALSE]
    
    # -----------------------------
    ## order the samples by base and tx
    samples2plot <- md[md$response %in% c("NR", "R"), ]
    samples2plot <- samples2plot$shortname[order(samples2plot$day, samples2plot$response)]
    
    ## gap in the heatmap 
    gaps_col <- c(max(grep("base_NR", samples2plot)), rep(max(grep("base", samples2plot)), 2), max(grep("tx_NR", samples2plot)))
    gaps_row <- NULL
    
    ## expression 
    expr <- expr_heat[ , samples2plot, drop = FALSE]
    
    labels_row <- paste0(expr_heat$label) 
    labels_col <- colnames(expr)
    
    pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "frequencies_", k, "_pheatmap3", suffix, ".pdf")))
    
    
    pvs_heat <- expr_heat[, adjpval_name_list, drop = FALSE]
    
    labels_col <- colnames(pvs_heat)
    gaps_col <- NULL
    
    
    pheatmap(pvs_heat, cellwidth = 60, cellheight = 24, color = c("grey50", "grey70", "grey90"), breaks = c(0, 0.05, 0.1, 1), legend_breaks = c(0, 0.05, 0.1, 1), legend = FALSE, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, display_numbers = TRUE, number_format = "%.02e", number_color = "black", fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "frequencies_", k, "_pheatmap3pvs", suffix, ".pdf")))
    
    
  }
  
  
  
  # ----------------------------------------
  # Plot coefficients NRvsR for base and tx (to show that they correlate)
  # ----------------------------------------
  
  
  if("adjp_NRvsR_basevstx" %in% colnames(pvs)){
    
    adjpval_name <- "adjp_NRvsR_basevstx"
    
    ggdf <- coeffs[, c("NRvsR_base", "NRvsR_tx")]
    
    limmin <- min(ggdf, na.rm = TRUE)
    limmax <- max(ggdf, na.rm = TRUE)
    
    ggdf$interaction <- factor(pvs[, adjpval_name] < 0.05, levels = c("FALSE", "TRUE"))
    
    
    ggp <- ggplot(data = ggdf, aes(x = NRvsR_base, y = NRvsR_tx, shape = interaction)) +
      geom_point(size = 3, alpha = 0.9) +
      geom_abline(intercept = 0, slope = 1) +
      coord_cartesian(xlim = c(limmin, limmax), ylim = c(limmin, limmax)) +
      theme_bw() +
      theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face="bold"))
    
    pdf(file.path(outdir, paste0(prefix, "frequencies_", k, "_coeffs", suffix, ".pdf")), w=6, h=5, onefile=TRUE)
    print(ggp)
    dev.off()
    
  }
  
  
  
}






sessionInfo()













################################
### 04_frequencies done!
################################