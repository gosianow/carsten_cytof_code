##############################################################################
## <<08_frequencies_merged.R>>

# BioC 3.3
# Created 22 Sep 2016
# Updated 22 Sep 2016

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
library(tools)

##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-merged_23_29'
freq_prefix='23m6_29m4_'
freq_outdir='050_frequencies'

path_metadata=c('/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.xlsx','/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_29_01.xlsx')
path_counts=c('/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01/050_frequencies/23_01_pca1_merging6_counts.xls','/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-29_01/050_frequencies/29_01_pca1_merging4_counts.xls')
data_name=c('data23','data29')

path_fun_models='/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_models_merged.R'

##############################################################################
# Read in the arguments
##############################################################################

# args <- (commandArgs(trailingOnly = TRUE))
# for (i in 1:length(args)) {
#   eval(parse(text = args[[i]]))
# }
# 
# print(args)

##############################################################################

if(!file.exists(rwd)) 
  dir.create(rwd, recursive = TRUE)
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

md <- lapply(1:length(data_name), function(i){
  
  path <- path_metadata[i]
  md <- read.xls(path, stringsAsFactors=FALSE)
  md$data <- data_name[i]
  md
  
})

md <- rbind.fill(md)


# add more info about samples
cond_split <- strsplit2(md$condition, "_")
colnames(cond_split) <- c("day", "response")

md[, c("day", "response")] <- cond_split
md$response <- factor(md$response, levels = c("NR", "R", "HD"))
md$day <- factor(md$day, levels = c("base", "tx"))
md$patient_id <- factor(md$patient_id)

md$data_day <- interaction(md$data, md$day, lex.order = TRUE)


### Colors 
colors <- unique(md[, c("condition", "color")])
colors$condition <- factor(colors$condition)
## replace _ with \n
levels(colors$condition) <- gsub("_", "\n", levels(colors$condition ))

color_groups <- colors$color
names(color_groups) <- colors$condition

color_samples <- md$color
names(color_samples) <- md$shortname


# ---------------------------------------
# Load cluster frequencies 
# ---------------------------------------

freq <- lapply(1:length(data_name), function(i){
  
  path <- path_counts[i]
  freq <- read.table(path, header = TRUE, sep = "\t", as.is = TRUE)
  freq
  
})

freq_out <- Reduce(function(...) merge(..., by = c("cluster", "label"), all=TRUE, sort = FALSE), freq)

freq_out <- freq_out[freq_out$label != "drop", , drop = FALSE]

if(any(!complete.cases(freq_out)))
  stop(paste0("Files: ", paste(basename(path_counts), collapse = ", "), " have different clusters!!!"))


prop_out <- freq_out
prop_out[md$shortname] <- t(t(freq_out[md$shortname]) / colSums(freq_out[md$shortname]))


labels <- freq_out[, c("cluster", "label")]
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))



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
## add data info
ggdf$data <- factor(md$data[mm])

## replace _ with \n
levels(ggdf$group) <- gsub("_", "\n", levels(ggdf$group))

## calculate mean and sd for the error bars on the plot
ggds <- ddply(ggdf, .(group, cluster, data), summarise, mean = mean(prop), sd = sd(prop))

## add more info about samples
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
  
  ggp[[i]] <- ggplot(df, aes(x = group, y = prop, color = group, shape = data)) +
    geom_point(size=2.5, position = position_jitterdodge(jitter.width = 3, jitter.height = 0, dodge.width = 0.7)) +
    geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean, ymax=mean), color='black', width=0.4, position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.7)) +
    geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean-sd, ymax=mean+sd), color='black', width=0.25, position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.7)) +
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
      legend.position = "right") +
    scale_color_manual(values = color_groups)
  
}

pdf(file.path(outdir, paste0(prefix, "frequencies_plot.pdf")), w=6, h=4, onefile=TRUE)
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
  
  ggp <- ggplot(df, aes(x = cluster, y = prop, color = group, shape = data)) +
    geom_point(size=2, alpha = 0.8, position = position_jitterdodge(jitter.width = 1, jitter.height = 0, dodge.width = 0.7)) +
    geom_errorbar(data=ds, aes(x=cluster, y=mean, ymin=mean, ymax=mean), width=0.4, position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.7)) +
    geom_errorbar(data=ds, aes(x=cluster, y=mean, ymin=mean-sd, ymax=mean+sd), width=0.25, position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.7)) +
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
  
  pdf(file.path(outdir, paste0(prefix, "frequencies_plot_", days[i] ,".pdf")), w = 10, h = 4)
  print(ggp)
  dev.off()
  
}


# ------------------------------------------------------------
# Test for frequency differences between groups
# ------------------------------------------------------------

source(path_fun_models)

models2fit <- c("glm_binomial_interglht", "glm_quasibinomial_interglht", "glmer_binomial_interglht", "glmmadmb_fixed_betabinomial_interglht", "glmmadmb_fixed_beta_interglht")

# x <- seq(0, 1, 0.01)
# plot(x, logit(x), type = "l")
# lines(x, asin(sqrt(x)), col = "blue")
# lines(x, asin(x), col = "red")
# lines(x, sqrt(x), col = "green")

for(k in models2fit){
  # k = "glmer_logit_interglht"
  
  switch(k,
    glm_binomial_interglht = {
      # Fit a GLM binomial with interactions + test contrasts with multcomp pckg
      fit_out <- fit_glm_interglht(data = freq_out, md, family = "binomial")
      
    }, 
    glm_quasibinomial_interglht = {
      # Fit a GLM quasibinomial with interactions + test contrasts with multcomp pckg
      fit_out <- fit_glm_interglht(data = freq_out, md, family = "quasibinomial")
      
    },
    glmer_binomial_interglht = {
      # Fit a GLMM binomial with interactions + test contrasts with multcomp pckg
      fit_out <- fit_glmer_interglht(data = freq_out, md, family = "binomial")
      
    },
    lmer_logit_interglht = {
      
      logit_freq_out <- freq_out
      logit_freq_out[md$shortname] <- logit(t(t(freq_out[md$shortname]) / colSums(freq_out[md$shortname])))
      ## Be carefull about Inf and -Inf for prop = 0, 1
      fit_out <- fit_lmer_interglht(data = logit_freq_out, md)
      
    },
    lmer_arcsinesqrt_interglht = {
      
      ass_freq_out <- freq_out
      ass_freq_out[md$shortname] <- asin(sqrt((t(t(freq_out[md$shortname]) / colSums(freq_out[md$shortname])))))
      
      fit_out <- fit_lmer_interglht(data = ass_freq_out, md)
      
    },
    glmmadmb_fixed_beta_interglht = {
      fit_out <- fit_glm_interglht(data = freq_out, md, family = "beta")
    },
    glmmadmb_fixed_betabinomial_interglht = {
      fit_out <- fit_glm_interglht(data = freq_out, md, family = "betabinomial")
    },
    glmmadmb_mixed_beta_interglht = {
      fit_out <- fit_glmer_interglht(data = freq_out, md, family = "beta")
    },
    glmmadmb_mixed_betabinomial_interglht = {
      fit_out <- fit_glmer_interglht(data = freq_out, md, family = "betabinomial")
    }
    
  )
  
  # ----------------------------------------
  # Extract p-values and coeffs
  # ----------------------------------------
  
  pvs <- data.frame(freq_out[, c("cluster", "label")], fit_out[["pvals"]])
  coeffs <- data.frame(freq_out[, c("cluster", "label")], fit_out[["coeffs"]])
  
  oo <- order(pvs$pval_NRvsR, decreasing = FALSE)
  pvs <- pvs[oo, ]
  coeffs <- coeffs[oo, ]
  
  ## save the results
  write.table(pvs, file=file.path(outdir, paste0(prefix, "frequencies_pvs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  write.table(coeffs, file=file.path(outdir, paste0(prefix, "frequencies_coeffs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  
  table(pvs$adjp_NRvsR < 0.05, useNA = "always")
  table(pvs$adjp_NRvsR_base < 0.05, useNA = "always")
  table(pvs$adjp_NRvsR_tx < 0.05, useNA = "always")
  table(pvs$adjp_NRvsR_basevstx < 0.05, useNA = "always")
  
  table(pvs$pval_NRvsR < 0.05, useNA = "always")
  table(pvs$pval_NRvsR_base < 0.05, useNA = "always")
  table(pvs$pval_NRvsR_tx < 0.05, useNA = "always")
  table(pvs$pval_NRvsR_basevstx < 0.05, useNA = "always")
  
  
  # ----------------------------------------
  # Plot a heatmap of significant cases
  # ----------------------------------------
  
  ### normalize the expression for base and tx separately
  ass_freq_out <- freq_out
  ass_freq_out[md$shortname] <- asin(sqrt((t(t(freq_out[md$shortname]) / colSums(freq_out[md$shortname])))))
  
  expr_norm <- ass_freq_out[, grep("cluster|label|_NR|_R", colnames(ass_freq_out))]
  th <- 2.5
  
  for(i in c("base", "tx")){
    # i = "base"
    expr_norm[, grep(i, colnames(expr_norm))] <- t(apply(expr_norm[, grep(i, colnames(expr_norm)), drop = FALSE], 1, function(x){ x <- (x-mean(x))/sd(x); x[x > th] <- th; x[x < -th] <- -th; return(x)}))
  }
  
  breaks = seq(from = -th, to = th, length.out = 101)
  legend_breaks = seq(from = -round(th), to = round(th), by = 1)
  
  ### add p-value info
  expr_all <- merge(pvs, expr_norm, by = c("cluster", "label"), all.x = TRUE, sort = FALSE)
  
  
  # -----------------------------
  ### Plot one heatmap with R vs NR
  
  adjpval_name <- "adjp_NRvsR"
  
  ## group the expression by cluster
  expr_all <- expr_all[order(expr_all[, adjpval_name]), , drop = FALSE]
  
  which_top_pvs <- expr_all[, adjpval_name] < 0.05 & !is.na(expr_all[, adjpval_name])
  which(which_top_pvs)
  
  if(sum(which_top_pvs) > 0) {
    
    expr_heat <- expr_all[which_top_pvs, , drop = FALSE]
    
    # -----------------------------
    ## order the samples by NR and R
    
    samples2plot <- md[md$response %in% c("NR", "R"), ]
    samples2plot <- samples2plot$shortname[order(samples2plot$response, samples2plot$day)]
    
    ## gap in the heatmap 
    gaps_col <- sum(grepl("_NR", samples2plot))
    gaps_row <- NULL
    
    ## expression scaled by row
    expr <- expr_heat[ , samples2plot, drop = FALSE]
    
    labels_row <- paste0(expr_heat$label, " (", sprintf( "%.02e", expr_heat[, adjpval_name]), ")") 
    labels_col <- colnames(expr)
    
    pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "frequencies_", k, "_pheatmap1", suffix, ".pdf")))
    
    
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
  
  for(i in c("base", "tx")){
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
  
  adjpval_name <- c("adjp_NRvsR", "adjp_NRvsR_base", "adjp_NRvsR_tx", "adjp_NRvsR_basevstx")
  
  ## group the expression by cluster
  expr_all <- expr_all[order(expr_all[, adjpval_name[2]], expr_all[, adjpval_name[3]], expr_all[, adjpval_name[4]], expr_all$label), , drop = FALSE]
  
  # which_top_pvs <- rowSums(expr_all[, adjpval_name] < 0.05, na.rm = TRUE) > 0 & rowSums(is.na(expr_all[, adjpval_name])) == 0
  # which_top_pvs <- rowSums(is.na(expr_all[, adjpval_name])) < length(adjpval_name)
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
    
    
    pvs_heat <- expr_heat[, adjpval_name, drop = FALSE]
    
    labels_col <- colnames(pvs_heat)
    gaps_col <- NULL
    
    
    pheatmap(pvs_heat, cellwidth = 60, cellheight = 24, color = c("grey50", "grey70", "grey90"), breaks = c(0, 0.05, 0.1, 1), legend_breaks = c(0, 0.05, 0.1, 1), legend = FALSE, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, display_numbers = TRUE, number_format = "%.02e", number_color = "black", fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "frequencies_", k, "_pheatmap3pvs", suffix, ".pdf")))
    
    
  }
  
  
  
  # ----------------------------------------
  # Plot coefficients NRvsR for base and tx (to show that they correlate)
  # ----------------------------------------
  
  adjpval_name <- "adjp_NRvsR_basevstx"
  
  ggdf <- coeffs[, c("NRvsR_base", "NRvsR_tx")]
  
  if(sum(complete.cases(ggdf)) > 0){
    
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