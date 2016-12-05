# source("/Users/gosia/Dropbox/UZH/carsten_cytof_code/10_patient_data.R")

# BioC 3.3
# Created 21 Oct 2016
# Updated 05 Dec 2016

library(gdata)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2)
library(limma) # for strsplit2
library(RColorBrewer)
library(pheatmap)
library(gtools) # for logit



rwd <- "/Users/gosia/Dropbox/UZH/carsten_cytof/Patientdata/"

setwd(rwd)

outdir <- "ck_analysis"
dir.create(outdir, recursive = TRUE)

prefix <- ""
suffix <- ""
out_name <- "patientdata"

path_fun_models <- '/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_models.R'
path_fun_formulas <- '/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_formulas_1dataset_2responses.R'

path_data <- "Patient_data_cytof_FACS_05_12_16.xlsx"

# ----------------------------------------------------------
# Read in the data
# ----------------------------------------------------------

data_orig <- read.xls(file.path(rwd, "ck_orig_files", path_data), as.is = TRUE, check.names = FALSE)

### Create metadata

md <- data.frame(patient_id = data_orig[, "patient ID"], day = data_orig[, "time point"], response = data_orig[, "response (RECIST, +/-)"], shortname = data_orig[, "shortname"], stringsAsFactors = FALSE)

md$patient_id <- factor(md$patient_id)
md$day <- factor(md$day, levels = c("before", "after"), labels = c("base", "tx"))
md$response <- factor(md$response, levels = c(0, 1), labels = c("NR", "R"))
md$group <- interaction(md$day, md$response, lex.order = TRUE, sep = "_")
md$shortname <- paste0(md$day, "_", md$shortname)


### DE analysis for the following variables

vars_cont <- colnames(data_orig[, c(17:20, 22, 23, 25:34, 38:46, 48:51)])


data_cont <- data_orig[, vars_cont, drop = FALSE]
data_cont <- apply(data_cont, 2, as.numeric)

data_cont <- data.frame(t(data_cont))
colnames(data_cont) <- md$shortname
data_cont$marker <- rownames(data_cont)
data_cont$cluster <- ""
data_cont$label <- ""
rownames(data_cont) <- NULL

color_groups <- c(base_NR = "#CC79A7", base_R = "#009E73", tx_NR = "#CC79A7", tx_R = "#009E73")


# ----------------------------------------------------------
# Plot features stratified by day and response
# ----------------------------------------------------------

ggdf <- melt(data_cont, id.vars = c("cluster", "label", "marker"), value.name = "expr", variable.name = "sample")

## add group info
mm <- match(ggdf$sample, md$shortname)
ggdf$day <- factor(md$day[mm])
ggdf$response <- factor(md$response[mm])
ggdf$group <- interaction(ggdf$day, ggdf$response, lex.order = TRUE, sep = "_")
ggdf$marker <- factor(ggdf$marker, levels = vars_cont)



# ------------------------------------

ggp <- ggplot(ggdf, aes(x = group, y = expr, color = group)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(size=2.5, shape = 17, width = 0.5, height = 0) +
  theme_bw() +
  ylab("") +
  xlab("") +
  facet_wrap(~ marker, scale = "free") +
  theme(axis.text.x = element_text(size=10, face="bold"), 
    axis.title.y = element_text(size=10, face="bold"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"), 
    axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
    legend.position = "none",
    strip.text = element_text(size = 10, hjust = 0), strip.background = element_blank()) +
  scale_color_manual(values = color_groups)

pdf(file.path(outdir, paste0("expr_", out_name, "_plot.pdf")), w=18, h=16, onefile=TRUE)
print(ggp)
dev.off()


# ----------------------------------------------------------
# Test for differences between NR and R
# ----------------------------------------------------------

source(path_fun_models)
source(path_fun_formulas)


exprc <- data_cont

models2fit <- c("lmer_interglht")


for(k in models2fit){
  # k = "lmer_interglht"
  print(k)
  
  switch(k,
    lm_interglht = {
      
      # Fit a LM with interactions + test contrasts with multcomp pckg
      fit_out <- fit_lm_interglht(data = exprc, md, method = "lm", formula = formula_lm, K = K, skippNAs = FALSE)
      
    }, 
    lmer_interglht = {
      
      # Fit a lmer with interactions + test contrasts with multcomp pckg
      fit_out <- fit_lmer_interglht(data = exprc, md, formula = formula_lmer, K = K, skippNAs = FALSE)
      
    },
    lmrob_interglht = {
      ## Problems with running lmrob!!!
      fit_out <- fit_lm_interglht(data = exprc, md, method = "lmrob", formula = formula_lm, K = K, skippNAs = FALSE)
    },
    rlm_interglht = {
      fit_out <- fit_lm_interglht(data = exprc, md, method = "rlm", formula = formula_lm, K = K, skippNAs = FALSE)
    }
    
  )
  
  
  # ----------------------------------------
  # Extract p-values and coeffs
  # ----------------------------------------
  
  pvs <- data.frame(exprc[, c("cluster", "label", "marker")], fit_out[["pvals"]])
  coeffs <- data.frame(exprc[, c("cluster", "label", "marker")], fit_out[["coeffs"]])
  
  oo <- order(pvs[, pval_name1], decreasing = FALSE)
  pvs <- pvs[oo, , drop = FALSE]
  coeffs <- coeffs[oo, , drop = FALSE]
  
  ## save the results
  write.table(pvs, file = file.path(outdir, paste0(prefix, "expr_", out_name, "_pvs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  write.table(coeffs, file = file.path(outdir, paste0(prefix, "expr_", out_name, "_coeffs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  
  
  # ----------------------------------------
  # Plot a heatmap of significant cases
  # ----------------------------------------
  
  ### normalize the expression
  expr_norm <- exprc[, c("cluster", "label", "marker", md[md$response != "HD", "shortname"])]
  th <- 2.5
  
  days <- levels(md$day)
  
  ### Normalized to mean = 0 and sd = 1 per day
  for(i in days){
    # i = "base"
    expr_norm[, md[md$response != "HD" & md$day == i, "shortname"]] <- t(apply(expr_norm[, md[md$response != "HD" & md$day == i, "shortname"], drop = FALSE], 1, function(x){ 
      
      if(sum(!is.na(x)) == 0)
        return(x)
      
      if(sum(!is.na(x)) < 2)
        return(x-mean(x, na.rm = TRUE))
      
      sdx <- sd(x, na.rm = TRUE)
      if(sdx == 0)
        x <- (x-mean(x, na.rm = TRUE))
      else 
        x <- (x-mean(x, na.rm = TRUE))/sdx
      
      x[x > th] <- th
      x[x < -th] <- -th
      
      return(x)}))
  }
  
  breaks = seq(from = -th, to = th, length.out = 101)
  legend_breaks = seq(from = -round(th), to = round(th), by = 1)
  
  ### add p-value info
  expr_all <- merge(pvs, expr_norm, by = c("cluster", "label", "marker"), all.x = TRUE, sort = FALSE)
  
  
  
  # -----------------------------
  ### Plot one heatmap with base and tx
  
  ## group the expression by cluster
  if(is.null(adjpval_name2)){
    which_top_pvs <- FALSE
  }else{
    adjpval_name <- adjpval_name2
    pval_name <- pval_name2
    expr_all <- expr_all[order(expr_all$label, expr_all[, pval_name]), , drop = FALSE]
    
    which_top_pvs <- expr_all[, adjpval_name] < 0.05 & !is.na(expr_all[, adjpval_name])
    which(which_top_pvs)
  }
  
  
  if(sum(which_top_pvs) > 0) {
    print("Plot pheatmap2")
    
    expr_heat <- expr_all[which_top_pvs, , drop = FALSE]
    
    # -----------------------------
    ## order the samples by base and tx
    samples2plot <- md[md$response %in% c("NR", "R"), ]
    samples2plot <- samples2plot$shortname[order(samples2plot$day, samples2plot$response)]
    
    ## gap in the heatmap 
    gaps_col <- c(max(grep("base_NR", samples2plot)), rep(max(grep("base", samples2plot)), 2), max(grep("tx_NR", samples2plot)))
    gaps_row <- unique(cumsum(table(expr_heat$label)))
    gaps_row <- gaps_row[gaps_row > 0]
    if(length(gaps_row) == 1) 
      gaps_row <- NULL
    
    ## expression 
    expr <- expr_heat[ , samples2plot, drop = FALSE]
    
    labels_row <- paste0(expr_heat$label, "/ ", expr_heat$marker, " (", sprintf( "%.02e", expr_heat[, adjpval_name]), ")") 
    labels_col <- colnames(expr)
    
    pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "expr_", out_name, "_", k, "_pheatmap2", suffix, ".pdf")))
    
    
  }
  
  
  # -----------------------------
  ### Plot heatmaps with R vs NR for base and tx
  
  for(i in levels(md$day)){
    # i = "tx"
    
    adjpval_name <- paste0("adjp_NRvsR_", i)
    pval_name <- paste0("pval_NRvsR_", i)
    
    ## group the expression by cluster
    expr_all <- expr_all[order(expr_all$label, expr_all[, pval_name]), , drop = FALSE]
    
    which_top_pvs <- expr_all[, adjpval_name] < 0.05 & !is.na(expr_all[, adjpval_name])
    which(which_top_pvs)
    
    if(sum(which_top_pvs) > 0){
      print(paste0("Plot pheatmap_", i))
      
      expr_heat <- expr_all[which_top_pvs, , drop = FALSE]
      
      # -----------------------------
      ## order the samples by NR and R
      
      samples2plot <- md[md$response %in% c("NR", "R"), ]
      samples2plot <- samples2plot$shortname[order(samples2plot$response, samples2plot$day)]
      samples2plot <- samples2plot[grep(i, samples2plot)]
      
      ## gap in the heatmap 
      gaps_col <- sum(grepl("_NR", samples2plot))
      gaps_row <- unique(cumsum(table(expr_heat$label)))
      gaps_row <- gaps_row[gaps_row > 0]
      if(length(gaps_row) == 1) 
        gaps_row <- NULL
      
      ## expression
      expr <- expr_heat[ , samples2plot, drop = FALSE]
      
      labels_row <- paste0(expr_heat$label, "/ ", expr_heat$marker, " (", sprintf( "%.02e", expr_heat[, adjpval_name]), ")") 
      labels_col <- colnames(expr)
      
      pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "expr_", out_name, "_", k, "_pheatmap_", i, suffix, ".pdf")))
      
      
    }
    
  }
  
  # -----------------------------
  ### Plot one heatmap with R vs NR + heatmap with p-values for NRvsR_base, NRvsR_tx and NRvsR_basevstx
  
  ## group the expression by cluster and order by adjpval
  for(i in length(pval_name_list):1){
    expr_all <- expr_all[order(expr_all[, pval_name_list[i]]), , drop = FALSE]
  }
  
  expr_all <- expr_all[order(expr_all$label), , drop = FALSE]
  
  # which_top_pvs <- rowSums(expr_all[, adjpval_name_list, drop = FALSE] < 0.05, na.rm = TRUE) > 0 & rowSums(is.na(expr_all[, adjpval_name_list, drop = FALSE])) < length(adjpval_name_list)
  which_top_pvs <- rep(TRUE, nrow(expr_all))
  which(which_top_pvs)
  
  if(sum(which_top_pvs) > 0) {
    print("Plot pheatmap3")
    
    expr_heat <- expr_all[which_top_pvs, , drop = FALSE]
    
    # -----------------------------
    ## order the samples by base and tx
    samples2plot <- md[md$response %in% c("NR", "R"), ]
    samples2plot <- samples2plot$shortname[order(samples2plot$day, samples2plot$response)]
    
    ## gap in the heatmap 
    gaps_col <- c(max(grep("base_NR", samples2plot)), rep(max(grep("base", samples2plot)), 2), max(grep("tx_NR", samples2plot)))
    gaps_row <- unique(cumsum(table(expr_heat$label)))
    gaps_row <- gaps_row[gaps_row > 0]
    if(length(gaps_row) == 1) 
      gaps_row <- NULL
    
    ## expression 
    expr <- expr_heat[ , samples2plot, drop = FALSE]
    
    labels_row <- paste0(expr_heat$label, "/ ", expr_heat$marker) 
    labels_col <- colnames(expr)
    
    pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "expr_", out_name, "_", k, "_pheatmap3", suffix, ".pdf")))
    
    pvs_heat <- expr_heat[, adjpval_name_list, drop = FALSE]
    
    labels_col <- colnames(pvs_heat)
    gaps_col <- NULL
    
    pheatmap(pvs_heat, cellwidth = 60, cellheight = 24, color = c("grey50", "grey70", "grey90"), breaks = c(0, 0.05, 0.1, 1), legend_breaks = c(0, 0.05, 0.1, 1), legend = FALSE, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, display_numbers = TRUE, number_format = "%.02e", number_color = "black", fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "expr_", out_name, "_", k, "_pheatmap3pvs", suffix, ".pdf")))
    
    
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
      geom_point(size = 3, alpha = 0.75) +
      geom_abline(intercept = 0, slope = 1) +
      coord_cartesian(xlim = c(limmin, limmax), ylim = c(limmin, limmax)) +
      theme_bw() +
      theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face="bold"))
    
    pdf(file.path(outdir, paste0(prefix, "expr_", out_name, "_", k, "_coeffs", suffix, ".pdf")), w=6, h=5, onefile=TRUE)
    print(ggp)
    dev.off()
    
    
  }
  
  
  
  
} # models2fit

























