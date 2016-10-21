# source("/Users/gosia/Dropbox/UZH/carsten_cytof_code/10_patient_data.R")

# BioC 3.3
# Created 21 Oct 2016
# Updated 21 Oct 2016

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
out_name <- ""

path_fun_models <- '/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_models.R'

# ----------------------------------------------------------
# Read in the data
# ----------------------------------------------------------


data_orig <- read.xls(file.path(rwd, "ck_orig_files", "Patient_data_cytof.xlsx"), stringsAsFactors=FALSE)

data_orig <- data_orig[order(factor(data_orig$TP), factor(data_orig$response)), ]
data_orig$shortname <- paste0("sample", 1:nrow(data_orig))

md <- data.frame(shortname = data_orig$shortname, patient_id = data_orig$ID, day = data_orig$TP, response = data_orig$response, stringsAsFactors = FALSE)

md$patient_id <- factor(md$patient_id)
md$day <- factor(md$day, levels = c("baseline", "TP1"), labels = c("base", "tx"))
md$response <- factor(md$response, levels = c(0, 1), labels = c("NR", "R"))
md$group <- interaction(md$day, md$response, lex.order = TRUE, sep = "_")

md$shortname <- paste0(md$day, "_", md$response, 1:nrow(md))
data_orig$shortname <- md$shortname


vars_cont <- c("leukocytes.count.", "LDH", "ANC", "ALC", "ANC.ALC", "S100", "hemoglobin", "hematocrit", "erytro", "MCV", "MCH", "MCHC", "RDW", "thrombocytes", "monocytes", "eosino", "baso", "IG.abs", "natrium", "kalium", "urea", "creatinin", "eGFR", "bilirubin", "protein", "albumin", "AST", "ALT", "GGT", "alk.phosph", "CRP", "TSH") 


data_cont <- data_orig[, vars_cont, drop = FALSE]
data_cont <- data.frame(t(data_cont))
colnames(data_cont) <- data_orig$shortname
data_cont$feature <- rownames(data_cont)


color_groups <- c(base_NR = "#CC79A7", base_R = "#009E73", tx_NR = "#CC79A7", tx_R = "#009E73")



# ----------------------------------------------------------
# Plot features stratified by day and response
# ----------------------------------------------------------

ggdf <- melt(data_cont, id.vars = c("feature"), value.name = "expr", variable.name = "samp")

## add group info
mm <- match(ggdf$samp, md$shortname)
ggdf$day <- factor(md$day[mm])
ggdf$response <- factor(md$response[mm])
ggdf$group <- interaction(ggdf$day, ggdf$response, lex.order = TRUE, sep = "_")
ggdf$feature <- factor(ggdf$feature, levels = vars_cont)

## calculate mean and sd for the error bars on the plot
ggds <- ddply(ggdf, .(feature, group), summarise, mean = mean(expr, na.rm = TRUE), sd = sd(expr, na.rm = TRUE))


# ------------------------------------

ggp <- ggplot(ggdf, aes(x = group, y = expr, color = group)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(size=2.5, shape = 17, width = 0.5, height = 0) +
  theme_bw() +
  ylab("") +
  xlab("") +
  facet_wrap(~ feature, scale = "free") +
  theme(axis.text.x = element_text(size=10, face="bold"), 
    axis.title.y = element_text(size=10, face="bold"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"), 
    axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
    legend.position = "none") +
  scale_color_manual(values = color_groups)

pdf(file.path(outdir, paste0("plot_features.pdf")), w=18, h=16, onefile=TRUE)
print(ggp)
dev.off()


# ----------------------------------------------------------
# Test for differences between NR and R
# ----------------------------------------------------------


source(path_fun_models)

levels(md$day)
levels(md$response)


if(identical(levels(md$day), c("base", "tx")) && identical(levels(md$response), c("NR", "R"))){
  ## create formulas
  formula_lm <- y ~ response + day + response:day
  formula_lmer <- y ~ response + day + response:day + (1|patient_id)
  
  data_tmp <- data.frame(y = 1:nrow(md), md)
  mm <- model.matrix(formula_lm, data = data_tmp)
  
  ## create contrasts
  contrast_names <- c("NRvsR", "NRvsR_base", "NRvsR_tx", "NRvsR_basevstx")
  k1 <- c(0, 1, 0, 0) # NR vs R in base
  k2 <- c(0, 1, 0, 1) # NR vs R in tx
  k0 <- (k1 + k2) / 2 # NR vs R
  k3 <- c(0, 0, 0, 1) # whether NR vs R is different in base and tx
  K <- matrix(c(k0, k1, k2, k3), nrow = 4, byrow = TRUE)
  rownames(K) <- contrast_names
  
  ### p-value for sorting the output
  pval_name <- "pval_NRvsR"
  ### p-value for plotting the pheatmap2
  adjpval_name <- "adjp_NRvsR"
  ### p-value for plotting the pheatmap3
  adjpval_name_list <- c("adjp_NRvsR", "adjp_NRvsR_base", "adjp_NRvsR_tx", "adjp_NRvsR_basevstx")
  
}else{
  stop("Metadata does not fit to any the models that are specified !!!")  
}


exprc <- data_cont
data_features <- c("feature")

models2fit <- c("lm_interglht", "rlm_interglht", "lmer_interglht", "test_wilcoxon")


for(k in models2fit){
  # k = "lm_interglht"
  print(k)
  
  ### p-value for sorting the output
  pval_name <- "pval_NRvsR"
  ### p-value for plotting the pheatmap2
  adjpval_name <- "adjp_NRvsR"
  ### p-value for plotting the pheatmap3
  adjpval_name_list <- c("adjp_NRvsR", "adjp_NRvsR_base", "adjp_NRvsR_tx", "adjp_NRvsR_basevstx")
  
  
  switch(k,
    lm_interglht = {
      
      # Fit a LM with interactions + test contrasts with multcomp pckg
      fit_out <- fit_lm_interglht(data = exprc, md, method = "lm", formula = formula_lm, K = K)
      
    }, 
    lmer_interglht = {
      
      # Fit a lmer with interactions + test contrasts with multcomp pckg
      fit_out <- fit_lmer_interglht(data = exprc, md, formula = formula_lmer, K = K)
      
    },
    lmrob_interglht = {
      ## Problems with running lmrob!!!
      fit_out <- fit_lm_interglht(data = exprc, md, method = "lmrob", formula = formula_lm, K = K)
    },
    rlm_interglht = {
      fit_out <- fit_lm_interglht(data = exprc, md, method = "rlm", formula = formula_lm, K = K)
    },
    test_wilcoxon = {
      fit_out <- test_wilcoxon(data = exprc, md)
      
      ### p-value for sorting the output
      pval_name <- "pval_NRvsR_base"
      ### p-value for plotting the pheatmap2
      adjpval_name <- "adjp_NRvsR_base"
      ### p-value for plotting the pheatmap3
      adjpval_name_list <- c("adjp_NRvsR_base", "adjp_NRvsR_tx")
      
      
    }
    
  )
  
  
  # ----------------------------------------
  # Extract p-values and coeffs
  # ----------------------------------------
  
  pvs <- data.frame(exprc[, data_features, drop = FALSE], fit_out[["pvals"]])
  coeffs <- data.frame(exprc[, data_features, drop = FALSE], fit_out[["coeffs"]])
  
  oo <- order(pvs[, pval_name], decreasing = FALSE)
  pvs <- pvs[oo, , drop = FALSE]
  coeffs <- coeffs[oo, , drop = FALSE]
  
  ## save the results
  write.table(pvs, file = file.path(outdir, paste0(prefix, "expr_", out_name, "_pvs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  write.table(coeffs, file = file.path(outdir, paste0(prefix, "expr_", out_name, "_coeffs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  
  
  # ----------------------------------------
  # Plot a heatmap of significant cases
  # ----------------------------------------
  
  ### normalize the expression
  expr_norm <- exprc[, c(data_features, md[md$response != "HD", "shortname"])]
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
  expr_all <- merge(pvs, expr_norm, by = data_features, all.x = TRUE, sort = FALSE)
  
  
  
  # -----------------------------
  ### Plot one heatmap with base and tx
  
  ## group the expression by cluster
  if(is.null(adjpval_name)){
    which_top_pvs <- FALSE
  }else{
    expr_all <- expr_all[order(expr_all[, adjpval_name]), , drop = FALSE]
    
    which_top_pvs <- expr_all[, adjpval_name] < 0.05 & !is.na(expr_all[, adjpval_name])
    which(which_top_pvs)
  }
  
  
  if(sum(which_top_pvs) > 0) {
    print("Plot pheatmap2")
    
    expr_heat <- expr_all[which_top_pvs, , drop = FALSE]
    
    # -----------------------------
    ## order the samples by base and tx
    samples2plot <- md[md$response %in% c("NR", "R"), ]
    samples2plot <- samples2plot[order(samples2plot$day, samples2plot$response), ]
    
    ## gap in the heatmap 
    gaps_col <- c(max(grep("base_NR", samples2plot$group)), rep(max(grep("base", samples2plot$day)), 2), max(grep("tx_NR", samples2plot$group)))
    gaps_row <- NULL
    
    ## expression 
    expr <- expr_heat[ , samples2plot$shortname, drop = FALSE]
    
    labels_row <- paste0(expr_all[which_top_pvs, data_features], " (", sprintf( "%.02e", expr_heat[, adjpval_name]), ")") 
    labels_col <- colnames(expr)
    
    pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "expr_", out_name, "_", k, "_pheatmap2", suffix, ".pdf")))
    
    
  }
  
  
  # -----------------------------
  ### Plot heatmaps with R vs NR for base and tx
  
  for(i in levels(md$day)){
    # i = "tx"
    
    adjpval_name <- paste0("adjp_NRvsR_", i)
    
    ## group the expression by cluster
    expr_all <- expr_all[order(expr_all[, adjpval_name]), , drop = FALSE]
    
    which_top_pvs <- expr_all[, adjpval_name] < 0.05 & !is.na(expr_all[, adjpval_name])
    which(which_top_pvs)
    
    if(sum(which_top_pvs) > 0){
      print(paste0("Plot pheatmap_", i))
      
      expr_heat <- expr_all[which_top_pvs, , drop = FALSE]
      
      # -----------------------------
      ## order the samples by NR and R
      
      samples2plot <- md[md$response %in% c("NR", "R"), ]
      samples2plot <- samples2plot[order(samples2plot$response, samples2plot$day), ]
      samples2plot <- samples2plot[grep(i, samples2plot$day), ]
      
      ## gap in the heatmap 
      gaps_col <- sum(grepl("NR", samples2plot$response))
      gaps_row <- NULL
      
      ## expression
      expr <- expr_heat[ , samples2plot$shortname, drop = FALSE]
      
      labels_row <- paste0(expr_all[which_top_pvs, data_features], " (", sprintf( "%.02e", expr_heat[, adjpval_name]), ")") 
      labels_col <- colnames(expr)
      
      pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "expr_", out_name, "_", k, "_pheatmap_", i, suffix, ".pdf")))
      
      
    }
    
  }
  
  # -----------------------------
  ### Plot one heatmap with R vs NR + heatmap with p-values for NRvsR_base, NRvsR_tx and NRvsR_basevstx
  
  ## group the expression by cluster and order by adjpval
  for(i in length(adjpval_name_list):1){
    expr_all <- expr_all[order(expr_all[, adjpval_name_list[i]]), , drop = FALSE]
  }
  
  # which_top_pvs <- rowSums(expr_all[, adjpval_name_list, drop = FALSE] < 0.05, na.rm = TRUE) > 0 & rowSums(is.na(expr_all[, adjpval_name_list, drop = FALSE])) < length(adjpval_name_list)
  which_top_pvs <- rep(TRUE, nrow(expr_all))
  which(which_top_pvs)
  
  if(sum(which_top_pvs) > 0) {
    print("Plot pheatmap3")
    
    expr_heat <- expr_all[which_top_pvs, , drop = FALSE]
    
    # -----------------------------
    ## order the samples by base and tx
    samples2plot <- md[md$response %in% c("NR", "R"), ]
    samples2plot <- samples2plot[order(samples2plot$day, samples2plot$response), ]
    
    ## gap in the heatmap 
    gaps_col <- c(max(grep("base_NR", samples2plot$group)), rep(max(grep("base", samples2plot$day)), 2), max(grep("tx_NR", samples2plot$group)))
    gaps_row <- NULL
    
    ## expression 
    expr <- expr_heat[ , samples2plot$shortname, drop = FALSE]
    
    labels_row <- paste0(expr_all[which_top_pvs, data_features]) 
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

























