# BioC 3.3
# Created 21 Oct 2016
# Updated 09 Dec 2016

##############################################################################
Sys.time()
##############################################################################

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

rwd <- "/Users/gosia/Dropbox/UZH/carsten_cytof/Patientdata/"
outdir <- "ck_analysis"

path_data <- "ck_orig_files/Patient_data_cytof_FACS_09_12_16_CK_GN.xlsx"
path_variables <- "ck_orig_files/variables_to_use.xlsx"

path_fun_models <- '/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_models.R'
path_fun_formulas <- '/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_formulas_1dataset_2responses_base.R'
path_fun_plot_heatmaps <- "/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_plot_heatmaps_for_sign_expr.R"



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

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

prefix <- ""
suffix <- ""
out_name <- "patientdata"


color_groups <- c(base_NR = "#CC79A7", base_R = "#009E73", tx_NR = "#CC79A7", tx_R = "#009E73")


# ----------------------------------------------------------
# Read in the data
# ----------------------------------------------------------

data_orig <- read.xls(path_data, as.is = TRUE, check.names = FALSE)

variables_orig <- read.xls(path_variables, as.is = TRUE, check.names = FALSE)


### Create metadata

md <- data.frame(patient_id = data_orig[, "ID"], day = data_orig[, "TP"], response = data_orig[, "RESPONSE"], shortname = data_orig[, "shortname"], stringsAsFactors = FALSE)

md$patient_id <- factor(md$patient_id)
md$day <- factor(md$day, levels = c("baseline", "TP1"), labels = c("base", "tx"))
md$response <- factor(md$response, levels = c(0, 1), labels = c("NR", "R"))
md$group <- interaction(md$day, md$response, lex.order = TRUE, sep = "_")
md$shortname <- paste0(md$day, "_", md$shortname)



###############################################################################
### Analysis for continous variables
###############################################################################


# ----------------------------------------------------------
### DE analysis for the following variables
# ----------------------------------------------------------

vars_cont <- variables_orig$variable[variables_orig$usage == 1 & variables_orig$continous == 1]

data_cont <- data_orig[, vars_cont, drop = FALSE]
data_cont <- apply(data_cont, 2, as.numeric)

data_cont <- data.frame(t(data_cont))
colnames(data_cont) <- md$shortname
data_cont$marker <- rownames(data_cont)
data_cont$cluster <- ""
data_cont$label <- ""
rownames(data_cont) <- NULL


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

pdf(file.path(outdir, paste0(prefix, "expr_", out_name, "_plot.pdf")), w=18, h=16, onefile=TRUE)
print(ggp)
dev.off()





# ----------------------------------------------------------
### normalize the expression
# ----------------------------------------------------------

exprc <- data_cont

expr_norm <- exprc[, c("cluster", "label", "marker", md[md$response != "HD", "shortname"])]
th <- 2.5

days <- levels(md$day)

### Normalized to mean = 0 and sd = 1 per day
for(i in days){
  # i = "tx"
  
  cols2use <- md[md$response != "HD" & md$day == i, "shortname"]
  
  expr_norm[, cols2use] <- t(apply(expr_norm[, cols2use, drop = FALSE], 1, function(x){ 
    
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




# ----------------------------------------------------------
# Test for differences between NR and R
# ----------------------------------------------------------

source(path_fun_models)
source(path_fun_formulas)
source(path_fun_plot_heatmaps)


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
  
  ### add p-value info
  expr_all <- merge(pvs, expr_norm, by = c("cluster", "label", "marker"), all.x = TRUE, sort = FALSE)
  
  prefix2 <- paste0(out_name, "_", k, "_")
  
  plot_heatmaps_for_sign_expr(expr_all = expr_all, md = md, FDR_cutoff = 0.05, pval_name2 = pval_name2, adjpval_name2 = adjpval_name2, pval_name_list = pval_name_list, adjpval_name_list = adjpval_name_list, breaks = breaks, legend_breaks = legend_breaks, outdir = outdir, prefix = prefix, prefix2 = prefix2, suffix = suffix)
  
  
}


###############################################################################
### Analysis for binary variables
###############################################################################


# ----------------------------------------------------------
### DE analysis for the following variables
# ----------------------------------------------------------

vars_bi <- variables_orig$variable[variables_orig$usage == 1 & variables_orig$continous == 0]


data_bi <- data_orig[, vars_bi, drop = FALSE]

data_bi <- apply(data_bi, 2, as.numeric)

data_bi <- data.frame(t(data_bi))
colnames(data_bi) <- md$shortname

data_bi$marker <- rownames(data_bi)
data_bi$cluster <- ""
data_bi$label <- ""
rownames(data_bi) <- NULL



# ----------------------------------------------------------
# Plot features stratified by day and response
# ----------------------------------------------------------

ggdf <- melt(data_bi, id.vars = c("cluster", "label", "marker"), value.name = "expr", variable.name = "sample")

## add group info
mm <- match(ggdf$sample, md$shortname)
ggdf$day <- factor(md$day[mm])
ggdf$response <- factor(md$response[mm])
ggdf$group <- interaction(ggdf$day, ggdf$response, lex.order = TRUE, sep = "_")
ggdf$marker <- factor(ggdf$marker, levels = vars_bi)

ggdf$expr <- factor(ggdf$expr, levels = c(1, 0))

ggdf <- ggdf[complete.cases(ggdf), ]

# ------------------------------------

ggp <- ggplot(ggdf, aes(x = group, color = group, fill = expr)) +
  geom_bar() +
  theme_bw() +
  ylab("") +
  xlab("") +
  facet_wrap(~ marker, scale = "fixed") +
  theme(axis.text.x = element_text(size=10, face="bold"), 
    axis.title.y = element_text(size=10, face="bold"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"), 
    axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
    legend.position = "right",
    strip.text = element_text(size = 10, hjust = 0), strip.background = element_blank()) +
  scale_color_manual(values = color_groups) +
  scale_fill_manual(values = c("gray50", "gray80"))

pdf(file.path(outdir, paste0(prefix, "frequencies_", out_name, "_plot.pdf")), w=8, h=3, onefile=TRUE)
print(ggp)
dev.off()



# ----------------------------------------------------------
### normalize the expression
# ----------------------------------------------------------

exprc <- data_bi

expr_norm <- exprc[, c("cluster", "label", "marker", md[md$response != "HD", "shortname"])]
th <- 2.5

days <- levels(md$day)

### Normalized to mean = 0 and sd = 1 per day
for(i in days){
  # i = "tx"
  
  cols2use <- md[md$response != "HD" & md$day == i, "shortname"]
  
  expr_norm[, cols2use] <- t(apply(expr_norm[, cols2use, drop = FALSE], 1, function(x){ 
    
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




# ----------------------------------------------------------
# Test for differences between NR and R
# ----------------------------------------------------------

source(path_fun_models)
source(path_fun_formulas)
source(path_fun_plot_heatmaps)


models2fit <- c("glmer_binomial_interglht_01")

freq_out <- data_bi


for(k in models2fit){
  # k = "glmer_binomial_interglht_01"
  print(k)
  
  switch(k,
    
    glmer_binomial_interglht_01 = {
      # Fit a GLMM binomial with interactions + test contrasts with multcomp pckg
      fit_out <- fit_glmer_interglht_01(data = freq_out, md, family = "binomial", formula = formula_glmer_binomial_01, K = K, skippNAs = FALSE)
      
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
  write.table(pvs, file = file.path(outdir, paste0(prefix, "frequencies_", out_name, "_pvs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  write.table(coeffs, file = file.path(outdir, paste0(prefix, "frequencies_", out_name, "_coeffs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  
  
  # ----------------------------------------
  # Plot a heatmap of significant cases
  # ----------------------------------------
  
  ### add p-value info
  expr_all <- merge(pvs, expr_norm, by = c("cluster", "label", "marker"), all.x = TRUE, sort = FALSE)
  
  prefix2 <- paste0(out_name, "_", k, "_")
  
  plot_heatmaps_for_sign_expr(expr_all = expr_all, md = md, FDR_cutoff = 0.05, pval_name2 = pval_name2, adjpval_name2 = adjpval_name2, pval_name_list = pval_name_list, adjpval_name_list = adjpval_name_list, breaks = breaks, legend_breaks = legend_breaks, outdir = outdir, prefix = prefix, prefix2 = prefix2, suffix = suffix)
  
  
}





sessionInfo()











