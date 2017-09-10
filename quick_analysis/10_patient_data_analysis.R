
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

outdir <- "../carsten_cytof/PD1_project/PD-1_PatientClinicalData/ck_analysis3"

path_data <- "../carsten_cytof/PD1_project/PD-1_PatientClinicalData/ck_orig_files3/2017_09_06_updated_info_CyTOF_patients.xlsx"
path_variables <- "../carsten_cytof/PD1_project/PD-1_PatientClinicalData/ck_orig_files3/variables_to_use3.xlsx"

### Analysis between NA and R from base 
path_fun_models <- '00_models.R'
path_fun_formulas <- '00_formulas_1dataset_2responses_base.R'
path_fun_plot_heatmaps <- "00_plot_heatmaps_for_sign_expr.R"



##############################################################################
# Read in the arguments
##############################################################################

# rm(list = ls())

# args <- (commandArgs(trailingOnly = TRUE))
# for (i in 1:length(args)) {
#   eval(parse(text = args[[i]]))
# }

# cat(paste0(args, collapse = "\n"), fill = TRUE)


##############################################################################

FDR_cutoff='05'
FDR_cutoff <- as.numeric(paste0("0.", FDR_cutoff))
FDR_cutoff


dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


suffix <- ""
out_name <- "patientdata"


color_groups <- c(base_NR = "#CC79A7", base_R = "#009E73", tx_NR = "#CC79A7", tx_R = "#009E73")
color_response <- c(NR = "#CC79A7", R = "#009E73")


###############################################################################
###
### Analysis of:
### CyTOF samples (baseline)
### FACS samples (baseline)
### CyTOF + FACS samples (baseline)
###
###############################################################################


# ----------------------------------------------------------
# Read in the data
# ----------------------------------------------------------

analysis_of <- "CyTOFandFACS"


if(analysis_of == "CyTOF"){
  
  prefix <- paste0(analysis_of, "_")
  
  data_orig <- read.xls(path_data, as.is = TRUE, check.names = FALSE)
  
  variables_orig <- read.xls(path_variables, as.is = TRUE, check.names = FALSE)
  
  data_orig <- data_orig[data_orig$use == analysis_of & data_orig$TP == "baseline", , drop = FALSE]
  
  
  ### Create metadata
  
  md <- data.frame(patient_id = data_orig[, "ID"], day = data_orig[, "TP"], response = data_orig[, "RESPONSE"], shortname = data_orig[, "shortname"], data = data_orig[, "batches"], stringsAsFactors = FALSE)
  
  md$patient_id <- factor(md$patient_id)
  md$day <- factor(md$day, levels = c("baseline", "TP1"), labels = c("base", "tx"))
  md$data <- factor(md$data, levels = c("1", "2"), labels = c("data23", "data29"))
  md$response <- factor(md$response, levels = c(0, 1), labels = c("NR", "R"))
  
  md$group <- interaction(md$day, md$response, lex.order = TRUE, sep = "_")
  md$shortname <- paste0(md$day, "_", md$shortname)
  md$data_day <- interaction(md$data, md$day, lex.order = TRUE, drop = TRUE)
  
  ### More adjustments
  # We use only the base samples
  md$day <- factor(md$day)
  
}


if(analysis_of == "FACS"){
  
  prefix <- paste0(analysis_of, "_")
  
  data_orig <- read.xls(path_data, as.is = TRUE, check.names = FALSE)
  
  variables_orig <- read.xls(path_variables, as.is = TRUE, check.names = FALSE)
  
  data_orig <- data_orig[data_orig$use == analysis_of & data_orig$TP == "baseline", , drop = FALSE]
  
  
  ### Create metadata
  
  md <- data.frame(patient_id = data_orig[, "ID"], day = data_orig[, "TP"], response = data_orig[, "RESPONSE"], shortname = data_orig[, "shortname"], data = data_orig[, "batches"], stringsAsFactors = FALSE)
  
  md$patient_id <- factor(md$patient_id)
  md$day <- factor(md$day, levels = c("baseline", "TP1"), labels = c("base", "tx"))
  md$data[is.na(md$data)] <- "3"
  md$data <- factor(md$data, levels = c("1", "2", "3"), labels = c("data23", "data29", "data30"))
  md$response <- factor(md$response, levels = c(0, 1), labels = c("NR", "R"))
  
  md$group <- interaction(md$day, md$response, lex.order = TRUE, sep = "_")
  md$shortname <- paste0(md$day, "_", md$shortname)
  md$data_day <- interaction(md$data, md$day, lex.order = TRUE, drop = TRUE)
  
  ### More adjustments
  # We use only the base samples
  md$day <- factor(md$day)
  
}

if(analysis_of == "CyTOFandFACS"){
  
  prefix <- paste0(analysis_of, "_")
  
  data_orig <- read.xls(path_data, as.is = TRUE, check.names = FALSE)
  
  variables_orig <- read.xls(path_variables, as.is = TRUE, check.names = FALSE)
  
  data_orig <- data_orig[data_orig$TP == "baseline", , drop = FALSE]
  
  
  ### Create metadata
  
  md <- data.frame(patient_id = data_orig[, "ID"], day = data_orig[, "TP"], response = data_orig[, "RESPONSE"], shortname = data_orig[, "shortname"], data = data_orig[, "batches"], stringsAsFactors = FALSE)
  
  md$patient_id <- factor(md$patient_id)
  md$day <- factor(md$day, levels = c("baseline", "TP1"), labels = c("base", "tx"))
  md$data[is.na(md$data)] <- "3"
  md$data <- factor(md$data, levels = c("1", "2", "3"), labels = c("data23", "data29", "data30"))
  md$response <- factor(md$response, levels = c(0, 1), labels = c("NR", "R"))
  
  md$group <- interaction(md$day, md$response, lex.order = TRUE, sep = "_")
  md$shortname <- paste0(md$day, "_", md$shortname)
  md$data_day <- interaction(md$data, md$day, lex.order = TRUE, drop = TRUE)
  
  ### More adjustments
  # We use only the base samples
  md$day <- factor(md$day)
  
}

###############################################################################
### Analysis for continous variables
###############################################################################


# ----------------------------------------------------------
### DE analysis for the following variables
# ----------------------------------------------------------

vars_cont <- variables_orig$variable[variables_orig$usage == 1 & variables_orig$continous == 1]
vars_cont

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
ggdf$data <- factor(md$data[mm])
ggdf$response <- factor(md$response[mm])
ggdf$group <- interaction(ggdf$day, ggdf$response, lex.order = TRUE, sep = "_")
ggdf$marker <- factor(ggdf$marker, levels = vars_cont)



# ------------------------------------

nr_marker <- nlevels(ggdf$marker)
nr_in_one_row <- 10
nrow <- ceiling(nr_marker/nr_in_one_row)
hh <- nrow * 2.5
ncol <- ceiling(nr_marker/nrow)
ww <- ncol * 2 + 1

ggp <- ggplot(ggdf, aes(x = group, y = expr, color = group, shape = data)) +
  geom_boxplot(width = 0.9, position = position_dodge(width = 0.95), outlier.colour = NA) +
  geom_point(size=2, alpha = 0.8, position = position_jitterdodge(jitter.width = 0.95, jitter.height = 0, dodge.width = 0.95)) +
  theme_bw() +
  ylab("Expression") +
  xlab("") +
  theme(axis.text.x = element_text(size=10, face="bold"),
    axis.title.y = element_text(size=10, face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
    axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
    legend.position = "right",
    strip.background = element_blank(),) +
  scale_color_manual(values = color_groups) +
  facet_wrap(~ marker, scales = "free", nrow = nrow)


pdf(file.path(outdir, paste0(prefix, "expr_", out_name, "_plot.pdf")), w = ww, h = hh)
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

### LMM when there are base and tx samples
# models2fit <- c("lmer_interglht")
### LM when there are only base samples
models2fit <- c("lm_interglht")


for(k in models2fit){
  # k = "lm_interglht"
  print(k)
  
  switch(k,
    lm_interglht = {
      
      # Fit a LM with interactions + test contrasts with multcomp pckg
      fit_out <- fit_lm_interglht(data = exprc, md, formula = formula_lm, K = K, skippNAs = FALSE)
      
    }, 
    lmer_interglht = {
      
      # Fit a lmer with interactions + test contrasts with multcomp pckg
      fit_out <- fit_lmer_interglht(data = exprc, md, formula = formula_lmer, K = K, skippNAs = FALSE)
      
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
  
  plot_heatmaps_for_sign_expr(expr_all = expr_all, md = md, FDR_cutoff = FDR_cutoff, pval_name2 = pval_name2, adjpval_name2 = adjpval_name2, pval_name_list = pval_name_list, adjpval_name_list = adjpval_name_list, breaks = breaks, legend_breaks = legend_breaks, color_response = color_response, outdir = outdir, prefix = prefix, prefix2 = prefix2, suffix = suffix)
  
  
}


###############################################################################
### Analysis for binary variables
###############################################################################


# ----------------------------------------------------------
### DE analysis for the following variables
# ----------------------------------------------------------

vars_bi <- variables_orig$variable[variables_orig$usage == 1 & variables_orig$continous == 0]
vars_bi

data_bi <- data_orig[, vars_bi, drop = FALSE]

check_bi <- apply(data_bi, 2, function(x){
  nlevels(factor(x))
})

data_bi <- data_bi[, check_bi == 2, drop = FALSE]


data_bi <- apply(data_bi, 2, function(x){
  as.numeric(factor(x)) -1
})

data_bi <- data.frame(t(data_bi))
colnames(data_bi) <- md$shortname

data_bi$marker <- rownames(data_bi)
data_bi$cluster <- ""
data_bi$label <- ""
rownames(data_bi) <- NULL
head(data_bi)


# ----------------------------------------------------------
# Plot features stratified by day and response
# ----------------------------------------------------------

ggdf <- melt(data_bi, id.vars = c("cluster", "label", "marker"), value.name = "expr", variable.name = "sample")

## add group info
mm <- match(ggdf$sample, md$shortname)
ggdf$day <- factor(md$day[mm])
ggdf$data <- factor(md$data[mm])
ggdf$response <- factor(md$response[mm])
ggdf$group <- interaction(ggdf$day, ggdf$response, lex.order = TRUE, sep = "_")
ggdf$marker <- factor(ggdf$marker, levels = vars_bi)

ggdf$expr <- factor(ggdf$expr, levels = c(1, 0))

ggdf <- ggdf[complete.cases(ggdf), ]

# ------------------------------------

nr_marker <- nlevels(ggdf$marker)
nr_in_one_row <- 5
nrow <- ceiling(nr_marker/nr_in_one_row)
hh <- nrow * 2.5
ncol <- ceiling(nr_marker/nrow)
ww <- ncol * 2 + 1

ggp <- ggplot(ggdf, aes(x = interaction(data, group), color = group, fill = expr)) +
  geom_bar() +
  theme_bw() +
  ylab("") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size=10, face="bold"), 
    axis.title.y = element_text(size=10, face="bold"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"), 
    axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
    legend.position = "right",
    strip.text = element_text(size = 10, hjust = 0), strip.background = element_blank()) +
  scale_color_manual(values = color_groups) +
  scale_fill_manual(values = c("gray50", "gray80")) +
  facet_wrap(~ marker, scale = "free", nrow = nrow)

pdf(file.path(outdir, paste0(prefix, "frequencies_", out_name, "_plot.pdf")), w = ww, h = hh)
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
# freq_out <- data_bi[freq_out <- data_bi$marker %in% c("prev.chemo", "bone.met"), , drop = FALSE]

for(k in models2fit){
  # k = "glmer_binomial_interglht_01"
  print(k)
  
  switch(k,
    
    glmer_binomial_interglht_01 = {
      # Fit a GLMM binomial with interactions + test contrasts with multcomp pckg
      fit_out <- fit_glmer_interglht_01(data = freq_out, md = md, family = "binomial", formula = formula_glmer_binomial_01, K = K, skippNAs = FALSE)
      
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
  
  plot_heatmaps_for_sign_expr(expr_all = expr_all, md = md, FDR_cutoff = FDR_cutoff, pval_name2 = pval_name2, adjpval_name2 = adjpval_name2, pval_name_list = pval_name_list, adjpval_name_list = adjpval_name_list, breaks = breaks, legend_breaks = legend_breaks, color_response = color_response, outdir = outdir, prefix = prefix, prefix2 = prefix2, suffix = suffix)
  
}





sessionInfo()











