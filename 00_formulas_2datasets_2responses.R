
### Remove HD from the metadata
md <- md[md$response %in% c("NR", "R"), , drop = FALSE]
md$response <- factor(md$response)


mm <- model.matrix(y ~ response + data_day + response:data_day, data = data.frame(y = 1, md))



if(identical(levels(md$data), c("data23", "data29")) && identical(levels(md$day), c("base", "tx")) && identical(levels(md$response), c("NR", "R"))){
  
  ## create formulas
  formula_lm <- y ~ response + data_day + response:data_day
  formula_lmer <- y ~ response + data_day + response:data_day + (1|patient_id)
  
  formula_glm_binomial <- cbind(y, total-y) ~ response + data_day + response:data_day
  formula_glm_beta <- y/total ~ response + data_day + response:data_day
  formula_glmer_binomial <- y/total ~ response + data_day + response:data_day + (1|patient_id)
  
  ## create contrasts
  contrast_names <- c("NRvsR", "NRvsR_base", "NRvsR_tx", "NRvsR_basevstx")
  k1 <- c(0, 1, 0, 0, 0, 0, 1/2, 0) # NR vs R in base
  k2 <- c(0, 1, 0, 0, 0, 1/2, 0, 1/2) # NR vs R in tx
  k0 <- (k1 + k2) / 2 # NR vs R
  k3 <- c(0, 0, 0, 0, 0, 1/2, 0, 1/2) # whether NR vs R is different in base and tx
  K <- matrix(c(k0, k1, k2, k3), nrow = 4, byrow = TRUE)
  rownames(K) <- contrast_names
  
  ### p-value for sorting the output
  pval_name1 <- "pval_NRvsR"
  ### p-value for plotting the pheatmap2
  adjpval_name2 <- "adjp_NRvsR"
  pval_name2 <- "pval_NRvsR"
  ### p-value for plotting the pheatmap3
  adjpval_name_list <- c("adjp_NRvsR", "adjp_NRvsR_base", "adjp_NRvsR_tx", "adjp_NRvsR_basevstx")
  pval_name_list <- c("pval_NRvsR", "pval_NRvsR_base", "pval_NRvsR_tx", "pval_NRvsR_basevstx")
  
}else{
  
  stop("Metadata does not fit to any the models that are specified !!!")
  
}







