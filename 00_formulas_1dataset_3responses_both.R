


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
