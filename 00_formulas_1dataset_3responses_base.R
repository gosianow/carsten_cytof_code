
if(identical(levels(md$day), "base") && identical(levels(md$response), c("NR", "R", "HD"))){
  ## create formulas
  formula_lm <- y ~ response
  formula_lmer <- y ~ response + (1|patient_id)
  
  formula_glm_binomial <- cbind(y, total-y) ~ response
  formula_glm_beta <- y/total ~ response
  formula_glmer_binomial <- y/total ~ response + (1|patient_id)
  formula_glmer_binomial_01 <- y ~ response + (1|patient_id)
  formula_glm_beta_01 <- y ~ response
  
  ## create contrasts
  contrast_names <- c("NRvsR_base")
  k1 <- c(0, 1, 0) # NR vs R in base
  K <- matrix(k1, nrow = 1, byrow = TRUE)
  rownames(K) <- contrast_names
  
  ### p-value for sorting the output
  pval_name1 <- "pval_NRvsR_base"
  ### p-value for plotting the pheatmap2
  adjpval_name2 <- NULL
  pval_name2 <- NULL
  ### p-value for plotting the pheatmap3
  adjpval_name_list <- "adjp_NRvsR_base"
  pval_name_list <- "pval_NRvsR_base"
  
  
}else{
  
  stop("Metadata does not fit to any the models that are specified !!!")  
  
}
