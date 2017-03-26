
### Remove HD and tx samples from metadata
md <- md[md$response %in% c("NR", "R"), , drop = FALSE]
md <- md[md$day == "base", , drop = FALSE]
md$day <- factor(md$day)
md$data_day <- factor(md$data_day)
md$response <- factor(md$response)


model.matrix( ~ response + data + response:data, data = md)
model.matrix( ~ response + data, data = md)


if(identical(levels(md$data), c("data23", "data29")) && identical(levels(md$day), c("base")) && identical(levels(md$response), c("NR", "R")) && identical(levels(md$data_day), c("data23.base", "data29.base"))){
  
  ## create formulas
  formula_lm <- y ~ response + data
  formula_lmer <- y ~ response + data + (1|patient_id)
  
  formula_glm_binomial <- cbind(y, total-y) ~ response + data
  formula_glm_beta <- y/total ~ response + data
  formula_glmer_binomial <- y/total ~ response + data + (1|patient_id)
  
  ## create contrasts
  contrast_names <- c("NRvsR_base")
  k1 <- c(0, 1, 0) # NR vs R in base
  K <- matrix(c(k1), nrow = 1, byrow = TRUE)
  rownames(K) <- contrast_names
  
  ### p-value for sorting the output
  pval_name1 <- "pval_NRvsR_base"
  ### p-value for plotting the pheatmap2
  adjpval_name2 <- NULL
  pval_name2 <- NULL
  ### p-value for plotting the pheatmap3
  adjpval_name_list <- c("adjp_NRvsR_base")
  pval_name_list <- c("pval_NRvsR_base")
  
}else{
  
  stop("Metadata does not fit to any the models that are specified !!!")
  
}







