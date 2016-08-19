library(lme4) # for fitting mixed models
library(multcomp)


# -----------------------------
### Run two-way ANOVA
# -----------------------------

fit_two_way_anova <- function(data, md){
  
  pvs_anova <- t(apply(data[, md$shortname], 1, function(y){
    # y <- data[1, md$shortname]
    
    ## there must be at least 10 proportions greater than 0
    if(sum(y > 0) < 10)
      return(rep(NA, 3))
    
    data_tmp <- data.frame(y = as.numeric(y), md[, c("day", "response")])
    
    res_tmp <- aov(y ~ day * response, data = data_tmp)
    summ_tmp <- summary(res_tmp)
    
    summ_tmp[[1]][1:3, "Pr(>F)"]
    
  }))
  
  movars <- c("day", "response", "day:response")
  colnames(pvs_anova) <- paste0("pval_", movars)
  pvs_anova <- data.frame(pvs_anova)
  
  ## get adjusted p-values
  adjp_anova <- data.frame(apply(pvs_anova, 2, p.adjust, method = "BH"))
  colnames(adjp_anova) <- paste0("adjp_", movars)
  
  pvs_anova_out <- data.frame(group = rownames(pvs_anova), pvs_anova, adjp_anova)
  
  return(pvs_anova_out)
  
}


# -----------------------------
### Fit a normal GLM
# -----------------------------

fit_glm_norm <- function(data, md){
  
  pvs_glm <- t(apply(data[, md$shortname], 1, function(y){
    # y <- data[1, md$shortname]
    
    ## there must be at least 10 proportions greater than 0
    if(sum(y > 0) < 10)
      return(rep(NA, 4))
    
    data_tmp <- data.frame(y = as.numeric(y), md[, c("day", "response")])
    
    res_tmp <- glm(y ~ response + day, data = data_tmp)
    
    sum_tmp <- summary(res_tmp)
    
    out <- as.numeric(sum_tmp$coefficients[, "Pr(>|t|)"])
    
    return(out)
    
  }))
  
  data_tmp <- data.frame(y = as.numeric(data[1, md$shortname]), md[, c("day", "response")])
  momat <- model.matrix(y ~ response + day, data = data_tmp)
  movars <- colnames(momat)
  
  colnames(pvs_glm) <- paste0("pval_", movars)
  pvs_glm <- data.frame(pvs_glm)
  
  
  ## get adjusted p-values
  
  adjp_glm <- data.frame(apply(pvs_glm, 2, p.adjust, method = "BH"))
  colnames(adjp_glm) <- paste0("adjp_", movars)
  
  ## save the results
  pvs_glm_out <- data.frame(group = rownames(pvs_glm), pvs_glm, adjp_glm)
  
  return(pvs_glm_out)
  
}


# -----------------------------
### Fit a normal GLM with inteactions
# -----------------------------

fit_glm_norm_inter <- function(data, md){
  
  
  pvs_glm <- t(apply(data[, md$shortname], 1, function(y){
    # y <- data[1, md$shortname]
    
    ## there must be at least 10 proportions greater than 0
    if(sum(y > 0) < 10)
      return(rep(NA, 4))
    
    data_tmp <- data.frame(y = as.numeric(y), md[, c("day", "response")])
    
    res_tmp <- glm(y ~ response + day + response:day, data = data_tmp)
    
    sum_tmp <- summary(res_tmp)
    
    out <- as.numeric(sum_tmp$coefficients[, "Pr(>|t|)"])
    
    return(out)
    
  }))
  
  data_tmp <- data.frame(y = as.numeric(data[1, md$shortname]), md[, c("day", "response")])
  momat <- model.matrix(y ~ response + day + response:day, data = data_tmp)
  movars <- colnames(momat)
  
  colnames(pvs_glm) <- paste0("pval_", movars)
  pvs_glm <- data.frame(pvs_glm)
  
  
  ## get adjusted p-values
  
  adjp_glm <- data.frame(apply(pvs_glm, 2, p.adjust, method = "BH"))
  colnames(adjp_glm) <- paste0("adjp_", movars)
  
  ## save the results
  pvs_glm_out <- data.frame(cluster = rownames(pvs_glm), pvs_glm, adjp_glm)
  
  return(pvs_glm_out)
  
}

