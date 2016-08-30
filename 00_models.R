library(lme4) # for fitting mixed models
library(multcomp) # for contrasts glht()


# -----------------------------
### Run two-way ANOVA
# -----------------------------

fit_two_way_anova <- function(data, md){
  
  pvs_anova <- t(apply(data[, md$shortname], 1, function(y){
    # y <- data[1, md$shortname]
    
    NAs <- is.na(y)
    
    ## there must be at least 10 proportions greater than 0
    if(sum(y[!NAs] > 0) < 10)
      return(rep(NA, 2))
    
    data_tmp <- data.frame(y = as.numeric(y)[!NAs], md[!NAs, c("day", "response")])
    
    res_tmp <- aov(y ~ day + response, data = data_tmp)
    summ_tmp <- summary(res_tmp)
    
    summ_tmp[[1]][1:2, "Pr(>F)"]
    
  }))
  
  movars <- c("day", "response")
  colnames(pvs_anova) <- paste0("pval_", movars)
  pvs_anova <- data.frame(pvs_anova)
  
  ## get adjusted p-values
  adjp_anova <- data.frame(apply(pvs_anova, 2, p.adjust, method = "BH"))
  colnames(adjp_anova) <- paste0("adjp_", movars)
  
  pvs_out <- data.frame(pvs_anova, adjp_anova)
  
  return(pvs_out)
  
}


# -----------------------------
### Fit a normal GLM
# -----------------------------

fit_glm_norm <- function(data, md){
  
  pvs <- t(apply(data[, md$shortname], 1, function(y){
    # y <- data[1, md$shortname]
    
    NAs <- is.na(y)
    
    data_tmp <- data.frame(y = as.numeric(y)[!NAs], md[!NAs, c("day", "response")])
    
    ## there must be at least 10 proportions greater than 0
    if(sum(y[!NAs] > 0) < 10){
      mm <- model.matrix(y ~ response + day, data = data_tmp)
      out <- rep(NA, ncol(mm))
      names(out) <- colnames(mm)
      return(out)
    }
    
    res_tmp <- glm(y ~ response + day, data = data_tmp)
    summ_tmp <- summary(res_tmp)
    
    out <- as.numeric(summ_tmp$coefficients[, "Pr(>|t|)"])
    names(out) <- rownames(summ_tmp$coefficients)
    
    return(out)
    
  }))
  
  movars <- colnames(pvs)
  
  colnames(pvs) <- paste0("pval_", movars)
  pvs <- data.frame(pvs)
  
  ## get adjusted p-values
  adjp <- data.frame(apply(pvs, 2, p.adjust, method = "BH"))
  colnames(adjp) <- paste0("adjp_", movars)
  
  ## save the results
  pvs_out <- data.frame(pvs, adjp)
  
  return(pvs_out)
  
}


# -----------------------------
### Fit a normal GLM - with response only
# -----------------------------

fit_glm_norm_resp <- function(data, md){
  
  pvs <- t(apply(data[, md$shortname], 1, function(y){
    # y <- data[1, md$shortname]
    
    NAs <- is.na(y)
    
    data_tmp <- data.frame(y = as.numeric(y)[!NAs], md[!NAs, c("response"), drop = FALSE])
    
    ## there must be at least 10 proportions greater than 0
    if(sum(y[!NAs] > 0) < 10){
      mm <- model.matrix(y ~ response, data = data_tmp)
      out <- rep(NA, ncol(mm))
      names(out) <- colnames(mm)
      return(out)
    }
    
    res_tmp <- glm(y ~ response, data = data_tmp)
    summ_tmp <- summary(res_tmp)
    
    out <- as.numeric(summ_tmp$coefficients[, "Pr(>|t|)"])
    names(out) <- rownames(summ_tmp$coefficients)
    
    return(out)
    
  }))
  
  movars <- colnames(pvs)
  
  colnames(pvs) <- paste0("pval_", movars)
  pvs <- data.frame(pvs)
  
  ## get adjusted p-values
  adjp <- data.frame(apply(pvs, 2, p.adjust, method = "BH"))
  colnames(adjp) <- paste0("adjp_", movars)
  
  ## save the results
  pvs_out <- data.frame(pvs, adjp)
  
  return(pvs_out)
  
}


# -----------------------------
### Fit a normal GLM with inteactions
# -----------------------------

fit_glm_norm_inter <- function(data, md){
  
  pvs <- t(apply(data[, md$shortname], 1, function(y){
    # y <- data[1, md$shortname]
    
    NAs <- is.na(y)
    
    data_tmp <- data.frame(y = as.numeric(y)[!NAs], md[!NAs, c("day", "response")])
    
    ## create contrasts
    contrast_names <- c("NRvsR", "NRvsR_base", "NRvsR_tx", "NRvsR_basevstx")
    k0 <- c(0, 1, 0, 0, 1/2, 0) # NR vs R
    k1 <- c(0, 1, 0, 0, 0, 0) # NR vs R in base
    k2 <- c(0, 1, 0, 0, 1, 0) # NR vs R in tx
    k3 <- c(0, 0, 0, 0, 1, 0) # whether NR vs R is different in base and tx
    K <- matrix(c(k0, k1, k2, k3), nrow = 4, byrow = TRUE)
    rownames(K) <- contrast_names
    
    ## there must be at least 10 proportions greater than 0
    if(sum(y[!NAs] > 0) < 10){
      out <- rep(NA, nrow(K))
      names(out) <- rownames(K)
      return(out)
    }
    
    model_tmp <- glm(y ~ response + day + response:day, data = data_tmp)
    
    summ_tmp <- summary(model_tmp)
    
    ## fit contrasts
    contr_tmp <- glht(model_tmp, linfct = K)
    summ_tmp <- summary(contr_tmp)
    
    out <- summ_tmp$test$pvalues
    names(out) <- contrast_names
    
    return(out)
    
  }))
  
  movars <- colnames(pvs)
  
  colnames(pvs) <- paste0("pval_", movars)
  pvs <- data.frame(pvs)
  
  ## get adjusted p-values
  adjp <- data.frame(apply(pvs, 2, p.adjust, method = "BH"))
  colnames(adjp) <- paste0("adjp_", movars)
  
  pvs_out <- data.frame(pvs, adjp)
  
  return(pvs_out)
  
}



# -----------------------------
### Fit a logit GLM
# -----------------------------


fit_glm_logit <- function(data, md){
  
  ntot <- colSums(data[, md$shortname, drop = FALSE])
  
  pvs <- t(apply(data[, md$shortname], 1, function(y){
    # y <- data[1, md$shortname]
    
    NAs <- is.na(y)
    
    data_tmp <- data.frame(y = as.numeric(y)[!NAs], total = as.numeric(ntot), md[!NAs, c("day", "response")])
    
    if(sum(y[!NAs] > 0) < 10){
      mm <- model.matrix(y ~ response + day, data = data_tmp)
      out <- rep(NA, ncol(mm))
      names(out) <- colnames(mm)
      return(out)
    }
    
    res_tmp <- glm(cbind(y, total-y) ~ response + day, family = binomial(link = "logit"), data = data_tmp)
    
    summ_tmp <- summary(res_tmp)
    
    out <- as.numeric(summ_tmp$coefficients[, "Pr(>|z|)"])
    names(out) <- rownames(summ_tmp$coefficients)
    
    return(out)
    
  }))
  
  movars <- colnames(pvs)
  
  colnames(pvs) <- paste0("pval_", movars)
  pvs <- data.frame(pvs)
  
  
  ## get adjusted p-values
  adjp <- data.frame(apply(pvs, 2, p.adjust, method = "BH"))
  colnames(adjp) <- paste0("adjp_", movars)
  
  ## save the results
  pvs_out <- data.frame(pvs, adjp)
  
  return(pvs_out)
  
}



# -----------------------------
### Fit a logit GLMM 
# -----------------------------

fit_glmm_logit <- function(data, md){
  
  ntot <- colSums(data[, md$shortname, drop = FALSE])
  
  pvs <- t(apply(data[, md$shortname], 1, function(y){
    # y <- data[1, md$shortname]

    NAs <- is.na(y)
    
    ## prepare binomial/Bernoulli data
    data_tmp <- data.frame(y = as.numeric(y)[!NAs], total = as.numeric(ntot), md[!NAs, c("day", "response", "patient_id")])
    
    ## there must be at least 10 proportions greater than 0
    if(sum(y[!NAs] > 0) < 10){
      mm <- model.matrix(y ~ response + day, data = data_tmp)
      out <- rep(NA, ncol(mm))
      names(out) <- colnames(mm)
      return(out)
    }
    
    res_tmp <- glmer(y/total ~ response + day + (1|patient_id), weights = total, family = binomial, data = data_tmp)
    
    summ_tmp <- summary(res_tmp)
    
    out <- as.numeric(summ_tmp$coefficients[, "Pr(>|z|)"])
    names(out) <- rownames(summ_tmp$coefficients)
    
    return(out)
    
  }))
  
  movars <- colnames(pvs)
  
  colnames(pvs) <- paste0("pval_", movars)
  pvs <- data.frame(pvs)
  
  
  ## get adjusted p-values
  adjp <- data.frame(apply(pvs, 2, p.adjust, method = "BH"))
  colnames(adjp) <- paste0("adjp_", movars)
  
  ## save the results
  pvs_out <- data.frame(pvs, adjp)
  
  return(pvs_out)
  
}













