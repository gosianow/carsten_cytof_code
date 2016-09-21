library(lme4) # for fitting mixed models
library(multcomp) # for contrasts glht()
library(glmmADMB) # for glmmadmb()
 
# ---------------------------------------------------------------------------------------
# normal models
# ---------------------------------------------------------------------------------------

# -----------------------------
### Fit a LM
# -----------------------------

fit_lm <- function(data, md){
  
  ### Fit the LM
  fit <- lapply(1:nrow(data), function(i){
    # i = 1
    
    y <- data[i, md$shortname]
    NAs <- is.na(y)
    data_tmp <- data.frame(y = as.numeric(y)[!NAs], md[!NAs, c("day", "response")])
    
    ## there must be at least 10 proportions greater than 0
    if(sum(y[!NAs] > 0) < 10){
      mm <- model.matrix(y ~ response + day, data = data_tmp)
      out <- matrix(NA, nrow = ncol(mm), ncol = 2)
      colnames(out) <- c("coeff", "pval")
      rownames(out) <- colnames(mm)
      return(out)
    }
    
    fit_tmp <- lm(y ~ response + day, data = data_tmp)
    summ_tmp <- summary(fit_tmp)
    
    out <- summ_tmp$coefficients[, c("Estimate", "Pr(>|t|)")]
    colnames(out) <- c("coeff", "pval")
    
    return(out)
    
  })
  
  ### Extract p-values
  pvals <- lapply(fit, function(x){
    x[, "pval"]
  })
  pvals <- do.call(rbind, pvals)
  
  ### Extract fitted coefficients
  coeffs <- lapply(fit, function(x){
    x[, "coeff"]
  })
  coeffs <- do.call(rbind, coeffs)
  
  
  movars <- colnames(pvals)
  
  colnames(pvals) <- paste0("pval_", movars)
  
  ## get adjusted p-values
  adjp <- apply(pvals, 2, p.adjust, method = "BH")
  colnames(adjp) <- paste0("adjp_", movars)
  
  out <- list(pvals = cbind(pvals, adjp), coeffs = coeffs)
  
  return(out)
  
}



# -----------------------------
### Fit a LM with inteactions
# -----------------------------

fit_lm_inter <- function(data, md){
  
  ### Fit the LM
  fit <- lapply(1:nrow(data), function(i){
    # i = 1
    
    y <- data[i, md$shortname]
    NAs <- is.na(y)
    data_tmp <- data.frame(y = as.numeric(y)[!NAs], md[!NAs, c("day", "response")])
    
    ## there must be at least 10 proportions greater than 0
    if(sum(y[!NAs] > 0) < 10){
      mm <- model.matrix(y ~ response + day + response:day, data = data_tmp)
      out <- matrix(NA, nrow = ncol(mm), ncol = 2)
      colnames(out) <- c("coeff", "pval")
      rownames(out) <- colnames(mm)
      return(out)
    }
    
    fit_tmp <- lm(y ~ response + day + response:day, data = data_tmp)
    summ_tmp <- summary(fit_tmp)
    
    out <- summ_tmp$coefficients[, c("Estimate", "Pr(>|t|)")]
    colnames(out) <- c("coeff", "pval")
    
    return(out)
    
  })
  
  ### Extract p-values
  pvals <- lapply(fit, function(x){
    x[, "pval"]
  })
  pvals <- do.call(rbind, pvals)
  
  ### Extract fitted coefficients
  coeffs <- lapply(fit, function(x){
    x[, "coeff"]
  })
  coeffs <- do.call(rbind, coeffs)
  
  movars <- colnames(pvals)
  
  colnames(pvals) <- paste0("pval_", movars)
  
  ## get adjusted p-values
  adjp <- apply(pvals, 2, p.adjust, method = "BH")
  colnames(adjp) <- paste0("adjp_", movars)
  
  out <- list(pvals = cbind(pvals, adjp), coeffs = coeffs)
  
  return(out)
  
}

# -----------------------------
### Fit a LM with inteactions + test contrasts with multcomp pckg
# -----------------------------

fit_lm_interglht <- function(data, md){
  
  ### Fit the LM
  fit <- lapply(1:nrow(data), function(i){
    # i = 1
    
    y <- data[i, md$shortname]
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
      out <- matrix(NA, nrow = nrow(K), ncol = 2)
      colnames(out) <- c("coeff", "pval")
      rownames(out) <- rownames(K)
      return(out)
    }
    
    fit_tmp <- lm(y ~ response + day + response:day, data = data_tmp)
    summary(fit_tmp)
    
    ## fit all contrasts at once
    # contr_tmp <- glht(fit_tmp, linfct = K)
    # summ_tmp <- summary(contr_tmp)
    # 
    # out <- cbind(summ_tmp$test$coefficients, summ_tmp$test$pvalues)
    # colnames(out) <- c("coeff", "pval")
    # out
    
    ## fit contrasts one by one
    out <- t(apply(K, 1, function(k){
      
      contr_tmp <- glht(fit_tmp, linfct = matrix(k, 1))
      summ_tmp <- summary(contr_tmp)
      
      out <- c(summ_tmp$test$coefficients, summ_tmp$test$pvalues)
      names(out) <- c("coeff", "pval")
      return(out)
      
    }))
    out
    
    return(out)
    
  })
  
  ### Extract p-values
  pvals <- lapply(fit, function(x){
    x[, "pval"]
  })
  pvals <- do.call(rbind, pvals)
  
  ### Extract fitted coefficients
  coeffs <- lapply(fit, function(x){
    x[, "coeff"]
  })
  coeffs <- do.call(rbind, coeffs)
  
  movars <- colnames(pvals)
  
  colnames(pvals) <- paste0("pval_", movars)
  
  ## get adjusted p-values
  adjp <- apply(pvals, 2, p.adjust, method = "BH")
  colnames(adjp) <- paste0("adjp_", movars)
  
  out <- list(pvals = cbind(pvals, adjp), coeffs = coeffs)
  
  return(out)
  
}


# -----------------------------
### Fit a lmer with inteactions + test contrasts with multcomp pckg
# -----------------------------

fit_lmer_interglht <- function(data, md){
  
  ### Fit the LM
  fit <- lapply(1:nrow(data), function(i){
    # i = 1
    
    y <- data[i, md$shortname]
    NAs <- is.na(y)
    data_tmp <- data.frame(y = as.numeric(y)[!NAs], md[!NAs, c("day", "response", "patient_id")])
    
    ## create contrasts
    contrast_names <- c("NRvsR", "NRvsR_base", "NRvsR_tx", "NRvsR_basevstx")
    k0 <- c(0, 1, 0, 0, 1/2, 0) # NR vs R
    k1 <- c(0, 1, 0, 0, 0, 0) # NR vs R in base
    k2 <- c(0, 1, 0, 0, 1, 0) # NR vs R in tx
    k3 <- c(0, 0, 0, 0, 1, 0) # whether NR vs R is different in base and tx
    K <- matrix(c(k0, k1, k2, k3), nrow = 4, byrow = TRUE)
    rownames(K) <- contrast_names
    
    ## there must be at least 10 values greater than 0
    if(sum(y[!NAs] > 0) < 10){
      out <- matrix(NA, nrow = nrow(K), ncol = 2)
      colnames(out) <- c("coeff", "pval")
      rownames(out) <- rownames(K)
      return(out)
    }
    
    fit_tmp <- lmer(y ~ response + day + response:day + (1|patient_id), data = data_tmp)
    
    ## get p-values
    ## The summary method for lmer objects doesn’t print p-values for Gaussian mixed models because the degrees of freedom of the t reference distribution are not obvious. However, one can rely on the asymptotic normal distribution for computing univariate p-values for the fixed effects using the cftest function from package multcomp. 
    
    summary(fit_tmp)
    cftest(fit_tmp)
    
    ## fit contrasts one by one
    out <- t(apply(K, 1, function(k){
      # k = K[2, ]
      
      contr_tmp <- glht(fit_tmp, linfct = matrix(k, 1))
      summ_tmp <- summary(contr_tmp)
      summ_tmp
      
      out <- c(summ_tmp$test$coefficients, summ_tmp$test$pvalues)
      names(out) <- c("coeff", "pval")
      return(out)
      
    }))
    out
    
    return(out)
    
  })
  
  ### Extract p-values
  pvals <- lapply(fit, function(x){
    x[, "pval"]
  })
  pvals <- do.call(rbind, pvals)
  
  ### Extract fitted coefficients
  coeffs <- lapply(fit, function(x){
    x[, "coeff"]
  })
  coeffs <- do.call(rbind, coeffs)
  
  movars <- colnames(pvals)
  
  colnames(pvals) <- paste0("pval_", movars)
  
  ## get adjusted p-values
  adjp <- apply(pvals, 2, p.adjust, method = "BH")
  colnames(adjp) <- paste0("adjp_", movars)
  
  out <- list(pvals = cbind(pvals, adjp), coeffs = coeffs)
  
  return(out)
  
}











# ---------------------------------------------------------------------------------------
# logit models
# ---------------------------------------------------------------------------------------

# -----------------------------
### Fit a logit GLM
# -----------------------------

fit_glm_logit <- function(data, md){
  
  ntot <- colSums(data[, md$shortname, drop = FALSE])
  
  ### Fit the GLM
  fit <- lapply(1:nrow(data), function(i){
    # i = 1
    
    y <- data[i, md$shortname]
    NAs <- is.na(y)
    data_tmp <- data.frame(y = as.numeric(y)[!NAs], total = as.numeric(ntot), md[!NAs, c("day", "response")])
    
    ## there must be at least 10 proportions greater than 0
    if(sum(y[!NAs] > 0) < 10){
      mm <- model.matrix(y ~ response + day, data = data_tmp)
      out <- matrix(NA, nrow = ncol(mm), ncol = 2)
      colnames(out) <- c("coeff", "pval")
      rownames(out) <- colnames(mm)
      return(out)
    }
    
    fit_tmp <- glm(cbind(y, total-y) ~ response + day, family = binomial(link = "logit"), data = data_tmp)
    summ_tmp <- summary(fit_tmp)
    
    out <- summ_tmp$coefficients[, c("Estimate", "Pr(>|z|)")]
    colnames(out) <- c("coeff", "pval")
    
    return(out)
    
  })
  
  ### Extract p-values
  pvals <- lapply(fit, function(x){
    x[, "pval"]
  })
  pvals <- do.call(rbind, pvals)
  
  ### Extract fitted coefficients
  coeffs <- lapply(fit, function(x){
    x[, "coeff"]
  })
  coeffs <- do.call(rbind, coeffs)
  
  
  movars <- colnames(pvals)
  
  colnames(pvals) <- paste0("pval_", movars)
  
  ## get adjusted p-values
  adjp <- apply(pvals, 2, p.adjust, method = "BH")
  colnames(adjp) <- paste0("adjp_", movars)
  
  out <- list(pvals = cbind(pvals, adjp), coeffs = coeffs)
  
  return(out)
  
}



# -----------------------------
### Fit a logit GLM with inteactions
# -----------------------------

fit_glm_logit_inter <- function(data, md){
  
  ntot <- colSums(data[, md$shortname, drop = FALSE])
  
  ### Fit the GLM
  fit <- lapply(1:nrow(data), function(i){
    # i = 1
    
    y <- data[i, md$shortname]
    NAs <- is.na(y)
    data_tmp <- data.frame(y = as.numeric(y)[!NAs], total = as.numeric(ntot), md[!NAs, c("day", "response")])
    
    ## there must be at least 10 proportions greater than 0
    if(sum(y[!NAs] > 0) < 10){
      mm <- model.matrix(y ~ response + day + response:day, data = data_tmp)
      out <- matrix(NA, nrow = ncol(mm), ncol = 2)
      colnames(out) <- c("coeff", "pval")
      rownames(out) <- colnames(mm)
      return(out)
    }
    
    fit_tmp <- glm(cbind(y, total-y) ~ response + day + response:day, family = binomial(link = "logit"), data = data_tmp)
    summ_tmp <- summary(fit_tmp)
    
    out <- summ_tmp$coefficients[, c("Estimate", "Pr(>|z|)")]
    colnames(out) <- c("coeff", "pval")
    
    return(out)
    
  })
  
  ### Extract p-values
  pvals <- lapply(fit, function(x){
    x[, "pval"]
  })
  pvals <- do.call(rbind, pvals)
  
  ### Extract fitted coefficients
  coeffs <- lapply(fit, function(x){
    x[, "coeff"]
  })
  coeffs <- do.call(rbind, coeffs)
  
  movars <- colnames(pvals)
  
  colnames(pvals) <- paste0("pval_", movars)
  
  ## get adjusted p-values
  adjp <- apply(pvals, 2, p.adjust, method = "BH")
  colnames(adjp) <- paste0("adjp_", movars)
  
  out <- list(pvals = cbind(pvals, adjp), coeffs = coeffs)
  
  return(out)
  
}


# -----------------------------
### Fit a logit GLM with inteactions + test contrasts with multcomp pckg
# -----------------------------

fit_glm_logit_interglht <- function(data, md, family = "binomial"){
  
  ntot <- colSums(data[, md$shortname, drop = FALSE])
  
  ### Fit the GLM
  fit <- lapply(1:nrow(data), function(i){
    # i = 1
    
    y <- data[i, md$shortname]
    NAs <- is.na(y)
    data_tmp <- data.frame(y = as.numeric(y)[!NAs], total = as.numeric(ntot), md[!NAs, c("day", "response")])
    
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
      out <- matrix(NA, nrow = nrow(K), ncol = 2)
      colnames(out) <- c("coeff", "pval")
      rownames(out) <- rownames(K)
      return(out)
    }
    
    
    switch(family,
      binomial = {
        fit_tmp <- glm(cbind(y, total-y) ~ response + day + response:day, family = binomial(link = "logit"), data = data_tmp)
        summary(fit_tmp)
      },
      quasibinomial = {
        fit_tmp <- glm(cbind(y, total-y) ~ response + day + response:day, family = quasibinomial(link = "logit"), data = data_tmp)
        summary(fit_tmp)
      }
    )
    
    ## fit all contrasts at once
    # contr_tmp <- glht(fit_tmp, linfct = K)
    # summ_tmp <- summary(contr_tmp)
    # 
    # out <- cbind(summ_tmp$test$coefficients, summ_tmp$test$pvalues)
    # colnames(out) <- c("coeff", "pval")
    # out
    
    ## fit contrasts one by one
    out <- t(apply(K, 1, function(k){
      
      contr_tmp <- glht(fit_tmp, linfct = matrix(k, 1))
      summ_tmp <- summary(contr_tmp)
      
      out <- c(summ_tmp$test$coefficients, summ_tmp$test$pvalues)
      names(out) <- c("coeff", "pval")
      return(out)
      
    }))
    out
    
    
    return(out)
    
  })
  
  ### Extract p-values
  pvals <- lapply(fit, function(x){
    x[, "pval"]
  })
  pvals <- do.call(rbind, pvals)
  
  ### Extract fitted coefficients
  coeffs <- lapply(fit, function(x){
    x[, "coeff"]
  })
  coeffs <- do.call(rbind, coeffs)
  
  movars <- colnames(pvals)
  
  colnames(pvals) <- paste0("pval_", movars)
  
  ## get adjusted p-values
  adjp <- apply(pvals, 2, p.adjust, method = "BH")
  colnames(adjp) <- paste0("adjp_", movars)
  
  out <- list(pvals = cbind(pvals, adjp), coeffs = coeffs)
  
  return(out)
  
}


# -----------------------------
### Fit a logit GLMM with inteactions + test contrasts with multcomp pckg
# -----------------------------


fit_glmer_logit_interglht <- function(data, md, method = "glmer"){
  
  ntot <- colSums(data[, md$shortname, drop = FALSE])
  
  ### Fit the GLM
  fit <- lapply(1:nrow(data), function(i){
    # i = 2
    
    y <- data[i, md$shortname]
    NAs <- is.na(y)
    data_tmp <- data.frame(y = as.numeric(y)[!NAs], total = as.numeric(ntot), md[!NAs, c("day", "response", "patient_id")])
    
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
      out <- matrix(NA, nrow = nrow(K), ncol = 2)
      colnames(out) <- c("coeff", "pval")
      rownames(out) <- rownames(K)
      return(out)
    }

    switch(method, 
      glmer = {
        fit_tmp <- glmer(y/total ~ response + day + response:day + (1|patient_id), weights = total, family = binomial, data = data_tmp)
        summary(fit_tmp)
        },
      glmmadmb_mixed_beta = {
        data_tmp$prop <- data_tmp$y/data_tmp$total
        fit_tmp <- glmmadmb(prop ~ response + day + response:day + (1|patient_id), family = "beta", data = data_tmp)
        summary(fit_tmp)
        },
      glmmadmb_mixed_betabinomial = {
        fit_tmp <- glmmadmb(cbind(y, total-y) ~ response + day + response:day + (1|patient_id), family = "betabinomial", data = data_tmp)
        summary(fit_tmp)
        },
      glmmadmb_fixed_beta = {
        fit_tmp <- glmmadmb(y/total ~ response + day + response:day, family = "beta", data = data_tmp)
        summary(fit_tmp)
      },
      glmmadmb_fixed_betabinomial = {
        fit_tmp <- glmmadmb(cbind(y, total-y) ~ response + day + response:day, family = "betabinomial", data = data_tmp)
        summary(fit_tmp)
        }
      )
    
    
    ## fit contrasts one by one
    out <- t(apply(K, 1, function(k){
      # k = K[2, ]
      
      contr_tmp <- glht(fit_tmp, linfct = matrix(k, 1))
      summ_tmp <- summary(contr_tmp)
      summ_tmp
      
      out <- c(summ_tmp$test$coefficients, summ_tmp$test$pvalues)
      names(out) <- c("coeff", "pval")
      return(out)
      
    }))
    out
    
    return(out)
    
  })
  
  ### Extract p-values
  pvals <- lapply(fit, function(x){
    x[, "pval"]
  })
  pvals <- do.call(rbind, pvals)
  
  ### Extract fitted coefficients
  coeffs <- lapply(fit, function(x){
    x[, "coeff"]
  })
  coeffs <- do.call(rbind, coeffs)
  
  movars <- colnames(pvals)
  
  colnames(pvals) <- paste0("pval_", movars)
  
  ## get adjusted p-values
  adjp <- apply(pvals, 2, p.adjust, method = "BH")
  colnames(adjp) <- paste0("adjp_", movars)
  
  out <- list(pvals = cbind(pvals, adjp), coeffs = coeffs)
  
  return(out)
  
}













