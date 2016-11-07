library(lme4) # for fitting mixed models
library(multcomp) # for contrasts glht()
library(glmmADMB) # for glmmadmb()
library(robustbase) # glmrob
library(robust) # glmRob
library(MASS) # rlm


# -----------------------------
# Wilcoxon / Mann-Whitney U test
# -----------------------------

test_wilcoxon <- function(data, md){
  
  ## Keep only the NR and R samples
  md <- md[md$response %in% c("NR", "R"), ]
  md$response <- factor(md$response)
  days <- levels(md$day)
  contrast_names <- paste0("NRvsR_", days)
  
  ### Fit the LM
  fit <- lapply(1:nrow(data), function(i){
    # i = 25
    y <- data[i, md$shortname]
    NAs <- is.na(y)
    data_tmp <- data.frame(y = as.numeric(y), md)[!NAs, , drop = FALSE]
    
    out <- matrix(NA, nrow = length(contrast_names), ncol = 2)
    colnames(out) <- c("coeff", "pval")
    rownames(out) <- contrast_names
    
    data_tmp$response <- factor(data_tmp$response)
    if(nlevels(data_tmp$response) < 2){
      return(out)
    }
    
    ## there must be at least 10 values different than 0
    if(sum(y[!NAs] != 0) < 10 || any(y[!NAs] %in% c(-Inf, Inf))){
      return(out)
    }
    
    for(j in 1:nlevels(data_tmp$day)){
      # j = 2
      
      test_tmp <- wilcox.test(formula = y ~ response, data = data_tmp, subset = data_tmp$day == days[j])
      out[paste0("NRvsR_", days[j]), "pval"] <- test_tmp$p.value
      
    }
    
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
  
  movars <- contrast_names
  
  colnames(pvals) <- paste0("pval_", movars)
  colnames(coeffs) <- movars
  
  ## get adjusted p-values
  adjp <- apply(pvals, 2, p.adjust, method = "BH")
  colnames(adjp) <- paste0("adjp_", movars)
  
  out <- list(pvals = cbind(pvals, adjp), coeffs = coeffs)
  
  return(out)
  
}



# ---------------------------------------------------------------------------------------
# normal models
# ---------------------------------------------------------------------------------------

# -----------------------------
### Fit a LM with inteactions + test contrasts with multcomp pckg
# -----------------------------

fit_lm_interglht <- function(data, md, method = "lm", formula, K){
  
  contrast_names <- rownames(K)
  
  ### Fit the LM
  fit <- lapply(1:nrow(data), function(i){
    # i = 1
    
    y <- data[i, md$shortname]
    NAs <- is.na(y)
    data_tmp <- data.frame(y = as.numeric(y), md)[!NAs, , drop = FALSE]
    
    mm <- model.matrix(formula, data = data_tmp)
    
    ## there must be at least 10 values different than 0
    out <- matrix(NA, nrow = length(contrast_names), ncol = 2)
    colnames(out) <- c("coeff", "pval")
    rownames(out) <- contrast_names
    
    if(sum(y[!NAs] != 0) < 10 || any(y[!NAs] %in% c(-Inf, Inf)))
      return(out)
    
    ## do not anlyse a cluster with NAs; for merged data it means such cluster was not present in all the datasets
    if(any(NAs))
      return(out)
      
      
    switch(method, 
      lm = {
        fit_tmp <- lm(formula, data = data_tmp)
        summary(fit_tmp)
      },
      lmrob = {
        fit_tmp <- robustbase::lmrob(formula, data = data_tmp, setting = "KS2014")
        summary(fit_tmp)
      },
      rlm = {
        fit_tmp <- rlm(formula, data = data_tmp)
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
    
    colnames(out) <- c("coeff", "pval")
    rownames(out) <- contrast_names
    
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
  
  movars <- contrast_names
  
  colnames(pvals) <- paste0("pval_", movars)
  colnames(coeffs) <- movars
  
  ## get adjusted p-values
  adjp <- apply(pvals, 2, p.adjust, method = "BH")
  colnames(adjp) <- paste0("adjp_", movars)
  
  out <- list(pvals = cbind(pvals, adjp), coeffs = coeffs)
  
  return(out)
  
}


# -----------------------------
### Fit a lmer with inteactions + test contrasts with multcomp pckg
# -----------------------------

fit_lmer_interglht <- function(data, md, formula, K){
  
  contrast_names <- rownames(K)
  
  ### Fit the LM
  fit <- lapply(1:nrow(data), function(i){
    # i = 1
    
    y <- data[i, md$shortname]
    NAs <- is.na(y)
    data_tmp <- data.frame(y = as.numeric(y)[!NAs], md[!NAs, ])
    
    ## there must be at least 10 values different than 0
    out <- matrix(NA, nrow = length(contrast_names), ncol = 2)
    colnames(out) <- c("coeff", "pval")
    rownames(out) <- contrast_names
    
    if(sum(y[!NAs] != 0) < 10 || any(y[!NAs] %in% c(-Inf, Inf)))
      return(out)
    
    ## do not anlyse a cluster with NAs; for merged data it means such cluster was not present in all the datasets
    if(any(NAs))
      return(out)
    
    fit_tmp <- lmer(formula, data = data_tmp)
    
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
  
  movars <- contrast_names
  
  colnames(pvals) <- paste0("pval_", movars)
  colnames(coeffs) <- movars
  
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
### Fit a logit GLM with inteactions + test contrasts with multcomp pckg
# -----------------------------

fit_glm_interglht <- function(data, md, family = "binomial", formula, K){
  
  ntot <- colSums(data[, md$shortname, drop = FALSE], na.rm = TRUE)
  contrast_names <- rownames(K)
  
  ### Fit the GLM
  fit <- lapply(1:nrow(data), function(i){
    # i = 1
    
    y <- data[i, md$shortname]
    NAs <- is.na(y)
    data_tmp <- data.frame(y = as.numeric(y), total = as.numeric(ntot), md)[!NAs, , drop = FALSE]
    
    mm <- model.matrix(formula, data = data_tmp)
    
    ## there must be at least 10 proportions greater than 0
    out <- matrix(NA, nrow = length(contrast_names), ncol = 2)
    colnames(out) <- c("coeff", "pval")
    rownames(out) <- contrast_names
    
    if(sum(y[!NAs] > 0) < 10)
      return(out)
    
    ## do not anlyse a cluster with NAs; for merged data it means such cluster was not present in all the datasets
    if(any(NAs))
      return(out)
    
    switch(family,
      binomial = {
        fit_tmp <- glm(formula, family = binomial(link = "logit"), data = data_tmp)
        summary(fit_tmp)
      },
      quasibinomial = {
        fit_tmp <- glm(formula, family = quasibinomial(link = "logit"), data = data_tmp)
        summary(fit_tmp)
      },
      beta = {
        fit_tmp <- NULL
        try(fit_tmp <- glmmadmb(formula, family = "beta", data = data_tmp), silent = TRUE)
        summary(fit_tmp)
      },
      betabinomial = {
        fit_tmp <- NULL
        try(fit_tmp <- glmmadmb(formula, family = "betabinomial", data = data_tmp), silent = TRUE)
        summary(fit_tmp)
      },
      binomial_rob = {
        fit_tmp <- glmrob(formula, family = binomial, data = data_tmp)
        summary(fit_tmp)
      }
      
    )
    
    if(is.null(fit_tmp))
      return(out)
    
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
  
  movars <- contrast_names
  
  colnames(pvals) <- paste0("pval_", movars)
  colnames(coeffs) <- movars
  
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


fit_glmer_interglht <- function(data, md, family = "binomial", formula, K){
  
  ntot <- colSums(data[, md$shortname, drop = FALSE], na.rm = TRUE)
  contrast_names <- rownames(K)
  
  ### Fit the GLM
  fit <- lapply(1:nrow(data), function(i){
    # i = 1
    print(i)
    
    y <- data[i, md$shortname]
    NAs <- is.na(y)
    data_tmp <- data.frame(y = as.numeric(y), total = as.numeric(ntot), md)[!NAs, , drop = FALSE]
    
    ## there must be at least 10 proportions greater than 0
    out <- matrix(NA, nrow = length(contrast_names), ncol = 2)
    colnames(out) <- c("coeff", "pval")
    rownames(out) <- contrast_names
    
    if(sum(y[!NAs] > 0) < 10)
      return(out)
    
    ## do not anlyse a cluster with NAs; for merged data it means such cluster was not present in all the datasets
    if(any(NAs))
      return(out)
    
    switch(family, 
      binomial = {
        fit_tmp <- glmer(formula, weights = total, family = binomial, data = data_tmp)
        summary(fit_tmp)
      },
      beta = {
        data_tmp$prop <- data_tmp$y/data_tmp$total
        fit_tmp <- NULL
        try(fit_tmp <- glmmadmb(formula, family = "beta", data = data_tmp), silent = TRUE)
        summary(fit_tmp)
      },
      betabinomial = {
        fit_tmp <- NULL
        try(fit_tmp <- glmmadmb(formula, family = "betabinomial", data = data_tmp), silent = TRUE)
        summary(fit_tmp)
      }
      
    )
    
    if(is.null(fit_tmp))
      return(out)
    
    ## fit contrasts one by one
    out <- t(apply(K, 1, function(k){
      # k = K[2, ]
      
      contr_tmp <- glht(fit_tmp, linfct = matrix(k, 1))
      summ_tmp <- summary(contr_tmp)
      summ_tmp
      
      out <- c(summ_tmp$test$coefficients, summ_tmp$test$pvalues)
      names(out) <- c("coeff", "pval")
      out
      
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
  
  movars <- contrast_names
  
  colnames(pvals) <- paste0("pval_", movars)
  colnames(coeffs) <- movars
  
  colnames(pvals) <- paste0("pval_", movars)
  
  ## get adjusted p-values
  adjp <- apply(pvals, 2, p.adjust, method = "BH")
  colnames(adjp) <- paste0("adjp_", movars)
  
  out <- list(pvals = cbind(pvals, adjp), coeffs = coeffs)
  
  return(out)
  
}













