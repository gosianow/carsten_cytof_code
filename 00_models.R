library(lme4) # for fitting mixed models
library(multcomp) # for contrasts glht()

# ---------------------------------------------------------------------------------------
# normal models
# ---------------------------------------------------------------------------------------

# -----------------------------
# Fit a LM with inteactions + test contrasts with multcomp pckg
# -----------------------------

fit_lm_interglht <- function(data, md, formula, K, skippNAs = TRUE){
  
  contrast_names <- rownames(K)
  
  ### Fit the LM
  fit <- lapply(1:nrow(data), function(i){
    # i = 1
    
    y <- data[i, md$shortname]
    NAs <- is.na(y)
    data_tmp <- data.frame(y = as.numeric(y), md)[!NAs, , drop = FALSE]
    
    mm <- model.matrix(formula, data = data_tmp)
    
    out <- matrix(NA, nrow = length(contrast_names), ncol = 2)
    colnames(out) <- c("coeff", "pval")
    rownames(out) <- contrast_names
    
    if(any(y[!NAs] %in% c(-Inf, Inf)) || sum(y[!NAs] == 0) > length(y[!NAs])/2)
      return(out)
    
    ## do not anlyse a cluster with NAs; for merged data it means such cluster was not present in all the datasets
    if(any(NAs) && skippNAs)
      return(out)

    fit_tmp <- lm(formula, data = data_tmp)
    summary(fit_tmp)
    
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
# Fit a lmer with inteactions + test contrasts with multcomp pckg
# -----------------------------

fit_lmer_interglht <- function(data, md, formula, K, skippNAs = TRUE){
  
  contrast_names <- rownames(K)
  
  ### Fit the LM
  fit <- lapply(1:nrow(data), function(i){
    # i = 1
    
    y <- data[i, md$shortname]
    NAs <- is.na(y)
    data_tmp <- data.frame(y = as.numeric(y)[!NAs], md[!NAs, ])
    
    out <- matrix(NA, nrow = length(contrast_names), ncol = 2)
    colnames(out) <- c("coeff", "pval")
    rownames(out) <- contrast_names
    
    if(any(y[!NAs] %in% c(-Inf, Inf)) || sum(y[!NAs] == 0) > length(y[!NAs])/2)
      return(out)
    
    ## do not anlyse a cluster with NAs; for merged data it means such cluster was not present in all the datasets
    if(any(NAs) && skippNAs)
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
# Fit a logit GLM with inteactions + test contrasts with multcomp pckg
# Response is a count of 1s
# -----------------------------

fit_glm_interglht <- function(data, md, family = "binomial", formula, K, skippNAs = TRUE){
  
  ntot <- colSums(data[, md$shortname, drop = FALSE], na.rm = TRUE)
  contrast_names <- rownames(K)
  
  ### Fit the GLM
  fit <- lapply(1:nrow(data), function(i){
    # i = 1
    
    y <- data[i, md$shortname]
    NAs <- is.na(y)
    data_tmp <- data.frame(y = as.numeric(y), total = as.numeric(ntot), md)[!NAs, , drop = FALSE]
    
    mm <- model.matrix(formula, data = data_tmp)
    
    out <- matrix(NA, nrow = length(contrast_names), ncol = 2)
    colnames(out) <- c("coeff", "pval")
    rownames(out) <- contrast_names
    
    if(sum(y[!NAs] == 0) > length(y[!NAs])/2)
      return(out)
    
    ## do not anlyse a cluster with NAs; for merged data it means such cluster was not present in all the datasets
    if(any(NAs) && skippNAs)
      return(out)
    
    switch(family,

      binomial = {
        fit_tmp <- glm(formula, family = binomial(link = "logit"), data = data_tmp)
        summary(fit_tmp)
      },
      quasibinomial = {
        fit_tmp <- glm(formula, family = quasibinomial(link = "logit"), data = data_tmp)
        summary(fit_tmp)
      }
      
    )
    
    if(is.null(fit_tmp))
      return(out)
    
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
# Fit a logit GLMM with inteactions + test contrasts with multcomp pckg
# Response is a count of 1s
# -----------------------------


fit_glmer_interglht <- function(data, md, family = "binomial", formula, K, skippNAs = TRUE){
  
  ntot <- colSums(data[, md$shortname, drop = FALSE], na.rm = TRUE)
  contrast_names <- rownames(K)
  
  ### Fit the GLM
  fit <- lapply(1:nrow(data), function(i){
    # i = 1
    print(i)
    
    y <- data[i, md$shortname]
    NAs <- is.na(y)
    data_tmp <- data.frame(y = as.numeric(y), total = as.numeric(ntot), md)[!NAs, , drop = FALSE]
    
    out <- matrix(NA, nrow = length(contrast_names), ncol = 2)
    colnames(out) <- c("coeff", "pval")
    rownames(out) <- contrast_names
    
    if(sum(y[!NAs] == 0) > length(y[!NAs])/2)
      return(out)
      
    ## do not anlyse a cluster with NAs; for merged data it means such cluster was not present in all the datasets
    if(any(NAs) && skippNAs)
      return(out)
    
    switch(family, 

      binomial = {
        fit_tmp <- glmer(formula, weights = total, family = binomial, data = data_tmp)
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

# -----------------------------
# Fit a logit GLMM with inteactions + test contrasts with multcomp pckg
# Response is 0s or 1s
# -----------------------------



fit_glmer_interglht_01 <- function(data, md, family = "binomial", formula, K, skippNAs = TRUE){
  
  contrast_names <- rownames(K)
  
  ### Fit the GLM
  fit <- lapply(1:nrow(data), function(i){
    # i = 1
    print(i)
    
    y <- data[i, md$shortname]
    NAs <- is.na(y)

    data_tmp <- data.frame(y = as.numeric(y), md)[!NAs, , drop = FALSE]
    
    out <- matrix(NA, nrow = length(contrast_names), ncol = 2)
    colnames(out) <- c("coeff", "pval")
    rownames(out) <- contrast_names
    
    ## do not anlyse a cluster with NAs; for merged data it means such cluster was not present in all the datasets
    if(any(NAs) && skippNAs)
      return(out)
    
    switch(family, 
      binomial = {
        fit_tmp <- glmer(formula, family = binomial, data = data_tmp)
        summary(fit_tmp)
      }      
    )
    
    if(is.null(fit_tmp))
      return(out)
    
    ## fit contrasts one by one
    out <- t(apply(K, 1, function(k){
      # k = K[1, ]
      
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
























