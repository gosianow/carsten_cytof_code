library(lme4) # for fitting mixed models
library(multcomp) # for contrasts glht()


# -----------------------------
### Run two-way ANOVA
# -----------------------------

# fit_two_way_anova <- function(data, md){
#   
#   pvals_anova <- t(apply(data[, md$shortname], 1, function(y){
#     # y <- data[1, md$shortname]
#     
#     NAs <- is.na(y)
#     
#     ## there must be at least 10 proportions greater than 0
#     if(sum(y[!NAs] > 0) < 10)
#       return(rep(NA, 2))
#     
#     data_tmp <- data.frame(y = as.numeric(y)[!NAs], md[!NAs, c("day", "response")])
#     
#     fit_tmp <- aov(y ~ day + response, data = data_tmp)
#     summ_tmp <- summary(fit_tmp)
#     
#     summ_tmp[[1]][1:2, "Pr(>F)"]
#     
#   }))
#   
#   movars <- c("day", "response")
#   colnames(pvals_anova) <- paste0("pval_", movars)
#   pvals_anova <- data.frame(pvals_anova)
#   
#   ## get adjusted p-values
#   adjp_anova <- data.frame(apply(pvals_anova, 2, p.adjust, method = "BH"))
#   colnames(adjp_anova) <- paste0("adjp_", movars)
#   
#   pvals_out <- data.frame(pvals_anova, adjp_anova)
#   
#   return(pvals_out)
#   
# }

# -----------------------------
### Fit a normal GLM - with response only, no day
# -----------------------------

fit_glm_norm_resp <- function(data, md){
  
  ### Fit the GLM
  fit <- lapply(1:nrow(data), function(i){
    # i = 1
    
    y <- data[i, md$shortname]
    NAs <- is.na(y)
    data_tmp <- data.frame(y = as.numeric(y)[!NAs], md[!NAs, c("response"), drop = FALSE])
    
    ## there must be at least 10 proportions greater than 0
    if(sum(y[!NAs] > 0) < 10){
      mm <- model.matrix(y ~ response, data = data_tmp)
      out <- matrix(NA, nrow = ncol(mm), ncol = 2)
      colnames(out) <- c("coeff", "pval")
      rownames(out) <- colnames(mm)
      return(out)
    }
    
    fit_tmp <- glm(y ~ response, data = data_tmp)
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
### Fit a normal GLM
# -----------------------------

fit_glm_norm <- function(data, md){
  
  ### Fit the GLM
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
    
    fit_tmp <- glm(y ~ response + day, data = data_tmp)
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
### Fit a normal GLM with inteactions
# -----------------------------

fit_glm_norm_inter <- function(data, md){
  
  ### Fit the GLM
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
    
    fit_tmp <- glm(y ~ response + day + response:day, data = data_tmp)
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
### Fit a normal GLM with inteactions + test contrasts with multcomp pckg
# -----------------------------

fit_glm_norm_interglht <- function(data, md){
  
  ### Fit the GLM
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
    
    fit_tmp <- glm(y ~ response + day + response:day, data = data_tmp)
    summ_tmp <- summary(fit_tmp)
    
    ## fit contrasts
    contr_tmp <- glht(fit_tmp, linfct = K)
    summ_tmp <- summary(contr_tmp)
    
    out <- cbind(summ_tmp$test$coefficients, summ_tmp$test$pvalues)
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
### Fit a logit GLM
# -----------------------------


fit_glm_logit <- function(data, md){
  
  ntot <- colSums(data[, md$shortname, drop = FALSE])
  
  pvals <- t(apply(data[, md$shortname], 1, function(y){
    # y <- data[1, md$shortname]
    
    NAs <- is.na(y)
    data_tmp <- data.frame(y = as.numeric(y)[!NAs], total = as.numeric(ntot), md[!NAs, c("day", "response")])
    
    if(sum(y[!NAs] > 0) < 10){
      mm <- model.matrix(y ~ response + day, data = data_tmp)
      out <- rep(NA, ncol(mm))
      names(out) <- colnames(mm)
      return(out)
    }
    
    fit_tmp <- glm(cbind(y, total-y) ~ response + day, family = binomial(link = "logit"), data = data_tmp)
    
    summ_tmp <- summary(fit_tmp)
    
    out <- as.numeric(summ_tmp$coefficients[, "Pr(>|z|)"])
    names(out) <- rownames(summ_tmp$coefficients)
    
    return(out)
    
  }))
  
  movars <- colnames(pvals)
  
  colnames(pvals) <- paste0("pval_", movars)
  pvals <- data.frame(pvals)
  
  ## get adjusted p-values
  adjp <- data.frame(apply(pvals, 2, p.adjust, method = "BH"))
  colnames(adjp) <- paste0("adjp_", movars)
  
  ## save the results
  pvals_out <- data.frame(pvals, adjp)
  
  return(pvals_out)
  
}



# -----------------------------
### Fit a logit GLMM 
# -----------------------------

fit_glmm_logit <- function(data, md){
  
  ntot <- colSums(data[, md$shortname, drop = FALSE])
  
  pvals <- t(apply(data[, md$shortname], 1, function(y){
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
    
    fit_tmp <- glmer(y/total ~ response + day + (1|patient_id), weights = total, family = binomial, data = data_tmp)
    
    summ_tmp <- summary(fit_tmp)
    
    out <- as.numeric(summ_tmp$coefficients[, "Pr(>|z|)"])
    names(out) <- rownames(summ_tmp$coefficients)
    
    return(out)
    
  }))
  
  movars <- colnames(pvals)
  
  colnames(pvals) <- paste0("pval_", movars)
  pvals <- data.frame(pvals)
  
  
  ## get adjusted p-values
  adjp <- data.frame(apply(pvals, 2, p.adjust, method = "BH"))
  colnames(adjp) <- paste0("adjp_", movars)
  
  ## save the results
  pvals_out <- data.frame(pvals, adjp)
  
  return(pvals_out)
  
}













