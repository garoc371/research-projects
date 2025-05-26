
# script for implementing tmle and AIPW estimator in limited data context
# the idea can be traced from Josey et al 2021, who suggested using EB weights
# in place of trial assignment odds.
# The implementation of AIPW should be straightforward as it only involves reweighting
# the residual by EB weights.
# The implementation of TMLE can be mroe involved as it requires the clever covariate
# in the ALD trial

maic <- function(X.EM) {
  # X.EM: centered S=1 effect modifiers
  # objective function to be minimized for standard method of moments MAIC
  # Objective function
  Q <- function(a1, X){
    return(sum(exp(X %*% a1)))
  }
  
  # Gradient function
  gradfn_mm <- function(a1, X){
    return(colSums(sweep(X, 1, exp(X %*% a1), "*")))
  }
  X.EM <- as.matrix(X.EM)
  N <- nrow(X.EM)
  K.EM <- ncol(X.EM)
  alpha <- rep(1,K.EM) # arbitrary starting point for the optimizer
  # objective function minimized using BFGS
  Q.min <- optim(fn=Q, X=X.EM, gr = gradfn_mm,par=alpha, method="BFGS")
  hat.alpha <- Q.min$par # finite solution is the logistic regression parameters
  log.hat.w <- rep(0, N)
  for (k in 1:K.EM) {
    log.hat.w <- log.hat.w + hat.alpha[k]*X.EM[,k]
  }
  hat.w <- exp(log.hat.w) # estimated weights
  return(list(hat.w = hat.w, hat.alpha = hat.alpha))
}

sim_target <- function(ipd, ald, n){
  # generate the parametric target population (correct marginals)
  # based on gaussian copula. 
  # the first two moments are sourced from ald
  # while the correlation structure comes from ipd
  
  x.EM <- dplyr::select(ipd, X1:X5) # AC effect modifiers 
  rho <- cor(x.EM, method = "spearman") 
  #  covariate simulation for comparator trial using copula package
  cop <- normalCopula(P2p(rho), 
                      dim=5, dispstr="un") # AC IPD pairwise correlations
  # sample covariates from approximate joint distribution using copula
  # we specify a sufficiently large sample to target the PATE
  mvd <- mvdc(copula=cop, margins=c("norm", "binom", "norm", "gamma", "binom"), # marginals
              # BC covariate means and standard deviations
              paramMargins=list(list(mean = ald$mean.X1, sd = ald$sd.X1),
                                list(size = 1, prob = ald$mean.X2),       
                                list(mean = ald$mean.X3, sd = ald$sd.X3),
                                list(shape = (ald$mean.X4/ald$sd.X4)^2, rate = ald$mean.X4/(ald$sd.X4^2)),
                                list(size =1, prob = ald$mean.X5)))
  # data frame of simulated covariates
  x_star <- as.data.frame(rMvdc(n, mvd))
  colnames(x_star) <- c("X1", "X2", "X3", "X4", "X5")
  
  x_star$trt <- NA; x_star$y <- NA
  return(x_star)
  
}

sim_target_misspecified <- function(ipd, ald, n){
  # generate the parametric target population
  # based on gaussian copula. 
  # deliberately mis-specify the marginal and use
  # Normal for all continuous covariates
  # the first two moments are sourced from ald
  # while the correlation structure comes from ipd
  
  x.EM <- dplyr::select(ipd, X1:X5) # AC effect modifiers 
  rho <- cor(x.EM, method = "spearman") 
  #  covariate simulation for comparator trial using copula package
  cop <- normalCopula(P2p(rho), 
                      dim=5, dispstr="un") # AC IPD pairwise correlations
  # sample covariates from approximate joint distribution using copula
  # we specify a sufficiently large sample to target the PATE
  mvd <- mvdc(copula=cop, margins=c("norm", "binom", "norm", "norm", "binom"), # marginals
              # BC covariate means and standard deviations
              paramMargins=list(list(mean = ald$mean.X1, sd = ald$sd.X1),
                                list(size = 1, prob = ald$mean.X2),       
                                list(mean = ald$mean.X3, sd = ald$sd.X3),
                                list(mean = ald$mean.X4, sd = ald$sd.X4),
                                list(size = 1, prob = ald$mean.X5)))
  # data frame of simulated covariates
  x_star <- as.data.frame(rMvdc(n, mvd))
  colnames(x_star) <- c("X1", "X2", "X3", "X4", "X5")
  
  x_star$trt <- NA; x_star$y <- NA
  
  return(x_star)
  
}

maic.wrapper <- function(data.AC, data.BC, resamples) { 
  # Inputs: data.AC - AC individual patient-level data; data.BC - BC aggregate-level data  
  # resamples - number of resamples for non-parametric bootstrap
  maic.boot <- function(data, indices) {
    dat <- data[indices,]
    x.EM <- dplyr::select(dat, X1:X5) # AC effect modifiers 
    theta <- data.BC[c("mean.X1", "mean.X2", "mean.X3", "mean.X4", "mean.X5")] %>% as.matrix()
    # center the AC effect modifiers on the BC means
    x.EM <- sweep(x.EM, 2, theta, "-")
    # MAIC weights estimated using standard method of moments
    hat.w <- maic(X.EM=x.EM)$hat.w # estimated weights
    # aess <- sum(hat.w)^2/sum(hat.w^2) # approximate effective sample size
    # fit weighted logistic regression model using glm
    outcome.fit <- glm(y~trt, family="quasibinomial", weights=hat.w, data=dat)
    # fitted treatment coefficient is marginal treatment effect for A vs. C
    hat.Delta.AC <- coef(outcome.fit)["trt"] 
    return(hat.Delta.AC)
  }
  # non-parametric bootstrap
  boot.object <- boot::boot(data=data.AC, statistic=maic.boot, R=resamples)
  # bootstrap mean of marginal A vs. C treatment effect estimate
  hat.Delta.AC <- mean(boot.object$t)
  # bootstrap variance of A vs. C treatment effect estimate   
  hat.var.Delta.AC <- var(boot.object$t)
  # B vs. C marginal treatment effect from reported event counts
  hat.Delta.BC <- with(data.BC, log(y.B.sum*(N.C-y.C.sum)/(y.C.sum*(N.B-y.B.sum))))
  # B vs. C marginal treatment effect variance using the delta method
  hat.var.Delta.BC <- with(data.BC, 1/y.C.sum+1/(N.C-y.C.sum)+1/y.B.sum+1/(N.B-y.B.sum))
  hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC # A vs. B
  hat.var.Delta.AB <- hat.var.Delta.AC + hat.var.Delta.BC  
  list(mean = hat.Delta.AB, var = hat.var.Delta.AB)
}

gcomp.bayes.wrapper <- function(ipd, ald, true_ald, ndraws) {
  
  # Helper function to compute effect for a given target population
  compute_effect <- function(target_pop, outcome.model, ald) {
    target_pop <- target_pop %>% dplyr::select(-y)
    data.trtA <- data.trtC <- target_pop
    data.trtA$trt <- 1
    data.trtC$trt <- 0
    
    y.star.A <- posterior_predict(outcome.model, newdata = data.trtA, draws = ndraws)
    y.star.C <- posterior_predict(outcome.model, newdata = data.trtC, draws = ndraws)
    
    hat.delta.AC <- qlogis(rowMeans(y.star.A)) - qlogis(rowMeans(y.star.C))
    hat.Delta.AC <- mean(hat.delta.AC)
    hat.var.Delta.AC <- var(hat.delta.AC)
    
    # Compute B vs. C effect
    hat.Delta.BC <- with(ald, log(y.B.sum * (N.C - y.C.sum) / (y.C.sum * (N.B - y.B.sum))))
    hat.var.Delta.BC <- with(ald, 1/y.C.sum + 1/(N.C - y.C.sum) + 1/y.B.sum + 1/(N.B - y.B.sum))
    
    # Compute final A vs. B effect
    delta_AB <- hat.Delta.AC - hat.Delta.BC
    var_AB <- hat.var.Delta.AC + hat.var.Delta.BC
    
    return(list(est = delta_AB, var = var_AB))
  }
  
  # Fit the outcome model
  outcome.model <- stan_glm(y ~ trt * (X1 + X2 + X3 + X4 + X5), 
                            family = binomial(link = "logit"), 
                            refresh = 0, data = ipd)
  
  # 1. True ALD population (oracle)
  oracle <- true_ald %>% dplyr::select(-id, -b0)
  effect_oracle <- compute_effect(oracle, outcome.model, ald = ald)
  
  # 2. ALD with correct marginals
  x_star_margin <- sim_target(ipd, ald, n = 5000)
  effect_margin <- compute_effect(x_star_margin, outcome.model, ald = ald)
  
  # 3. ALD with mis-specified marginals
  x_star_misspecified <- sim_target_misspecified(ipd, ald, n = 5000)
  effect_misspecified <- compute_effect(x_star_misspecified, outcome.model, ald = ald)
  
  # Return results
  return(list(
    est_oracle = effect_oracle$est, 
    var_oracle = effect_oracle$var,
    est_margin = effect_margin$est, 
    var_margin = effect_margin$var,
    est_misspecified = effect_misspecified$est, 
    var_misspecified = effect_misspecified$var
  ))
}


aipw.ml <- function(data_sample, data_target, ald, models) {
  
  if (is.na(models$m1_model)[1] || is.na(models$m0_model)[1]) {
    return(NA)
  }
  
  # assign sample indicator
  
  data_sample$S <- 1
  data_target$S <- 0
  
  # Combine datasets
  
  data_combined <- rbind(data_sample, data_target)
  n.sample <- nrow(data_sample)
  
  x.EM <- dplyr::select(data_sample, X1:X5) # AC effect modifiers 
  theta <- ald[c("mean.X1", "mean.X2", "mean.X3", "mean.X4", "mean.X5")] %>% as.matrix()
  x.EM <- sweep(x.EM, 2, theta, "-") # center the AC effect modifiers on the BC means

  # Calculate weights
  w <- maic(X.EM = x.EM)$hat.w
  w <- w/sum(w) * n.sample
  
  # Fit treatment model
  pi_model <- glm(trt ~ 1, family = binomial(), data = data_sample)
  pi_pred <- predict(pi_model, newdata = data_sample, type = "response")
  
  wt <- w / (data_sample$trt*pi_pred + (1-data_sample$trt)*(1-pi_pred))
  

  m1_pred <- predict(models$m1_model, newdata = data_combined, type = "response")
  m0_pred <- predict(models$m0_model, newdata = data_combined, type = "response")
  
  # AIPW estimates
  AIPW_1 <- sum(data_sample$trt * wt * (data_sample$y - m1_pred[data_combined$S == 1]))
  AIPW_0 <- sum((1 - data_sample$trt) * wt * (data_sample$y - m0_pred[data_combined$S == 1]))
  
  norm_1 <- sum(data_sample$trt * wt)^-1
  norm_0 <- sum((1 - data_sample$trt) * wt)^-1
  
  p1_aipw <- norm_1 * AIPW_1 + mean(m1_pred[data_combined$S == 0])
  p0_aipw <- norm_0 * AIPW_0 + mean(m0_pred[data_combined$S == 0])
  
  # Log odds ratio
  lor_aipw <- log(p1_aipw / (1 - p1_aipw)) - log(p0_aipw / (1 - p0_aipw))
  
  return(lor_aipw)
}


### AIPW with maximum-likelihood estimation and bootstrapping
AIPW.ml.wrapper <- function(ipd, ald, true_ald, resamples = 1000, ymod, boot = TRUE) {
  # Generate target populations
  target_data <- sim_target(ipd = ipd, ald = ald, n = 5000)
  target_data_misspecified <- sim_target_misspecified(ipd = ipd, ald = ald, n = 5000)
  oracle <- true_ald %>% dplyr::select(-id, -b0)
  
  # Function to fit outcome models
  fit_outcome_models <- function(data, ymod) {
    fit_model <- function(formula, subset_data) {
      tryCatch(
        glm(formula, family = binomial(), data = subset_data),
        warning = function(w) {
          if (grepl("did not converge|fitted probabilities 0 or 1", w$message)) 
            return(NA)
        else { 
          glm(formula, family = binomial(), data = subset_data)
        }
      })
    }
    list(
      m1_model = fit_model(ymod, data[data$trt == 1, ]),
      m0_model = fit_model(ymod, data[data$trt == 0, ])
    )
  }
  
  # Function to compute AIPW estimate
  compute_aipw <- function(data, target, ymodel) {
    aipw.ml(data, target, ymodel, ald = ald)
  }
  
  # Compute point estimates
  ymodel <- fit_outcome_models(ipd, ymod = ymod)
  results_oracle <- compute_aipw(ipd, oracle, ymodel)
  results_margin <- compute_aipw(ipd, target_data, ymodel)
  results_misspecified <- compute_aipw(ipd, target_data_misspecified, ymodel)
  
  # Calculate B vs C log-odds ratio and variance
  hat.Delta.BC <- with(ald, log(y.B.sum*(N.C-y.C.sum)/(y.C.sum*(N.B-y.B.sum))))
  hat.var.Delta.BC <- with(ald, 1/y.C.sum + 1/(N.C-y.C.sum) + 1/y.B.sum + 1/(N.B-y.B.sum))
  
  # Calculate final estimates
  delta_AB_oracle <- results_oracle - hat.Delta.BC
  delta_AB_margin <- results_margin - hat.Delta.BC
  delta_AB_misspecified <- results_misspecified - hat.Delta.BC
  
  if (!boot) {
    return(list(
      est_oracle = delta_AB_oracle,
      est_margin = delta_AB_margin,
      est_misspecified = delta_AB_misspecified
    ))
  }
  
  # Bootstrap
  boot_results <- replicate(resamples, {
    boot_sample <- ipd[sample(nrow(ipd), replace = TRUE), ]
    boot_ymodel <- fit_outcome_models(boot_sample, ymod = ymod)
    c(
      compute_aipw(boot_sample, oracle, boot_ymodel),
      compute_aipw(boot_sample, target_data, boot_ymodel),
      compute_aipw(boot_sample, target_data_misspecified, boot_ymodel)
    )
  })
  
  # Calculate variances
  var_AB_oracle <- var(boot_results[1,], na.rm = TRUE) + hat.var.Delta.BC
  var_AB_margin <- var(boot_results[2,], na.rm = TRUE) + hat.var.Delta.BC
  var_AB_misspecified <- var(boot_results[3,], na.rm = TRUE) + hat.var.Delta.BC
  
  return(list(
    est_oracle = delta_AB_oracle, var_oracle = var_AB_oracle,
    est_margin = delta_AB_margin, var_margin = var_AB_margin,
    est_misspecified = delta_AB_misspecified, var_misspecified = var_AB_misspecified
  ))
}

weighted.gcomp.ml.wrapper <- function(ipd, ald, true_ald, resamples, boot = T) {
  # Generate target populations
  target_data <- sim_target(ipd = ipd, ald = ald, n = 5000)
  target_data_misspecified <- sim_target_misspecified(ipd = ipd, ald = ald, n = 5000)
  oracle <- true_ald %>% dplyr::select(-id, -b0)
  
  # Function to fit outcome models (done once per bootstrap iteration)
  weighted_outcome_model <- function(data, ald) {
    x.EM <- dplyr::select(data, X1:X5)
    theta <- ald[c("mean.X1", "mean.X2", "mean.X3", "mean.X4", "mean.X5")] %>% as.matrix()
    x.EM <- sweep(x.EM, 2, theta, "-")
    hat.w <- maic(X.EM = x.EM)$hat.w
    
    model <- tryCatch({
      fit <- glm(y ~ trt*(X1 + X2 + X3 + X4 + X5), data=data, family=quasibinomial(), weights = hat.w)
    },  
    warning = function(w) {
      # Only return NA for non-convergence warnings
      if (grepl("glm.fit: algorithm did not converge", w$message)||
          grepl("fitted probabilities numerically 0 or 1 occurred", w$message)) {
        return(NA)
      }
      # For other warnings, print the message but return the fit
      message("Warning in GLM fit: ", w$message)
      return(glm(y ~ trt*(X1 + X2 + X3 + X4 + X5), data=data, family=binomial, weights = hat.w))
    })
    
    return(model)
  }
  
  # Function to calculate weighted G-computation estimate
  weighted.gcomp <- function(outcome_model, target_data) {
    if (is.na(outcome_model)[1]) {  # Check only the first element
      return(NA)
    }
    
    data.trtA <- data.trtC <- target_data
    data.trtA$trt <- 1
    data.trtC$trt <- 0
    
    hat.mu.A.i <- predict(outcome_model, type="response", newdata=data.trtA)
    hat.mu.C.i <- predict(outcome_model, type="response", newdata=data.trtC)
    hat.mu.A <- mean(hat.mu.A.i)
    hat.mu.C <- mean(hat.mu.C.i)
    
    log(hat.mu.A/(1-hat.mu.A)) - log(hat.mu.C/(1-hat.mu.C))
  }
  
  # Compute point estimates
  outcome_model <- weighted_outcome_model(ipd, ald)
  results_oracle <- weighted.gcomp(outcome_model, oracle)
  results_margin <- weighted.gcomp(outcome_model, target_data)
  results_misspecified <- weighted.gcomp(outcome_model, target_data_misspecified)
  
  # Calculate B vs C log-odds ratio and variance
  hat.Delta.BC <- with(ald, log(y.B.sum*(N.C-y.C.sum)/(y.C.sum*(N.B-y.B.sum))))
  hat.var.Delta.BC <- with(ald, 1/y.C.sum + 1/(N.C-y.C.sum) + 1/y.B.sum + 1/(N.B-y.B.sum))
  
  # Calculate final estimates
  delta_AB_oracle <- results_oracle - hat.Delta.BC
  delta_AB_margin <- results_margin - hat.Delta.BC
  delta_AB_misspecified <- results_misspecified - hat.Delta.BC
  
  if (!boot) {
    return(list(
      est_oracle = delta_AB_oracle,
      est_margin = delta_AB_margin,
      est_misspecified = delta_AB_misspecified
    ))
  }
  
  # Bootstrap
  boot_results <- replicate(resamples, {
    boot_sample <- ipd[sample(nrow(ipd), replace = TRUE), ]
    boot_model <- weighted_outcome_model(boot_sample, ald)
    c(
      weighted.gcomp(boot_model, oracle),
      weighted.gcomp(boot_model, target_data),
      weighted.gcomp(boot_model, target_data_misspecified)
    )
  })
  
  # Calculate variances
  var_AB_oracle <- var(boot_results[1,], na.rm = TRUE) + hat.var.Delta.BC
  var_AB_margin <- var(boot_results[2,], na.rm = TRUE) + hat.var.Delta.BC
  var_AB_misspecified <- var(boot_results[3,], na.rm = TRUE) + hat.var.Delta.BC
  
  return(list(
    est_oracle = delta_AB_oracle, var_oracle = var_AB_oracle,
    est_margin = delta_AB_margin, var_margin = var_AB_margin,
    est_misspecified = delta_AB_misspecified, var_misspecified = var_AB_misspecified
  ))
}






