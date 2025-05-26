# TMLE estimation function
tmle_estimate <- function(data_sample, data_target, ald, y_model) {
  
  # assign sample indicator
  data_sample$S <- 1
  data_target$S <- 0
  
  # Combine datasets
  data_combined <- rbind(data_sample, data_target)
  n.dat <- nrow(data_combined)
  n.sample <- nrow(data_sample)
  n.target <- nrow(data_target)
  
  x.EM <- dplyr::select(data_sample, X1:X5)
  theta <- ald[c("mean.X1", "mean.X2", "mean.X3", "mean.X4", "mean.X5")] %>% as.matrix()
  x.EM <- sweep(x.EM, 2, theta, "-")
  
  # Calculate weights
  mom.rslt <- maic(X.EM = x.EM)
  w.s1 <- mom.rslt$hat.w
  w.s1 <- w.s1/sum(w.s1)*n.sample
  x.s0 <- dplyr::select(data_target, X1:X5)
  w.s0 <- exp(as.matrix(x.s0) %*% mom.rslt$hat.alpha)
  w.s0 <- w.s0/sum(w.s0)*n.target
  w <- c(w.s1, w.s0)
  
  # Calculate P(S=1)
  p_s1 <- mean(data_combined$S == 1)
  
  # Fit treatment model
  pi_model <- glm(trt ~ 1, family = binomial(), data = data_sample)
  pi_pred <- predict(pi_model, newdata = data_combined, type = "response")
  
  # Calculate clever covariates
  g0w <- (1 - pi_pred) / w
  g1w <- pi_pred / w
  
  h0w <- ((1 - data_combined$trt) * I(data_combined$S == 1)) / g0w
  h1w <- (data_combined$trt * I(data_combined$S == 1)) / g1w 
  
  data_new0 <- data_new1 <- data_combined
  data_new0$trt <- 0
  data_new1$trt <- 1
  
  # Initial prediction using the pre-fitted model
  q <- cbind(predict(y_model, type = "link", newdata = data_combined),
             predict(y_model, type = "link", newdata = data_new0),
             predict(y_model, type = "link", newdata = data_new1))
  
  # Updating step with warning handling
  update_model <- tryCatch({
    glm(y ~ -1 + offset(q[, 1]) + h0w + h1w, 
        family = binomial(), data = data_combined,
        subset = S == 1)
  }, warning = function(w) {
    if (grepl("glm.fit: algorithm did not converge", w$message) ||
        grepl("fitted probabilities numerically 0 or 1 occurred", w$message)) {
      return(NA)
    }
    glm(y ~ -1 + offset(q[, 1]) + h0w + h1w, 
        family = binomial(), data = data_combined,
        subset = S == 1)
  })
  
  if (is.na(update_model)[1]) return(NA)
  
  epsilon <- coef(update_model)
  
  # Update initial prediction
  q1 <- plogis(q + cbind(
    epsilon[1] * h0w + epsilon[2] * h1w,
    epsilon[1] / g0w,
    epsilon[2] / g1w
  ))
  
  # Calculate TMLE estimates
  p1_tmle <- mean(q1[, 3][data_combined$S == 0])
  p0_tmle <- mean(q1[, 2][data_combined$S == 0])
  
  # Log odds ratio
  lor_tmle <- log(p1_tmle / (1 - p1_tmle)) - log(p0_tmle / (1 - p0_tmle))
  
  return(lor_tmle)
}

# Wrapper function for TMLE with bootstrapping
TMLE.ml.wrapper <- function(ipd, ald, true_ald, resamples = 1000, ymod, boot = T) {
  # Generate target populations
  target_data <- sim_target(ipd = ipd, ald = ald, n = 5000)
  target_data_misspecified <- sim_target_misspecified(ipd = ipd, ald = ald, n = 5000)
  oracle <- true_ald %>% dplyr::select(-id, -b0)
  
  # Function to fit outcome model with warning handling
  fit_y_model <- function(data) {
    tryCatch({
      glm(as.formula(ymod), family = binomial(), data = data)
    }, warning = function(w) {
      if (grepl("glm.fit: algorithm did not converge", w$message) ||
          grepl("fitted probabilities numerically 0 or 1 occurred", w$message)) {
        return(NA)
      }
      glm(as.formula(ymod), family = binomial(), data = data)
    })
  }
  
  # Calculate point estimates
  y_model <- fit_y_model(ipd)
  if (is.na(y_model)[1]) {
    point_estimates <- list(
      est_oracle = NA,
      est_margin = NA,
      est_misspecified = NA
    )
  } else {
    point_estimates <- list(
      est_oracle = tmle_estimate(ipd, oracle, ald, y_model),
      est_margin = tmle_estimate(ipd, target_data, ald, y_model),
      est_misspecified = tmle_estimate(ipd, target_data_misspecified, ald, y_model)
    )
  }
  
  # Calculate B vs. C log odds ratio
  hat.Delta.BC <- with(ald, log(y.B.sum*(N.C-y.C.sum)/(y.C.sum*(N.B-y.B.sum))))
  hat.var.Delta.BC <- with(ald, 1/y.C.sum+1/(N.C-y.C.sum)+1/y.B.sum+1/(N.B-y.B.sum))
  
  # Adjust point estimates
  for (name in names(point_estimates)) {
    point_estimates[[name]] <- point_estimates[[name]] - hat.Delta.BC
  }
  
  if (!boot) {
    return(point_estimates)
  }
  
  # Bootstrap
  results <- list(
    oracle = numeric(resamples),
    margin = numeric(resamples),
    misspecified = numeric(resamples)
  )
  
  cat("Starting bootstrap...\n")
  for (i in 1:resamples) {
    boot_sample <- ipd[sample(nrow(ipd), replace = TRUE), ]
    y_model <- fit_y_model(boot_sample)
    
    if (is.na(y_model)[1]) {
      results$oracle[i] <- results$margin[i] <- results$misspecified[i] <- NA
    } else {
      results$oracle[i] <- tmle_estimate(boot_sample, oracle, ald, y_model)
      results$margin[i] <- tmle_estimate(boot_sample, target_data, ald, y_model)
      results$misspecified[i] <- tmle_estimate(boot_sample, target_data_misspecified, ald, y_model)
    }
  }
  
  # Function to safely calculate variance, excluding infinite values
  safe_var <- function(x) {
    x_finite <- x[is.finite(x)]
    if (length(x_finite) == 0) return(NA)
    var(x_finite, na.rm = TRUE)
  }
  
  # Calculate variances
  variances <- lapply(results, safe_var)
  
  # Prepare final results
  final_results <- list()
  for (name in c("oracle", "margin", "misspecified")) {
    est_name <- paste0("est_", name)
    var_name <- paste0("var_", name)
    final_results[[est_name]] <- point_estimates[[est_name]]
    final_results[[var_name]] <- variances[[name]] + hat.var.Delta.BC
  }
  
  return(final_results)
}