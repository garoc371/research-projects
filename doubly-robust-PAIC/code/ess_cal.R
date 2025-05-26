# 0.39 corresponds to 49.9% redution in ESS
# 0.45 corresponds to 30.4% reducction in ESS
# 0.265 corresponds to 80.1% reduction in ESS

#
# dispersion = 5;
# 0.25 corresponds to 82.6% reduction in ESS
# 0.375 corresponds to 55.3% reduction in ESS
# 0.45 corresponds to 33.1% reduction in ESS


## based on the exploration below, we select target location for X1 at 0.25 and 0.5,
## this should roughly correspond to 55% and 76% ESS reduction
## the target variance for X1 have negligible impact for the three values we consider
## we set X1_var = 0.3 for this simulation study.


## generating trial and target sample from distinct source populations
if(!require("dplyr")) {install.packages("dplyr"); library(dplyr)}
# package to sample/simulate the covariates from a multivariate normal
if(!require("simstudy")) {install.packages("simstudy"); library(simstudy)}
source("gen_pop.R")
library(doParallel)

N_sim = 2000
N_AC <- c(100,200,600)
b_trt <- log(0.25)
b_EM <- -log(0.6)
b_1 <- -log(0.8)

event_rate = 0.4


optim.function <- function(param, y_prob) {
  b_0 <- param
  fit <- sum(((1 / (1 + exp(-b_0))) - y_prob)^2)
  return(fit)
}


b_0 <- optim(par=0,fn=optim.function,y_prob=event_rate,
             method="Brent",lower=-2,upper=2)$par



loss_trial_target <- function(meanX_trial, X1_bar_target, X1_var_target, N = 600){
  
  # Generate trial data
  IPD.trial <- gen_pop_trial(N, b_0 = b_0, b_trt = b_trt, b_EM = b_EM, b_1 = b_1, 
                             meanX = meanX_trial, sdX = 0.3, rhoX = 0.4, event_rate = 0.4)
  
  # Generate target data
  IPD.target <- sim_pop_mix(N, b_0 = b_0, X1_bar = X1_bar_target, X1_var = X1_var_target, event_rate = 0.4)
  
  # Calculate mean covariates for the target population
  ALD.target <- IPD.target %>%
    summarise(across(X1:X5, mean)) %>%
    as.matrix()
  
  # Calculate differences between trial covariates and target means
  EM <- IPD.trial %>%
    dplyr::select(X1:X5) %>%
    sweep(., 2, ALD.target, "-")
  
  # Calculate weights and effective sample size
  w <- maic(EM)$hat.w
  w <- w/sum(w)*N
  ess <- sum(w)^2/sum(w^2)
  loss <- 1 - ess/N
  
  return(loss)
}


n_cores <- 6 
cl <- makeCluster(n_cores)
registerDoParallel(cl)



calculate_mc_ess_reduction <- function(meanX_trial, X1_bar_target_grid, X1_var_target_grid,
                                       sample_sizes = 100, n_reps = 1000) {
  grid <- expand.grid(X1_bar_target = X1_bar_target_grid,
                         X1_var_target = X1_var_target_grid,
                         meanX_trial = meanX_trial)
  
  
  

  results <- foreach(i = 1:nrow(grid), .combine = rbind, .packages = c("dplyr")) %dopar% {
    set.seed(i)  # Ensure reproducibility
    losses <- replicate(n_reps, loss_trial_target(meanX_trial, grid$X1_bar_target[i], grid$X1_var_target[i]))
    data.frame(
      X1_bar_target = grid$X1_bar_target[i],
      X1_var_target = grid$X1_var_target[i],
      meanX_trial = grid$meanX_trial[i],
      mean_loss = mean(losses),
      min_loss = min(losses),
      max_loss = max(losses),
      q10_loss = quantile(losses, probs = 0.1),
      q90_loss = quantile(losses, probs = 0.9)
    )
  }
    
    stopCluster(cl)
    
    return(results)
}

# Usage
meanX_trial <- c(0, 0.1,0.2)
X1_bar_target_grid <- seq(0, 0.5, 0.05)
X1_var_target_grid <- c(0.2, 0.3, 0.4)
sample_sizes <- 100

mc_results <- calculate_mc_ess_reduction(meanX_trial, X1_bar_target_grid, X1_var_target_grid, sample_sizes)

# Plot results
library(ggplot2)
ggplot(mc_results, aes(x = X1_bar_target, y = mean_loss, color = factor(sample_size))) +
  geom_line() +
  geom_line(aes(y = q10_loss), linetype = "dashed") +
  geom_line(aes(y = q10_loss), linetype = "dashed") +
  theme_minimal() +
  labs(title = "Average ESS Reduction by Target Population Mean",
       x = "Target Population Mean (X1_bar)",
       y = "ESS Reduction",
       color = "Sample Size",
       fill = "Sample Size") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  facet_grid(X1_var_target~meanX_trial)

View(mc_results)

