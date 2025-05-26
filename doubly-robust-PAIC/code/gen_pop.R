## generating trial and target sample from distinct source populations
if(!require("dplyr")) {install.packages("dplyr"); library(dplyr)}
# package to sample/simulate the covariates from a multivariate normal
if(!require("simstudy")) {install.packages("simstudy"); library(simstudy)}

if(!require("purrr")) {install.packages("purrr"); library(purrr)}


# the location of the target X1_BC is set to correspond to roughly 55% and 75%
# average ESS reduction

N_sim = 2000
N_sample <- c(100,200,600)
b_trt <- log(0.25)
b_EM <- -log(0.6)
b_1 <- -log(0.8)
sdX_AC <- 0.6
X1_var_target <- 0.3
X1_BC <- c(0.2, 0.5) 
event_rate <- 0.4
meanX_trial <- 0



optim.function <- function(param, y_prob) {
  b_0 <- param
  fit <- sum(((1 / (1 + exp(-b_0))) - y_prob)^2)
  return(fit)
}


b_0 <- optim(par=0,fn=optim.function,y_prob=event_rate,
             method="Brent",lower=-2,upper=2)$par


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

calculate_ess_loss <- function(IPD.target, IPD.trial) {
  ALD.target <- IPD.target %>%
    summarise(across(X1:X5, mean)) %>%
    as.matrix()
  
  # Calculate differences between trial covariates and target means
  EM <- IPD.trial %>%
    dplyr::select(X1:X5) %>%
    sweep(., 2, ALD.target, "-")
  
  # Calculate weights and effective sample size
  w <- maic(EM)$hat.w
  w <- w/sum(w)*nrow(IPD.trial)
  ess <- sum(w)^2/sum(w^2)
  loss <- 1 - ess/nrow(IPD.trial)
  
  return(loss)
}


gen_pop_trial <- function(N, b_0, b_trt, b_EM, b_1, meanX, sdX, rhoX = 0, event_rate){
  
  # 5 baseline covariates
  rho <- matrix(rhoX, nrow=5, ncol=5) # set correlation matrix
  diag(rho) <- rep(1, 5)
  N_active <- round(N*0.5) # number of patients under active treatment
  N_control <- N - N_active # number of patients under control
  sd.vec <-rep(sdX, 5) # vector of standard deviations
  cor2cov <- function(R, S) {
    # function to compute covariance matrix from correlation matrix R and vector
    # of standard deviations S. covariance matrix required as input for mvrnorm
    sweep(sweep(R, 1, S, "*"),2,S,"*")
  }
  cov.mat <- cor2cov(rho, sd.vec) # covariance matrix
  # simulate correlated continuous covariates using multivariate normal
  # patients under active treatment
  X_active <- as.data.frame(MASS::mvrnorm(n=N_active,mu=rep(meanX,5), Sigma=cov.mat))
  # patients under control treatment
  X_control <- as.data.frame(MASS::mvrnorm(n=N_control,mu=rep(meanX,5), Sigma=cov.mat))  
  # all patients
  X <- rbind(X_active, X_control)
  colnames(X) <- c("X1","X2","X3","X4", "X5") 
  
  # Transform X2 to binary using probit
  X$X2 <- as.integer(X$X2 > 0)
  
  # Transform X5 to binary using inverse logit (sigmoid)
  X$X5 <- as.integer(1 / (1 + exp(-X$X5)) > 0.5)
  
  
  # treatment assignment (1: active; 0: control)
  trt <- c(rep(1,N_active),rep(0,N_control)) 
  # generate binary outcomes using logistic regression
  LP <- b_0 + b_1*(X$X1+X$X2+X$X3+X$X4+X$X5) + b_trt*trt + 
    b_EM*trt*(X$X1+X$X2+X$X3+X$X4+X$X5) # linear predictor
  yprob <- 1/ (1+exp(-LP)) # binary outcome probability
  y <- rbinom(n=N, size=1, prob=yprob) # binary outcome
  return(as.data.frame(cbind(X, trt, y)))
}

sim_pop_mix <- function(n, b_0, X1_bar, X1_var, event_rate){
  
  
  def1 <- defData(varname = "X1", dist = "normal", formula = X1_bar, variance = X1_var)
  def1 <- defData(def1, varname = "X2", dist = "binary", 
                  formula = "-0.5 + X1 * 0.5", link = "logit")
  def1 <- defData(def1, varname = "X3", dist = "normal",
                  formula = "-0.5*X2", variance = 1)
  def1 <- defData(def1, varname = "X4", dist = "gamma", 
                  formula = "abs(0.2*X1 + 0.3*X3)", variance = 0.5)
  def1 <- defData(def1, varname = "X5", dist = "binary",
                  formula = "0.1*X1 - 0.1*X2 + 0.05*X3", link = "logit")
  
  def1 <- defData(def1, varname = "trt", dist = "trtAssign",
                  formula = "1;1")
  
  def1 <- defData(def1,varname = "b0", dist = "nonrandom", formula = b_0)
  
  def1 <- defData(def1, varname = "y", dist = "binary",
                  formula = "b0 + trt*(log(0.25) - log(0.6)*(X1 + X2 + X3 + X4 + X5)) - 
                    log(0.8)*(X1 + X2 + X3 + X4 + X5) ", link = "logit")
  
  d <- genData(n, def1)
  
  return(d)
}



# Generate data for all scenarios
scenarios <- expand.grid(N_sample = N_sample, X1_BC = X1_BC)

for (i in 1:nrow(scenarios)) {
  
  file_id <- paste0("N_sample_", scenarios$N_sample[i], "_X1_BC_", scenarios$X1_BC[i])
  
  IPD.AC <- replicate(N_sim, gen_pop_trial(scenarios$N_sample[i], b_0, b_trt, b_EM, b_1, 
                                               meanX = meanX_trial, sdX = sdX_AC, 
                                               rhoX = 0.4, event_rate = event_rate),
                      simplify = FALSE)
  IPD.BC <- replicate(N_sim, sim_pop_mix(scenarios$N_sample[i], b_0 = b_0, X1_bar = scenarios$X1_BC[i], 
                                             X1_var = X1_var_target, event_rate = event_rate),
                      simplify = FALSE)
  
  # Calculate ESS reduction for each pair
  loss_values <- map2_dbl(IPD.AC, IPD.BC, calculate_ess_loss)
  
  
  # Summarize BC IPD as ALD
  ALD.BC <- lapply(1:N_sim, function(j) {
    as.data.frame(cbind(
      # aggregate the data for the BC trial 
      summarise(IPD.BC[[j]], mean.X1=mean(X1), mean.X2=mean(X2), mean.X3=mean(X3),
                mean.X4=mean(X4), mean.X5=mean(X5), sd.X1=sd(X1), sd.X2=sd(X2), sd.X3=sd(X3), sd.X4=sd(X4), sd.X5=sd(X5)),
      # summarize the outcomes for the BC trial (treatment B)
      filter(IPD.BC[[j]], trt == 1) %>%
        summarise(y.B.sum=sum(y), y.B.bar = mean(y), N.B = n()),
      # summarize the outcomes for the BC trial (treatment C)
      filter(IPD.BC[[j]], trt == 0) %>%
        summarise(y.C.sum=sum(y), y.C.bar = mean(y), N.C = n())))    
  } )
  
  
  
  # Save the data
  save(IPD.AC, file = paste0("pop_data/IPD_AC_", file_id, ".RData"))
  save(IPD.BC, file = paste0("pop_data/IPD_BC_", file_id, ".RData"))
  save(ALD.BC, file = paste0("pop_data/ALD_BC_", file_id, ".RData"))
  
  # Save loss summary for this scenario
  loss_summary <- data.frame(
    N_sample = scenarios$N_sample[i],
    X1_bar_target = scenarios$X1_BC[i],
    mean_loss = mean(loss_values),
    median_loss = median(loss_values),
    min_loss = min(loss_values),
    max_loss = max(loss_values)
  )
  
  save(loss_summary, file = paste0("pop_data/ESS_loss/scenario_", file_id, ".RData"))
  
  cat("Completed scenario:", i, "of", nrow(scenarios), "\n")
  
}

