rm(list=ls())

# Load packages
# package for data manipulation
if(!require("dplyr")) {install.packages("dplyr"); library(dplyr)}
# package to sample/simulate the covariates from a multivariate normal
if(!require("MASS")) {install.packages("MASS"); library(MASS)}
N_sim = 5000
N_AC <- c(100,200, 600)
N_BC <- 600
b_trt <- log(0.25)
b_EM <- -log(0.6)
b_1 <- -log(0.8)
sdX = 0.4  # sd for both trials in the Normal scenario
meanX_BC = 0.6
rhoX = 0.25
meanX_AC = c(0.45,0.375,0.25)
event_rate = 0.4

param.combinations <- expand.grid(N_AC=N_AC, meanX_AC=meanX_AC)
pc <- param.combinations

scenarios <- nrow(pc) # number of simulation scenarios

optim.function <- function(param, y_prob) {
  b_0 <- param
  fit <- sum(((1 / (1 + exp(-b_0))) - y_prob)^2)
  return(fit)
}


b_0 <- optim(par=0,fn=optim.function,y_prob=event_rate,
             method="Brent",lower=-2,upper=2)$par


# generate normally distributed covariates w/ \rho = 0 or 0.25
# fix ALD_X = 0.6, ALD_sd = 0.4, varying IPD mean & sd
# rho = 0, meanX = 0.45, sd = 0.4 -- 51% reduction in ESS;
#          meanX = 0.395, sd = 0.4  -- 71.5%
#          meanX = 0.35, sd = 0.4 -- 83.2% reduction in ESS
# rho = 0.25, meanX = 0.375, sd = 0.4 -- 54.8% reduction in ESS
#             meanX = 0.315 -- 70.5% reduction
#             meanX = 0.25 -- 82.7% reduction
sim_pop_norm <- function(N, b_0, b_trt, b_EM, b_1, meanX, sdX, rhoX = 0, event_rate){
  
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
  # treatment assignment (1: active; 0: control)
  trt <- c(rep(1,N_active),rep(0,N_control)) 
  # generate binary outcomes using logistic regression
  LP <- b_0 + b_1*(X$X1+X$X2+X$X3+X$X4+X$X5) + b_trt*trt + 
    b_EM*trt*(X$X1+X$X2+X$X3+X$X4+X$X5) # linear predictor
  yprob <- 1/ (1+exp(-LP)) # binary outcome probability
  y <- rbinom(n=N, size=1, prob=yprob) # binary outcome
  return(as.data.frame(cbind(X, trt, y)))
}

for (i in 1:scenarios) {
  # simulate IPD covariates and outcome for A vs. C trial (S=1)
  IPD.AC <- replicate(n=N_sim, expr=sim_pop_norm(pc$N_AC[i], b_0, b_trt, b_EM, 
                                             b_1, pc$meanX_AC[i], sdX, rhoX, 
                                             event_rate),
                      simplify=FALSE)
  # simulate IPD covariates and outcome for B vs. C trial (S=2)
  IPD.BC <- replicate(n=N_sim, expr=sim_pop_norm(N_BC, b_0, b_trt, b_EM, 
                                             b_1, meanX_BC, sdX, rhoX, 
                                             event_rate),
                      simplify=FALSE)
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
  file.id <- paste0("N_AC", pc$N_AC[i], "meanX_AC", pc$meanX_AC[i]) 
  save(IPD.AC, file=paste0("pop_data/IPD_AC_", file.id, ".RData"))
  save(IPD.BC, file=paste0("pop_data/IPD_BC_", file.id, ".RData"))
  save(ALD.BC, file=paste0("pop_data/ALD_BC_", file.id, ".RData"))  
} 
