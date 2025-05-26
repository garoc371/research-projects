# Load packages
# package for data manipulation


if(!require("doSNOW")) {install.packages("doSNOW"); library(doSNOW)}
if(!require("rstanarm")) {install.packages("rstanarm"); library(rstanarm)}
if(!require("polyapost")) {install.packages("polyapost"); library(polyapost)}
if(!require("boot")) {install.packages("boot"); library(boot)}
if(!require("copula")) {install.packages("copula"); library(copula)}
if(!require("dplyr")) {install.packages("dplyr"); library(dplyr)}

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


bucher <- function(data.AC, data.BC){
  m1 <- glm(y ~trt, data = data.AC, family = binomial)
  hat.Delta.AC = summary(m1)$coefficients["trt",1]
  hat.var.Delta.AC = (summary(m1)$coefficients["trt",2])^2
  hat.Delta.BC <- with(data.BC, log(y.B.sum*(N.C-y.C.sum)/(y.C.sum*(N.B-y.B.sum))))
  # B vs. C marginal treatment effect variance using the delta method
  hat.var.Delta.BC <- with(data.BC, 1/y.C.sum+1/(N.C-y.C.sum)+1/y.B.sum+1/(N.B-y.B.sum))
  hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC # A vs. B
  hat.var.Delta.AB <- hat.var.Delta.AC + hat.var.Delta.BC
  list(mean = hat.Delta.AB, var = hat.var.Delta.AB)
}


# Function to estimate MAIC weights 
maic <- function(X.EM) {
  # X.EM: centered S=1 effect modifiers
  # objective function to be minimized for standard method of moments MAIC
  Q <- function(alpha, X) {
    return(sum(exp(X %*% alpha)))
  }
  X.EM <- as.matrix(X.EM)
  N <- nrow(X.EM)
  K.EM <- ncol(X.EM)
  alpha <- rep(1,K.EM) # arbitrary starting point for the optimizer
  # objective function minimized using BFGS
  Q.min <- optim(fn=Q, X=X.EM, par=alpha, method="BFGS")
  hat.alpha <- Q.min$par # finite solution is the logistic regression parameters
  log.hat.w <- rep(0, N)
  for (k in 1:K.EM) {
    log.hat.w <- log.hat.w + hat.alpha[k]*X.EM[,k]
  }
  hat.w <- exp(log.hat.w) # estimated weights
  return(hat.w)
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
    hat.w <- maic(X.EM=x.EM) # estimated weights
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





CE_gcomp <- function(IPD, ALD, ndraws){
  # data.AC IPD for own trial; data.BC: ALD for opponent trial
  # ndraws: number of draws from dirichlet distribution
  # Rubin's BB is based on a Dir(0_n) prior. 
  # Dir(1_n) corresponds to each oberservations are trted equally likely
  # k is the parameter for the Dirichlet prior, 
  # under a classic Rubin's implementation, k = 0
  # alpha is the concentration parameter. Default alpha = 1
  
  x.EM <- dplyr::select(IPD, X1:X5) # AC effect modifiers 
  theta <- ALD[c("mean.X1", "mean.X2", "mean.X3", "mean.X4", "mean.X5")] %>% as.matrix()
  # center the AC effect modifiers on the BC means
  x.EM <- sweep(x.EM,2,theta,"-")
  # MAIC weights estimated using standard method of moments
  hat.w <- maic(X.EM=x.EM) # estimated weights
  vec_len <- length(hat.w)
  n_weights <- hat.w/sum(hat.w)*vec_len
  print(n_weights)
  outcome.fit <- stan_glm(y~trt*(X1+X2+X3+X4+X5), 
                          family=binomial(link = "logit"), 
                          data=IPD, refresh = 0)
  
  d0 = d1 = IPD
  d0$trt = 0; d1$trt = 1
  delta.AC = numeric(ndraws)
  
  for (i in seq_len(ndraws)){
    y0 = posterior_epred(outcome.fit, newdata = d0, type = "response", draws = 1)
    y1 = posterior_epred(outcome.fit, newdata = d1, type = "response", draws = 1)
    w = MCMCpack::rdirichlet(1, n_weights)
    
    y0_w = y0 %*% t(w)
    y1_w = y1 %*% t(w)
    
    delta.AC[i] = qlogis(y1_w) - qlogis(y0_w)
  }
  
  
  hat.Delta.AC = mean(delta.AC)
  hat.var.Delta.AC = var(delta.AC)
  hat.Delta.BC <- with(ALD, log(y.B.sum*(N.C-y.C.sum)/(y.C.sum*(N.B-y.B.sum))))
  # B vs. C marginal treatment effect variance using the delta method
  hat.var.Delta.BC <- with(ALD, 1/y.C.sum+1/(N.C-y.C.sum)+1/y.B.sum+1/(N.B-y.B.sum))
  hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC # A vs. B
  hat.var.Delta.AB <- hat.var.Delta.AC + hat.var.Delta.BC 
  
  gc()
  rm(outcome.fit, d0,d1)
  
  list(mean = hat.Delta.AB, var = hat.var.Delta.AB)
  
}

gcomp.bayes.wrapper <- function(ipd, ald, ndraws) {
  # Inputs: data.AC - AC individual patient-level data; data.BC - BC aggregate-level data;
  # n.chains, burnin, iters - MCMC info
  # N_star - size of simulated BC pseudo-population (high for small Monte Carlo error)
  # matrix of pairwise correlations between IPD covariates
  x.EM <- dplyr::select(ipd, X1:X5) # AC effect modifiers 
  rho <- cor(x.EM, method = "spearman") 
  #  covariate simulation for comparator trial using copula package
  cop <- normalCopula(P2p(rho), 
                      dim=5, dispstr="un") # AC IPD pairwise correlations
  # sample covariates from approximate joint distribution using copula
  n <- nrow(ipd)
  mvd <- mvdc(copula=cop, margins=c("norm", "norm", "norm", "norm", "norm"), # marginals
              # BC covariate means and standard deviations
              paramMargins=list(list(mean = ald$mean.X1, sd = ald$sd.X1),
                                list(mean = ald$mean.X2, sd = ald$sd.X2),       
                                list(mean = ald$mean.X3, sd = ald$sd.X3),
                                list(mean = ald$mean.X4, sd = ald$sd.X4),
                                list(mean = ald$mean.X5, sd = ald$sd.X5)))
  # data frame of simulated covariates
  x_star <- as.data.frame(rMvdc(n, mvd))
  colnames(x_star) <- c("X1", "X2", "X3", "X4", "X5")  
  # outcome logistic regression fitted to original data using MCMC  
  outcome.model <- stan_glm(y~trt*(X1 + X2 + X3 + X4 + X5), 
                            family=binomial(link = "logit"), 
                            refresh = 0, data = ipd) # run Stan model
  # counterfactual datasets
  data.trtA <- data.trtC <- x_star
  # intervene on received treatment while keeping set covariates fixed
  data.trtA$trt <- 1 # generate dataset where everyone receives treatment A
  data.trtC$trt <- 0 # generate dataset where all observations receive C  
  # draw counterfactual binary responses from posterior predictive distribution
  # matrix of posterior predictive draws under treatment A
  y.star.A <- posterior_predict(outcome.model, newdata=data.trtA, draws = ndraws)
  # matrix of posterior predictive draws under treatment C
  y.star.C <- posterior_predict(outcome.model, newdata=data.trtC, draws = ndraws)
  # compute marginal log-odds ratio for A vs. C for each MCMC sample
  # by transforming from probability to linear predictor scale  
  hat.delta.AC <- qlogis(rowMeans(y.star.A)) - qlogis(rowMeans(y.star.C)) 
  hat.Delta.AC <- mean(hat.delta.AC) # average over samples
  hat.var.Delta.AC <- var(hat.delta.AC) # sample variance
  # B vs. C from reported aggregate event counts, e.g. in contingency table
  hat.Delta.BC <- with(ald, log(y.B.sum*(N.C-y.C.sum)/(y.C.sum*(N.B-y.B.sum))))
  # B vs. C variance using the delta method 
  hat.var.Delta.BC <- with(ald, 1/y.C.sum+1/(N.C-y.C.sum)+1/y.B.sum+1/(N.B-y.B.sum))
  hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC # treatment effect for A vs. B
  hat.var.Delta.AB <- hat.var.Delta.AC + hat.var.Delta.BC # A vs. B variance
  list(mean = hat.Delta.AB, var = hat.var.Delta.AB)
}

mis.gcomp.bayes.wrapper <- function(ipd, ald, ndraws) {
  # Inputs: data.AC - AC individual patient-level data; data.BC - BC aggregate-level data;
  # n.chains, burnin, iters - MCMC info
  # N_star - size of simulated BC pseudo-population (high for small Monte Carlo error)
  # matrix of pairwise correlations between IPD covariates
  x.EM <- dplyr::select(ipd, X1:X5) # AC effect modifiers 
  rho <- cor(x.EM, method = "spearman") 
  #  covariate simulation for comparator trial using copula package
  cop <- normalCopula(P2p(rho), 
                      dim=5, dispstr="un") # AC IPD pairwise correlations
  # sample covariates from approximate joint distribution using copula
  n <- nrow(ipd)
  mvd <- mvdc(copula=cop, margins=c("gamma", "gamma", "gamma", "gamma", "gamma"), # marginals
              # BC covariate means and standard deviations
              paramMargins=list(list(shape = (ald$mean.X1/ald$sd.X1)^2, rate = ald$mean.X1/(ald$sd.X1^2)),
                                list(shape = (ald$mean.X2/ald$sd.X2)^2, rate = ald$mean.X2/(ald$sd.X2^2)),       
                                list(shape = (ald$mean.X3/ald$sd.X3)^2, rate = ald$mean.X3/(ald$sd.X3^2)),
                                list(shape = (ald$mean.X4/ald$sd.X4)^2, rate = ald$mean.X4/(ald$sd.X4^2)),
                                list(shape = (ald$mean.X5/ald$sd.X5)^2, rate = ald$mean.X5/(ald$sd.X5^2))))
  # data frame of simulated covariates
  x_star <- as.data.frame(rMvdc(n, mvd))
  colnames(x_star) <- c("X1", "X2", "X3", "X4", "X5")  
  # outcome logistic regression fitted to original data using MCMC  
  outcome.model <- stan_glm(y~trt*(X1 + X2 + X3 + X4 + X5), 
                            family=binomial(link = "logit"), 
                            refresh = 0, data = ipd) # run Stan model
  # counterfactual datasets
  data.trtA <- data.trtC <- x_star
  # intervene on received treatment while keeping set covariates fixed
  data.trtA$trt <- 1 # generate dataset where everyone receives treatment A
  data.trtC$trt <- 0 # generate dataset where all observations receive C  
  # draw counterfactual binary responses from posterior predictive distribution
  # matrix of posterior predictive draws under treatment A
  y.star.A <- posterior_predict(outcome.model, newdata=data.trtA, draws = ndraws)
  # matrix of posterior predictive draws under treatment C
  y.star.C <- posterior_predict(outcome.model, newdata=data.trtC, draws = ndraws)
  # compute marginal log-odds ratio for A vs. C for each MCMC sample
  # by transforming from probability to linear predictor scale  
  hat.delta.AC <- qlogis(rowMeans(y.star.A)) - qlogis(rowMeans(y.star.C)) 
  hat.Delta.AC <- mean(hat.delta.AC) # average over samples
  hat.var.Delta.AC <- var(hat.delta.AC) # sample variance
  # B vs. C from reported aggregate event counts, e.g. in contingency table
  hat.Delta.BC <- with(ald, log(y.B.sum*(N.C-y.C.sum)/(y.C.sum*(N.B-y.B.sum))))
  # B vs. C variance using the delta method 
  hat.var.Delta.BC <- with(ald, 1/y.C.sum+1/(N.C-y.C.sum)+1/y.B.sum+1/(N.B-y.B.sum))
  hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC # treatment effect for A vs. B
  hat.var.Delta.AB <- hat.var.Delta.AC + hat.var.Delta.BC # A vs. B variance
  list(mean = hat.Delta.AB, var = hat.var.Delta.AB)
}


# simulated patient-level (AC) and aggregate-level (BC) datasets for all scenarios
IPD.AC.all <- vector(mode="list", nrow(pc))
ALD.BC.all <- vector(mode="list", nrow(pc))
# load data
for (i in 1:scenarios) {
  file.id <- paste0("N_AC", pc$N_AC[i], "meanX_AC", pc$meanX_AC[i]) 
  load(paste0("pop_data/normal/IPD_AC_", file.id, ".RData")) # load AC patient-level data
  load(paste0("pop_data/normal/ALD_BC_", file.id, ".RData")) # load BC aggregate-level data
  IPD.AC.all[[i]] <- IPD.AC
  ALD.BC.all[[i]] <- ALD.BC
}

n_cores <- 32
cl <- makeCluster(n_cores, type = "SOCK")
registerDoSNOW(cl)

for (i in 1:scenarios){
  IPD.AC <- IPD.AC.all[[i]]
  ALD.BC <- ALD.BC.all[[i]]
  file.id <- paste0("N_AC", pc$N_AC[i], "meanX_AC", pc$meanX_AC[i])
  
  bucher.results <- foreach(j=1:N_sim, .errorhandling = "pass") %dopar% {
    results <- bucher(IPD.AC[[j]], ALD.BC[[j]])
    
    return(results)
  }
  
  save(bucher.results, file = paste0("results/normal3/bucher_", file.id, ".RData"))
  
  
  
  maic.results <- foreach(j=1:N_sim, .errorhandling = "pass",
                          .packages=c("boot","dplyr")) %dopar% {
                            results <- maic.wrapper(IPD.AC[[j]], ALD.BC[[j]],
                                                    1000)
                            return(results)
                          }
  
  save(maic.results, file = paste0("results/normal3/maic_", file.id, ".RData"))
  
  CE.results <- foreach(j=1:N_sim, .errorhandling = "pass",
                        .packages=c("rstanarm","dplyr")) %dopar%{
                          
                          results <- CE_gcomp(IPD.AC[[j]], ALD.BC[[j]],
                                              1000)
                          
                          return(results)
                        }
  save(CE.results, file = paste0("results/normal3/ce_", file.id, ".RData"))
  
  bayes_gcomp.results <- foreach(j=1:N_sim, .errorhandling = "pass",
                                 .packages=c("copula", "rstanarm","dplyr")) %dopar%{
                                   results <- gcomp.bayes.wrapper(IPD.AC[[j]], ALD.BC[[j]],
                                                                  1000)
                                   
                                   return(results)
                                 }
  
  save(bayes_gcomp.results, file = paste0("results/normal3/bayes_gcomp_", file.id, ".RData"))
  
  mis_bayes_gcomp.results <- foreach(j=1:N_sim, .errorhandling = "pass",
                                     .packages=c("copula", "rstanarm","dplyr")) %dopar%{
                                       results <- mis.gcomp.bayes.wrapper(IPD.AC[[j]], ALD.BC[[j]],
                                                                          1000)
                                       
                                       return(results)
                                     }
  
  save(mis_bayes_gcomp.results, file = paste0("results/normal3/mis_bayes_gcomp_", file.id, ".RData"))
}

stopCluster(cl)
