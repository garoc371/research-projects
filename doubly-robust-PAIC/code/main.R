# Load packages
# package for data manipulation


if(!require("doSNOW")) {install.packages("doSNOW"); library(doSNOW)}
if(!require("boot")) {install.packages("boot"); library(boot)}
if(!require("copula")) {install.packages("copula"); library(copula)}
if(!require("dplyr")) {install.packages("dplyr"); library(dplyr)}
if(!require("here")) {install.packages("here"); library(here)}
if(!require("rstanarm")) {install.packages("rstanarm"); library(rstanarm)}


source("aipw_marginalization.R")
source("tmle_marginalization.R")

N_sim = 2000
N_sample <- 600
b_trt <- log(0.25)
b_EM <- -log(0.6)
b_1 <- -log(0.8)
sdX_AC <- 0.6
X1_var_target <- 0.3
X1_BC <- c(0.5)
event_rate <- 0.4
meanX_trial <- 0

param.combinations <- expand.grid(N_sample = N_sample, X1_BC = X1_BC)
pc <- param.combinations

scenarios <- nrow(pc) # number of simulation scenarios




# simulated patient-level (AC) and aggregate-level (BC) datasets for all scenarios
IPD.AC.all <- vector(mode="list", nrow(pc))
IPD.BC.all <- vector(mode = "list", nrow(pc))
ALD.BC.all <- vector(mode="list", nrow(pc))
# load data
for (i in 1:scenarios) {
  file.id <- paste0("N_sample_", pc$N_sample[i], "_X1_BC_", pc$X1_BC[i])
  load(paste0("pop_data/IPD_AC_", file.id, ".RData")) # load AC patient-level data
  load(paste0("pop_data/ALD_BC_", file.id, ".RData")) # load BC aggregate-level data
  load(paste0("pop_data/IPD_BC_", file.id, ".RData")) # load BC individual-level data
  IPD.AC.all[[i]] <- IPD.AC
  ALD.BC.all[[i]] <- ALD.BC
  IPD.BC.all[[i]] <- IPD.BC
}


n_cores <- 32
cl <- makeCluster(n_cores, type = "SOCK")
registerDoSNOW(cl)



for (i in 1:scenarios){
  IPD.AC <- IPD.AC.all[[i]]
  ALD.BC <- ALD.BC.all[[i]]
  IPD.BC <- IPD.BC.all[[i]]
  file.id <- paste0("N_sample_", pc$N_sample[i], "_X1_BC_", pc$X1_BC[i])
  
  maic.results <- foreach(j=1:N_sim, .errorhandling = "pass",
                          .packages=c("boot","dplyr")) %dopar% {
                            results <- maic.wrapper(IPD.AC[[j]], ALD.BC[[j]],
                                                    1000)
                            return(results)
                          }
  
  save(maic.results, file = paste0("results3/maic_", file.id, ".RData"))
  
  bayes_gcomp.results <- foreach(j=1:N_sim, .errorhandling = "pass",
                                 .packages=c("copula", "rstanarm","dplyr")) %dopar%{
                                   results <- gcomp.bayes.wrapper(IPD.AC[[j]], ALD.BC[[j]],
                                                                  IPD.BC[[j]],1000)
                                   
                                   return(results)
                                 }
  
  save(bayes_gcomp.results, file = paste0("results3/bayes_gcomp_", file.id, ".RData"))
  

  # AIPW fit each treatment arm separately
  aipw_model <- y ~ X1 + X2 + X3 + X4 + X5
  aipw.ml.results <- foreach(j=1:N_sim, .errorhandling = "pass",
                             .packages=c("copula","dplyr")) %dopar%{
                               results <- AIPW.ml.wrapper(IPD.AC[[j]], ALD.BC[[j]],
                                                          IPD.BC[[j]],1000, ymod = aipw_model, boot = T)

                               return(results)
                             }

  save(aipw.ml.results, file = paste0("results3/aipw_ml_", file.id, ".RData"))

  
  # weighted regression implementation of DR estimator
  weighted.dr.results <- foreach(j=1:N_sim, .errorhandling = "pass",
                                 .packages=c("copula","dplyr")) %dopar%{
                                   results <- weighted.gcomp.ml.wrapper(IPD.AC[[j]], ALD.BC[[j]],
                                                              IPD.BC[[j]],1000, boot = T)
                                   
                                   return(results)
                                 }
  
  save(weighted.dr.results, file = paste0("results3/weighted_dr_", file.id, ".RData"))
  
  # TMLE
  tmle_model <- y ~ trt*(X1 + X2 + X3 + X4 + X5)
  TMLE.results <- foreach(j=1:N_sim, .errorhandling = "pass",
                          .packages=c("copula","dplyr")) %dopar%{
                            results <- TMLE.ml.wrapper(IPD.AC[[j]], ALD.BC[[j]],
                                                       IPD.BC[[j]],1000, ymod = tmle_model, boot = T)
                            
                            return(results)
                          }
  
  save(TMLE.results, file = paste0("results3/tmle_", file.id, ".RData"))

}

stopCluster(cl)