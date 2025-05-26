if(!require("dplyr")) {install.packages("dplyr"); library(dplyr)}
# package to sample/simulate the covariates from a multivariate normal
if(!require("simstudy")) {install.packages("simstudy"); library(simstudy)}
N_sim = 5000
N_AC <- c(100,200,600)
N_BC <- 600
b_trt <- log(0.25)
b_EM <- -log(0.6)
b_1 <- -log(0.8)
X1_AC = 0
var1_AC = 0.2
var1_BC = 0.4
X1_BC = c(0.275,0.4,0.6)
event_rate = 0.4

param.combinations <- expand.grid(N_AC=N_AC, X1_BC=X1_BC)
pc <- param.combinations

scenarios <- nrow(pc) # number of simulation scenarios

optim.function <- function(param, y_prob) {
  b_0 <- param
  fit <- sum(((1 / (1 + exp(-b_0))) - y_prob)^2)
  return(fit)
}


b_0 <- optim(par=0,fn=optim.function,y_prob=event_rate,
             method="Brent",lower=-2,upper=2)$par


# fix IPD_X1 = 0, IPD_var = 0.2, ALD_var = 0.4
# varying ALD_X1 = 0.4/0.5/0.6 to get ESS reduction of 55%, 70% and 84.7%
# generate binary data with fixed event rate
# n: size of the trial
# event_rate: baseline event rate in the control group
# X1_bar: standardized mean of X1, used to generate other dependent covariates
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

for (i in 1:scenarios) {
  # simulate IPD covariates and outcome for A vs. C trial (S=1)
  IPD.AC <- replicate(n=N_sim, expr=sim_pop_mix(pc$N_AC[i], event_rate = event_rate,
                                                b_0 = b_0, X1_bar = X1_AC, X1_var = var1_AC),
                      simplify=FALSE)
  # simulate IPD covariates and outcome for B vs. C trial (S=2)
  IPD.BC <- replicate(n=N_sim, expr=sim_pop_mix(N_BC, event_rate = event_rate,
                                                 b_0 = b_0, X1_bar=pc$X1_BC[i], X1_var=var1_BC),
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
  file.id <- paste0("N_AC", pc$N_AC[i], "X1_BC", pc$X1_BC[i]) 
  save(IPD.AC, file=paste0("pop_data/non_norm/IPD_AC_", file.id, ".RData"))
  save(IPD.BC, file=paste0("pop_data/non_norm/IPD_BC_", file.id, ".RData"))
  save(ALD.BC, file=paste0("pop_data/non_norm/ALD_BC_", file.id, ".RData"))  
} 
