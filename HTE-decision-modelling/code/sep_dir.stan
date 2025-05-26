// fit a monotonic spline model to each arm seperately

data {
  int<lower=1> N; // total number of observations
  array[N] int Y; // response variable
  int<lower=1> m; // number of basis
  // data for spline 
  // spline basis function matrices for control group
  matrix[N, m] S;
}


parameters {
  real Intercept; // temporary intercept for centered predictors
  // parameters for spline s(age_std, by = trt, bs = "cs")0
  // scalar for the spline coefficient
  real<lower = 0> alpha;
  simplex[m] gamma; // unscaled spline coefficient 
}

transformed parameters {
  // actual spline coefficients for control group
  vector[m] s_0;
  
  

  // compute actual spline coefficients
  s_0 = alpha * gamma;
  
}

model {
  // define yhat, and ZS in the model block
  vector[N] yhat;

  
  
  yhat = Intercept + S*s_0;

  alpha ~ exponential(1);
  Intercept ~ normal(0,3);
  
  gamma ~ dirichlet(rep_vector(.1,m));

  // AR(1) prior for spline coefficient
  
  for (i in 1:N) {
    Y[i] ~ bernoulli(inv_logit(yhat[i])); // the response
  }
}

generated quantities {
  
}