#include functions.stan

data {
  // number of observations
  int<lower=0> n;
  // length for each observation
  int<lower=0> p;
  // n obs each is px1
  vector[p] y[n];
  
}

transformed data{
  // number of states = p_ar
  int<lower=0> m = 1;
  
  // observation intercept
  vector[p] d = rep_vector(0, p);
  
  matrix[p,m] Z = [[1]];
  
  // observation covariance
  matrix[p,p] H = [[0]];
  
  // state intercept
  vector[m] c = rep_vector(0, m);
  
  // state covariance
  matrix[m,1] R = [[1]];
  
}

parameters {
  real phi_1;
  real <lower = 0> var_zeta;
  // expected value for the initial state
  vector[m] a1;
}

transformed parameters{
  
  matrix[m,m] T = [[phi_1]];
  
  matrix[1,1] Q = [[var_zeta]];
  
  // variance for the initial state
  matrix[m,m] P1 = var_zeta * arima_stationary_cov(T, R);
  
}

model {
  target += ssm_constant_lpdf (y| d, Z, H, c, T, R, Q, a1, P1);
  
}

