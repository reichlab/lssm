#include functions.stan

// parameter estimation/fitting
data {
  // number of observations
  int<lower=0> n;
  // length for each observation
  int<lower=0> p;
  // n obs each is px1
  vector[p] y[n];
  // number of lags
  int<lower=0> p_ar;
  
}

transformed data{
  // number of states = p_ar
  int<lower=0> m = p_ar;
  
  // observation intercept
  vector[p] d = rep_vector(0, p);
  
  // Z in 3.19
  matrix[p,m] Z = append_col([[1]],to_matrix(rep_row_vector(0, m-1)));
  
  // observation covariance
  matrix[p,p] H = [[0]];
  
  // state intercept
  vector[m] c = rep_vector(0, m);
  
  // state covariance
  matrix[m,1] R = append_row([[1]], to_matrix(rep_vector(0, m-1)));
  
}

parameters {
  
  vector[p_ar] phi;
  
  real <lower = 0> var_zeta;
  // expected value for the initial state
  vector[m] a1;
}

transformed parameters{
  
  // T in 3.20
  matrix[m,m] T = append_col(phi,
  append_row(
    diag_matrix(rep_vector(1, m-1)), to_matrix(rep_row_vector(0, m-1))));
  
  matrix[1,1] Q = [[var_zeta]];
  
  // variance for the initial state
  matrix[m,m] P1 = stationary_cov(T, quad_form_sym(Q, R '));
  
}

model {
  target += ssm_constant_lpdf (y| d, Z, H, c, T, R, Q, a1, P1);
  
}
