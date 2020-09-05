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
  // number of previous noise terms
  int<lower=0> q_ma;
  
  int <lower=1> steps_ahead;
  
  // parameters
  vector[p_ar] phi;
  
  // lag coefficients
  vector[q_ma] theta;
  
  real <lower = 0> var_zeta;
  
  // r = max(p_ar, q_ma+1)
  int <lower=0> r;
  
  // number of states
  // m = r
  int<lower=0> m;
  
  // expected value for the initial state
  vector[m] a1;
  
}

transformed data{
  
  // observation intercept
  vector[p] d = rep_vector(0, p);
  
  // Z in 3.19
  matrix[p,m] Z = append_col([[1]],to_matrix(rep_row_vector(0, m-1)));
  
  // observation covariance
  matrix[p,p] H = [[0]];
  
  // state intercept
  vector[m] c = rep_vector(0, m);
  
  // transformed parameters
  
  vector [r-1] dummy_theta = append_row(theta, rep_vector(0, r-1-q_ma));
  
  vector [r] dummy_phi = append_row(phi,rep_vector(0, r-p_ar));
  
  matrix[1,1] Q = [[var_zeta]];
  
  // state covariance R in 3.20
  matrix[m,1] R = append_row([[1]], to_matrix(dummy_theta));
  
  // T in 3.20
  matrix[m,m] T = append_col(dummy_phi,
  append_row(
    diag_matrix(rep_vector(1, m-1)), to_matrix(rep_row_vector(0, m-1))));
  
  
  // variance for the initial state
  matrix[m,m] P1 = var_zeta * stationary_cov(T, quad_form_sym(Q, R '));

  
}

parameters {
  
  
}

generated quantities{
  matrix[p,p+1] results;
  results = predict (y, d, Z, H, c, T, R, Q, a1, P1,steps_ahead);
  
}


