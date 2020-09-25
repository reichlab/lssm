#include functions.stan

data {
  // number of observations
  int<lower=0> n;
  // length for each observation
  int<lower=0> p;
  // n obs each is px1
  vector[p] y[n];
  // number of lags in non-seasonal part
  int<lower=0> p_ar;
  // number of previous noise terms in non-seasonal part
  int<lower=0> q_ma;
  // number of lags in seasonal part
  int<lower=0> P_ar;
  // number of previous noise terms in seasonal part
  int<lower=0> Q_ma;
  // time frequency in seasonal part
  int<lower=0> ts_frequency;
  // boolean for optional state intercept phi_0
  // 1 is to include intercept
  int<lower = 0, upper = 1> include_state_intercept;
  
  // boolean for optional observation intercept d_0
  // 1 is to include intercept
  int<lower = 0, upper = 1> include_obs_intercept;
  
  // boolean for optional constraints for AR process
  // 1 is to add constraints
  int<lower = 0, upper = 1> stationary;
}

transformed data{
  int <lower=0> r = max(p_ar+P_ar*ts_frequency, q_ma+Q_ma*ts_frequency+1);
  
  // number of statess
  int<lower=0> m = r;
  
  // Z in 3.19
  matrix[p,m] Z = append_col([[1]],to_matrix(rep_row_vector(0, m-1)));
  
  // observation covariance
  matrix[p,p] H = [[0]];
  
}

parameters {
  
  // optional state intercept
  vector[include_state_intercept] phi_0;
  
  // optional observation intercept
  vector[include_obs_intercept] d_0;
  
  vector[p_ar] unconstrained_phi;
  
  vector[P_ar] unconstrained_phi_seasonal;

  vector[q_ma] unconstrained_theta;
  
  vector[Q_ma] unconstrained_theta_seasonal;
  
  real <lower = 0> var_zeta;
  
  // expected value for the initial state
  vector[m] a1;
  
}

transformed parameters{
  
  // state intercept
  vector[m] c;
  
  // observation intercept
  vector[p] d;
  
  matrix[1,1] Q = [[var_zeta]];
  
  matrix[m,1] R;
  
  matrix[m,m] T;
  
  matrix[m,m] P1;
  
  // initialize dummy_phi to 0s
  vector [r] dummy_phi = rep_vector(0, r);
  
  // initialize dummy_theta to 0s
  vector [r-1] dummy_theta = rep_vector(0, r-1);
  
  vector[p_ar] phi;
  
  vector[P_ar] phi_seasonal;

  vector[q_ma] theta = constrain_stationary(unconstrained_theta) ;
  
  vector[Q_ma] theta_seasonal = constrain_stationary(unconstrained_theta_seasonal);
  
  if (stationary) {
    phi = constrain_stationary(unconstrained_phi);
    phi_seasonal = constrain_stationary(unconstrained_phi_seasonal);
  }
  else {
    phi = unconstrained_phi;
    phi_seasonal = unconstrained_phi_seasonal;
  }
    
  
  if (p_ar > 0) {
    dummy_phi[1:p_ar] = phi;
  }
  
  if (P_ar > 0) {
    for (i in 1:P_ar){
      dummy_phi[i*ts_frequency] = phi_seasonal[i];
    }
  }
  
  if (p_ar > 0 && P_ar > 0){
    for (i in 1:p_ar){
      for (j in 1: P_ar){
        dummy_phi[i+ j*ts_frequency] = -phi[i]* phi_seasonal[j];
      }
    }
  }

  
  // fill in theta's
  if (q_ma > 0) {
    dummy_theta[1:q_ma] = theta;
  }  
  
  if (Q_ma > 0) {
    for (i in 1:Q_ma){
      dummy_theta[i*ts_frequency] = theta_seasonal[i];
    }
  }
  
  if (q_ma > 0 && Q_ma > 0) {
    for (i in 1:q_ma){
      for (j in 1: Q_ma){
        dummy_theta[i+ j*ts_frequency] = -theta[i]* theta_seasonal[j];
      }
    }
  }

  
  // state covariance R in 3.20
  R = append_row([[1]], to_matrix(dummy_theta));
  
  // T in 3.20
  T = append_col(dummy_phi,
  append_row(
    diag_matrix(rep_vector(1, m-1)), to_matrix(rep_row_vector(0, m-1))));
  
  
  // variance for the initial state
  P1 = var_zeta * stationary_cov(T, quad_form_sym(Q, R '));
  
  if (include_state_intercept){
    c = append_row(phi_0,rep_vector(0, m-1));
  }
  else {
    c = rep_vector(0, m);
  }
    
  if (include_obs_intercept) {
    d = append_row(d_0,rep_vector(0, p-1));
  }
  else {
    d = rep_vector(0, p);
  }
    
  
}

model {
  target += ssm_constant_lpdf (y| d, Z, H, c, T, R, Q, a1, P1);
  
}
