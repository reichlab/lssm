functions {
#include /lib/functions.stan
#include /lib/sarima_utils.stan
}

data {
  // number of observations
  int<lower=0> n;
  // length for each observation
  int<lower=0> p;
  // n obs each is px1
  vector[p] y[n];
  // number of lags in non-seasonal auto-regressive part
  int<lower=0> p_ar;
  // number of terms in non-seasonal moving average part
  int<lower=0> q_ma;
  // number of lags in seasonal auto-regressive part
  int<lower=0> P_ar;
  // number of terms in seasonal moving average part
  int<lower=0> Q_ma;
  
  // time frequency in seasonal part
  int<lower=0> ts_frequency;
  
  // boolean for optional state intercept phi_0
  // 1 is to include intercept
  int<lower = 0, upper = 1> include_intercept;

  // boolean for optional constraints for AR process
  // 1 is to add constraints
  int<lower = 0, upper = 1> stationary;

  // horizon targeted for estimation
  // 0 means maximum likelihood estimation
  // horizon > 0 calculates the log score of the predictive distribution at the
  // specified horizon
  int<lower = 0> horizon;
}

transformed data {
  int <lower=0> r =
    max(p_ar + P_ar * ts_frequency, q_ma + Q_ma * ts_frequency + 1);
  
  // number of states
  int<lower=0> m = r;
  
  // observation matrices
  matrix[p, 1 + m + p] observation_matrices =
    sarima_build_observation_matrices(m);
  
  // observation intercept
  vector[p] d = observation_matrices[, 1];
  
  // state to observation mean
  matrix[p, m] Z = observation_matrices[, (1 + 1):(1 + m)];
  
  // observation covariance
  matrix[p, p] H = observation_matrices[, (1 + m + 1):(1 + m + p)];
}

parameters {
  // optional state intercept
  vector[include_intercept] phi_0;
  
  // non-seasonal AR coefficients
  vector[p_ar] unconstrained_phi;
  
  // seasonal AR coefficients
  vector[P_ar] unconstrained_phi_seasonal;

  // non-seasonal MA coefficients
  vector[q_ma] unconstrained_theta;
  
  // seasonal MA coefficients
  vector[Q_ma] unconstrained_theta_seasonal;
  
  // variance for MA noise terms
  real <lower = 0> var_zeta;
}

transformed parameters{
  // model parameters after constraining to stationary and invertible
  vector[p_ar] phi;
  vector[P_ar] phi_seasonal;
  
//  vector[q_ma] theta = constrain_stationary(unconstrained_theta) ;
//  vector[Q_ma] theta_seasonal = constrain_stationary(unconstrained_theta_seasonal);
  vector[q_ma] theta = unconstrained_theta;
  vector[Q_ma] theta_seasonal = unconstrained_theta_seasonal;
  
  // state process matrices
  vector[m] c;
  matrix[m,m] T;
  matrix[m,1] R;
  matrix[1,1] Q;
  vector[m] a1;
  matrix[m,m] P1;
  
  if (stationary) {
    phi = constrain_stationary(unconstrained_phi);
    phi_seasonal = constrain_stationary(unconstrained_phi_seasonal);
  } else {
    phi = unconstrained_phi;
    phi_seasonal = unconstrained_phi_seasonal;
  }
//  print("modified function");
//  print("Stationary = ");
//  print(stationary);
//  print("phi = ");
//  print(phi);
//  print("phi_seasonal = ");
//  print(phi_seasonal);
//  print("theta = ");
//  print(theta);
//  print("theta_seasonal = ");
//  print(theta_seasonal);

  // state matrices
  {
    matrix[m, 1 + m + 1 + 1 + 1 + m] state_matrices = sarima_build_state_matrices(
      p_ar, q_ma, P_ar, Q_ma, ts_frequency, include_intercept,
      phi_0, phi, phi_seasonal, theta, theta_seasonal, var_zeta);
    
    // state intercept
    c = state_matrices[, 1];
  
    // update matrices for state
    T = state_matrices[, (1 + 1):(1 + m)];
    R = state_matrices[, (1 + m + 1):(1 + m + 1)];
    
    // covariance of state noise
    Q = state_matrices[1:1, (1 + m + 1 + 1):(1 + m + 1 + 1)];
  
    // expected value of state at time 1
    a1 = state_matrices[, (1 + m + 1 + 1 + 1)];
  
    // covariance of state at time 1
    P1 = state_matrices[, (1 + m + 1 + 1 + 1 + 1):(1 + m + 1 + 1 + 1 + m)];
  }
}

model {
  if(horizon == 0) {
    real temp;
//    print("calculating log likelihood...");
    temp = ssm_constant_lpdf(y | d, Z, H, c, T, R, Q, a1, P1);
//    print("log likelihood = ");
//    print(temp);
    target += ssm_constant_lpdf(y | d, Z, H, c, T, R, Q, a1, P1);
  } else {
    target += ssm_constant_forecast_lpdf(y | d, Z, H, c, T, R, Q, a1, P1, horizon);
  }
}
