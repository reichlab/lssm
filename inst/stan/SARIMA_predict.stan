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
  int <lower=0> ts_frequency;
  
  // boolean for optional state intercept phi_0
  // 1 is to include intercept
  int<lower = 0, upper = 1> include_intercept;
  
  // number of steps forward to predict
  int <lower=1> horizon;
  
  // What type of prediction to compute:
  // marginal distribution at each forecast horizon (joint = 0),
  // or joint across all horizons (joint = 1)
  int<lower=0, upper=1> joint;

  // optional state intercept
  vector[include_intercept] phi_0;
  
  // parameters
  vector[p_ar] phi;
  
  vector[P_ar] phi_seasonal;

  vector[q_ma] theta;
  
  vector[Q_ma] theta_seasonal;
  
  real <lower = 0> var_zeta;
  
  // r = max(p_ar+P_ar*ts_frequency, q_ma_+Q_ma*ts_frequency+1)
  int <lower=0> r;
  
  // number of states
  // m = r
  int<lower=0> m;
  
  // expected value for the initial state
//  vector[m] a1;
}

transformed data {
  // observation matrices
  matrix[p, 1 + m + p] observation_matrices =
    sarima_build_observation_matrices(m);
  
  // observation intercept
  vector[p] d = observation_matrices[, 1];
  
  // state to observation mean
  matrix[p, m] Z = observation_matrices[, (1 + 1):(1 + m)];
  
  // observation covariance
  matrix[p, p] H = observation_matrices[, (1 + m + 1):(1 + m + p)];

//  // model parameters after constraining to stationary and invertible
//  vector[p_ar] phi;
//  vector[P_ar] phi_seasonal;
  
//  vector[q_ma] theta = constrain_stationary(unconstrained_theta) ;
//  vector[Q_ma] theta_seasonal = constrain_stationary(unconstrained_theta_seasonal);
  
  // state process matrices
  vector[m] c;
  matrix[m,m] T;
  matrix[m,1] R;
  matrix[1,1] Q;
  vector[m] a1;
  matrix[m,m] P1;
  
  // shape of result
  int result_nrow;
  int result_ncol;
  int result_length;

//  if (stationary) {
//    phi = constrain_stationary(unconstrained_phi);
//    phi_seasonal = constrain_stationary(unconstrained_phi_seasonal);
//  } else {
//    phi = unconstrained_phi;
//    phi_seasonal = unconstrained_phi_seasonal;
//  }

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

  // shape of result
  if(joint == 0) {
    // marginal distribution parameters separately per horizon
    result_nrow = p;
    result_ncol = p + 1;
    result_length = horizon;
  } else {
    // joint distribution for all horizons
    result_nrow = horizon * p;
    result_ncol = horizon * p + 1;
    result_length = 1;
  }
}

parameters {
  
}

generated quantities {
  matrix[result_nrow, result_ncol] forecasts[result_length];
  if(joint == 0) {
    forecasts = predict(y, d, Z, H, c, T, R, Q, a1, P1, horizon);
  } else {
    forecasts[1] = ssm_constant_joint_predict(y, d, Z, H, c, T, R, Q, a1, P1, horizon);
  }
}
