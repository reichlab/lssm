functions{
#include /lib/functions.stan
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
  
  // number of statess
  int<lower=0> m = r;
  
  // Z in 3.19
  matrix[p, m] Z = append_col([[1]], to_matrix(rep_row_vector(0, m - 1)));
  
  // observation covariance
  matrix[p, p] H = [[0]];

  // expected value for the state at time 0
  vector[m] a0 = rep_vector(0, m);
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
  
  // expected value for the initial state
//  vector[m] a1;
}

transformed parameters{
  // state intercept
  vector[m] c;
  
  // observation intercept
  vector[p] d = rep_vector(0, p);
  
  // covariance of state noise
  matrix[1,1] Q = [[var_zeta]];
  
  // update matrices for state
  matrix[m,m] T;
  matrix[m,1] R;

  // expected value of state at time 1
  vector[m] a1;

  // covariance of state at time 1
  matrix[m,m] P1;
  
  // initialize dummy_phi to 0s
  vector [r] dummy_phi = rep_vector(0, r);
  
  // initialize dummy_theta to 0s
  vector [r-1] dummy_theta = rep_vector(0, r-1);
  
  // constrain parameters to stationary and nvertible
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
  
  // set AR coefficient values
  if (p_ar > 0) {
    dummy_phi[1:p_ar] = phi;
  }
  
  if (P_ar > 0) {
    for (i in 1:P_ar) {
      dummy_phi[i*ts_frequency] = phi_seasonal[i];
    }
  }
  
  if (p_ar > 0 && P_ar > 0) {
    for (i in 1:p_ar) {
      for (j in 1: P_ar) {
        dummy_phi[i+ j*ts_frequency] = -phi[i]* phi_seasonal[j];
      }
    }
  }
  

  // set MA coefficient values
  if (q_ma > 0) {
    dummy_theta[1:q_ma] = theta;
  }
  
  if (Q_ma > 0) {
    for (i in 1:Q_ma) {
      dummy_theta[i*ts_frequency] = theta_seasonal[i];
    }
  }
  
  if (q_ma > 0 && Q_ma > 0) {
    for (i in 1:q_ma) {
      for (j in 1: Q_ma) {
        dummy_theta[i+ j*ts_frequency] = -theta[i]* theta_seasonal[j];
      }
    }
  }
  

  // state covariance R in 3.20 of Durbin and Koopman
  R = append_row([[1]], to_matrix(dummy_theta));
  
  // T in 3.20 of Durbin and Koopman
  T = append_col(dummy_phi,
    append_row(
      diag_matrix(rep_vector(1, m - 1)), to_matrix(rep_row_vector(0, m - 1))));
  
  // variance for the state at time 0
  P1 = var_zeta * stationary_cov(T, quad_form_sym(Q, R'));
  
  // intercept for state transitions
  if (include_intercept) {
    c = append_row(phi_0, rep_vector(0, m - 1));
  } else {
    c = rep_vector(0, m);
  }

  // expected value and covariance of forecast distribution for state at time 1
  a1 = ssm_update_predicted_a(c, T, a0);
  P1 = ssm_update_predicted_P(P1, T, quad_form_sym(Q, R'));
}

model {
  if(horizon == 0) {
    target += ssm_constant_lpdf(y | d, Z, H, c, T, R, Q, a1, P1);
  } else {
    target += ssm_constant_forecast_lpdf(y | d, Z, H, c, T, R, Q, a1, P1, horizon);
  }
}
