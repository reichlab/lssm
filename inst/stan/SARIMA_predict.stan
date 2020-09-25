#include functions.stan

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
  vector[m] a1;
}

transformed data {
  // Z in 3.19 of Durbin and Koopman
  matrix[p, m] Z = append_col([[1]], to_matrix(rep_row_vector(0, m - 1)));
  
  // observation covariance
  matrix[p, p] H = [[0]];

  // state intercept
  vector[m] c;
  
  // observation intercept
  vector[p] d = rep_vector(0, p);
  
  // below are transformed parameters from model
  
  // covariance of state noise
  matrix[1, 1] Q = [[var_zeta]];
  
  // update matrices for state
  matrix[m, m] T;
  matrix[m, 1] R;
  
  // covariance of initial state
  matrix[m, m] P1;
  
  // initialize dummy_phi to 0s
  vector[r] dummy_phi = rep_vector(0, r);
  
  // initialize dummy_theta to 0s
  vector[r-1] dummy_theta = rep_vector(0, r-1);
  
  // shape of result
  int result_nrow;
  int result_ncol;
  int result_length;

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
    diag_matrix(rep_vector(1, m-1)), to_matrix(rep_row_vector(0, m-1))));
  
  // variance for the initial state
  P1 = var_zeta * stationary_cov(T, quad_form_sym(Q, R'));
  
  if (include_intercept) {
    c = append_row(phi_0, rep_vector(0, m - 1));
  } else{
    c = rep_vector(0, m);
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
