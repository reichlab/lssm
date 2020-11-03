  /**
  ---
  function: sarima_build_observation_matrices
  args:
  - name: m
    description: Length of state vector
  returns: A $p \times 1 + m + p$ matrix containing d, Z, and H
  ---
  
  */
  matrix sarima_build_observation_matrices(int m) {
    int p = 1;
    
    // intercept
    vector[p] d = rep_vector(0, p);

    // Z in 3.19
    matrix[p, m] Z = append_col([[1]], to_matrix(rep_row_vector(0, m - 1)));
    
    // observation covariance
    matrix[p, p] H = [[0]];
    
    // return value
    matrix[p, 1 + m + p] result = append_col(append_col(to_matrix(d), Z), H);
    
    return result;
  }
  
  
  /**
  ---
  function: sarima_build_state_matrices
  args:
  - name: m
    description: Length of state vector
  returns: A $p \times 1 + m + p$ matrix containing d, Z, and H
  ---
  
  */
  matrix sarima_build_state_matrices(
    int p_ar, int q_ma, int P_ar, int Q_ma, int ts_frequency,
    int include_intercept,
    vector phi_0, vector phi, vector phi_seasonal,
    vector theta, vector theta_seasonal,
    real var_zeta) {
    // number of states
    int r =
      max(p_ar + P_ar * ts_frequency, q_ma + Q_ma * ts_frequency + 1);
    int m = r;
    
    // intercept
    vector[m] c;
    
    // covariance of state noise
    matrix[1,1] Q = [[var_zeta]];
    
    // update matrices for state
    matrix[m,m] T;
    matrix[m,1] R;
  
    // expected value for the state at time 0
    vector[m] a0 = rep_vector(0, m);
    
    // expected value of state at time 1
    vector[m] a1;
  
    // covariance of state at time 1
    matrix[m,m] P1;
    
    // return value
    matrix[m, 1 + m + 1 + 1 + 1 + m] result;
    
    // initialize dummy_phi to 0s
    vector [r] dummy_phi = rep_vector(0, r);
    
    // initialize dummy_theta to 0s
    vector [r-1] dummy_theta = rep_vector(0, r-1);
    
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

    // return value
    result[, 1] = c;
    result[, (1 + 1):(1 + m)] = T;
    result[, (1 + m + 1):(1 + m + 1)] = R;
    result[1, 1 + m + 1 + 1] = Q[1, 1];
    result[, (1 + m + 1 + 1 + 1)] = a1;
    result[, (1 + m + 1 + 1 + 1 + 1):(1 + m + 1 + 1 + 1 + m)] = P1;

    return result;
  }
