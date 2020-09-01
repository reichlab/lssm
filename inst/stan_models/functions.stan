functions{
  // Kronecker product
  // @param matrix A: An m×n matrix
  // @param matrix B: A p×q matrix
  //
  // @return matrix: An mp×nq matrix
  matrix kronecker_prod(matrix A, matrix B) {
    matrix[rows(A) * rows(B), cols(A) * cols(B)] C;
    int m;
    int n;
    int p;
    int q;
    m = rows(A);
    n = cols(A);
    p = rows(B);
    q = cols(B);
    for (i in 1:m) {
      for (j in 1:n) {
        int row_start;
        int row_end;
        int col_start;
        int col_end;
        row_start = (i - 1) * p + 1;
        row_end = (i - 1) * p + p;
        col_start = (j - 1) * q + 1;
        col_end = (j - 1) * q + 1;
        C[row_start:row_end, col_start:col_end] = A[i, j] * B;
      }
    }
    return C;
  }
  
  // Convert vector to a matrix (column-major)
  // @param vector v: An n×m vector
  // @param int m: Number of rows in the vector
  // @param int n: Number of columns in the vector
  //
  // @return matrix: A m×n matrix containting the elements from v
  matrix to_matrix_colwise(vector v, int m, int n) {
    matrix[m, n] res;
    for (j in 1:n) {
     for (i in 1:m) {
       res[i, j] = v[(j - 1) * m + m];
     }
    }
    return res;
  }
  
  // Find the covariance of the stationary distribution of an ARMA model
  // @param matrix T: The m×m transition matrix
  // @param matrix R: The m×q system disturbance selection matrix
  //
  // @return matrix: An  m×m matrix with the stationary covariance matrix
  matrix arima_stationary_cov(matrix T, matrix R) {
    matrix[rows(T), cols(T)] Q0;
    matrix[rows(T) * rows(T), rows(T) * rows(T)] TT;
    vector[rows(T) * rows(T)] RR;
    int m;
    int m2;
    m = rows(T);
    m2 = m * m;
    RR = to_vector(tcrossprod(R));
    TT = kronecker_prod(T, T);
    Q0 = to_matrix_colwise((diag_matrix(rep_vector(1.0, m2)) - TT) \ RR, m, m);
    return Q0;
  }
  
  // Update the expected value of the predicted state
  // @param vector a: An m×1 vector with the predicted state a_t
  // @param vector c: An m×1 vector with the system intercept c_t
  // @param matrix T: An m×m matrix with the transition matrix T_t
  // @param vector v: A p×1 vector with the forecast error v_t
  // @param matrix K: An m×p matrix with the Kalman gain K_t
  //
  // @return vector: A m×1 vector with the predicted state at t+1, a_t+1
  vector ssm_filter_update_a(vector a, vector c, matrix T, vector v, matrix K) {
    vector[num_elements(a)] a_new;
    a_new = T * a + K * v + c;
    return a_new;
  }
  
  // Ensure a matrix is symmetrix
  // @param x: An n×n matrix
  //
  // @return An n×n symmetric matrix:  0.5(x+x′)
  matrix to_symmetric_matrix(matrix x) {
    return 0.5 * (x + x ');
  }
  
  // Update the variance of the predicted state
  // @param matrix P: An m×m vector with the variance of the prected state P_t
  // @param matrix Z: A p×m matrix with the design matrix Z_t
  // @param matrix T: An m×m matrix with the transition matrix T_t
  // @param matrix RQR: A m×m matrix with the system covariance matrix R_tQ_tR′_t
  // @param matrix K: An m×p matrix with the Kalman gain K_t
  //
  // @return matrix: An m×1 vector with the predicted state at t+1, a_t+1
  matrix ssm_filter_update_P(matrix P, matrix Z, matrix T,
                           matrix RQR, matrix K) {
    matrix[rows(P), cols(P)] P_new;
    P_new = to_symmetric_matrix(T * P * (T - K * Z)' + RQR);
    return P_new;
  }
  
  // Calculate the number of unique elements in a symmetric matrix
  // @param 
  //
  // @return int: The number of unique elements
  int symmat_size(int n) {
    int sz;
    // This calculates it iteratively because Stan gives a warning
    // with integer division.
    sz = 0;
    for (i in 1:n) {
      sz = sz + i;
    }
    return sz;
  }
  
  // Indexes of the return values of the Kalman filter functions: ssm_filter
  // @param int m: The number of states
  // @param int p: The size of the observation vector y_t
  //
  // @return int[,] : A 6×3 integer array containing
  // the indexes of the return values of the Kalman filter.
  int[,] ssm_filter_idx(int m, int p) {
    int sz[6, 3];
    // loglike
    sz[1, 1] = 1;
    // v
    sz[2, 1] = p;
    // Finv
    sz[3, 1] = symmat_size(p);
    // K
    sz[4, 1] = m * p;
    // a
    sz[5, 1] = m;
    // P
    sz[6, 1] = symmat_size(m);
    // Fill in start and stop points
    sz[1, 2] = 1;
    sz[1, 3] = sz[1, 2] + sz[1, 1] - 1;
    for (i in 2:6) {
      sz[i, 2] = sz[i - 1, 3] + 1;
      sz[i, 3] = sz[i, 2] + sz[i, 1] - 1;
    }
    return sz;
  }
  
  // Update the precision of the forcast error 
  // @param matrix P: An m×m vector with the variance of the prected state P_t
  // @param matrix Z: A p×m matrix with the design matrix Z_t
  // @param matrix H: A p×p matrix with the observation covariance matrix H_t
  //
  // @return matrix: A p×p vector with F_t
  matrix ssm_filter_update_F(matrix P, matrix Z, matrix H) {
    matrix[rows(H), cols(H)] F;
    F = quad_form(P, Z') + H;
    return F;
  }
  
  // Update the precision of the forcast error (inversed)
  // @param matrix P: An m×m vector with the variance of the prected state P_t
  // @param matrix Z: A p×m matrix with the design matrix Z_t
  // @param matrix H: A p×p matrix with the observation covariance matrix H_t
  //
  // @return matrix: A p×p vector with  inversed F_t
  matrix ssm_filter_update_Finv(matrix P, matrix Z, matrix H) {
    matrix[rows(H), cols(H)] Finv;
    Finv = inverse(ssm_filter_update_F(P, Z, H));
    return Finv;
  }
  
  // Update the Kalman gain K_t
  // @param matrix P: An m×m vector with the variance of the prected state P_t
  // @param matrix Z: A p×m matrix with the design matrix Z_t
  // @param matrix T: An m×m matrix with the transition matrix T_t
  // @param matrix Finv: A p×p matrix
  //
  // @return matrix: An m×p matrix with the Kalman gain K_t
  matrix ssm_filter_update_K(matrix P, matrix Z, matrix T, matrix Finv) {
    matrix[cols(Z), rows(Z)] K;
    K = T * P * Z' * Finv;
    return K;
  }
 
  // Calculate the log-likelihood for a period
  // @param vector v: A p×1 matrix with the forecast error v_t
  // @param matrix Finv: A p×p matrix with variance of the forecast 
  // error inv(F_t)
  //
  // @return real: An m×m matrix L_t
  real ssm_filter_update_ll(vector v, matrix Finv) {
    real ll;
    int p;
    p = num_elements(v);
    // det(A^{-1}) = 1 / det(A) -> log det(A^{-1}) = - log det(A)
    ll = (- 0.5 *
          (p * log(2 * pi())
          - log_determinant(Finv)
          + quad_form(Finv, v)
        ));
    return ll;
  }
  
  // Check if two matrices are approximately equal
  // @param matrix A: An m×n matrix
  // @param matrix B: An m×n matrix
  // @param real: The relative tolerance for convergence
  //
  // @return int: If converged, then 1, else 0.
  int ssm_check_matrix_equal(matrix A, matrix B, real tol) {
    real eps;
    eps = max(to_vector(A - B)) / max(to_vector(A));
    if (eps < tol) {
      return 1;
    } else {
      return 0;
    }
  }
  
  // Update the forcast error
  // 
  vector ssm_filter_update_v(vector y, vector a, vector d, matrix Z) {
    vector[num_elements(y)] v;
    v = y - Z * a - d;
    return v;
  }
  
  // Log-likelihood of a Time-Invariant Linear Gaussian State Space Model
  // @param vector[]: y Observations y_t. An array of size n of p×1 vectors
  // @param vector d: Observation intercept d_t. An array of p×1 vectors
  // @param matrix Z: Design matrix Z_t. An array of p×m matrices.
  // @param matrix H: Observation covariance matrix H_t. An array of p×p matrices
  // @param vector c: State intercept c_t.An array of m×1 vectors
  // @param matrix T: Transition matrix T_t. An array of m×m matrices
  // @param matrix R: State covariance selection matrix R_t. An array of p×q matrices
  // @param matrix Q: State covariance matrix Q_t. An array of q×q matrices
  // @param vector a1: Expected value of the intial state, a1=E(α1). An m×1 matrix
  // @param matrix P1: Variance of the initial state, P1=Var(α1). An m×m matrix
  //
  // @return real The log-likelihood, p(y1:n|d,Z,H,c,T,R,Q), marginalized over the latent states
  real ssm_constant_lpdf(vector[] y,
                      vector d, matrix Z, matrix H,
                      vector c, matrix T, matrix R, matrix Q,
                      vector a1, matrix P1) {
    real ll;
    int n;
    int m;
    int p;

    n = size(y); // number of obs
    m = cols(Z);
    p = rows(Z);
    {
      vector[n] ll_obs;
      vector[m] a;
      matrix[m, m] P;
      vector[p] v;
      matrix[p, p] Finv;
      matrix[m, p] K;
      matrix[m, m] RQR;
      // indicator for if the filter has converged
      // This only works for time-invariant state space models
      int converged;
      matrix[m, m] P_old;
      real tol;
      converged = 0;
      tol = 1e-7;

      RQR = quad_form(Q, R);
      a = a1;
      P = P1;
      for (t in 1:n) {
        v = ssm_filter_update_v(y[t], a, d, Z);
        if (converged < 1) {
          Finv = ssm_filter_update_Finv(P, Z, H);
          K = ssm_filter_update_K(P, Z, T, Finv);
        }
        ll_obs[t] = ssm_filter_update_ll(v, Finv);
        // don't save a, P for last iteration
        if (t < n) {
          a = ssm_filter_update_a(a, c, T, v, K);
          // check for convergence
          // should only check for convergence if there are no missing values
          if (converged < 1) {
            P_old = P;
            P = ssm_filter_update_P(P, Z, T, RQR, K);
            converged = ssm_check_matrix_equal(P, P_old, tol);
          }
        }
      }
      ll = sum(ll_obs);
    }
    return ll;
  }

  matrix predict (vector[] y, vector d, matrix Z, matrix H,
                      vector c, matrix T, matrix R, matrix Q,
                      vector a1, matrix P1,
                      int steps_ahead){
    int m = cols(Z);
    int n = size(y);
    int p = rows(Z);
    
    vector[m] a;
    matrix[p, p] F;
    matrix[m, m] P;
    {
      vector[p] v;
      matrix[p, p] Finv;
      matrix[m, p] K;
      matrix[m, m] RQR;

      RQR = quad_form(Q, R);
      a = a1;
      P = P1;
      for (t in 1:(n+steps_ahead-1)){
        v = ssm_filter_update_v(y[t], a, d, Z);
        Finv = ssm_filter_update_Finv(P, Z, H);
        K = ssm_filter_update_K(P, Z, T, Finv);
        a = ssm_filter_update_a(a, c, T, v, K);
        P = ssm_filter_update_P(P, Z, T, RQR, K);
      }
      F = ssm_filter_update_F(P, Z, H);
    }
    
    // return p x (p+1)  first col: Z*a 
    return append_col(to_matrix(Z*a), F);
  }
  
}
