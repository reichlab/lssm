// parameter estimation/fitting
data {
  // number of observations
  int<lower=0> n;
  
  // length for each observation
  int<lower=0> p;
  
  // n obs each is px1
  vector[p] y[n];
  
  // separate mean vector of length p for each observation
  vector[p] mu[n];
  
  // Cholesky factor for covariance of each observation
  matrix[p, p] L[n];
}

parameters {
  
}

generated quantities{
  real results[n];
  for(i in 1:n) {
    results[i] = multi_normal_cholesky_lpdf(y[i] | mu[i], L[i]);
  }
}
