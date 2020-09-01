# simulated data
#phi_1 = 0.8
fake_data <-arima.sim(n = 200, model = list(ar = 0.8))
ts.plot(fake_data, main = "Time Series Plot of Our Fake Data")

library(rstan)
rstan_options(auto_write = TRUE)
# fit
model <- stan_model("inst/stan_models/AR1.stan")

stan_data <- list(
  n = length(fake_data),
  p = 1,
  y = array(data=fake_data, dim = c(200,1))
)

estimate <- optimizing(model, stan_data)


# predict
model <- stan_model("inst/stan_models/predict.stan")

stan_data <- list(
  n = length(fake_data),
  p = 1,
  y = array(data=fake_data, dim = c(200,1)),
  steps_ahead = 1,
  m = 1,
  phi_1 = estimate$par["phi_1"],
  var_zeta =  estimate$par["var_zeta"],
  a1 = estimate$par["a1[1]"]
)
dim(stan_data$a1) = 1

prediction = optimizing(model, stan_data)
