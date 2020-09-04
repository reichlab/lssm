library(lssm)
library(ggplot2)

# simulated data
#phi_1 = 0.8
fake_data <- arima.sim(n = 200, model = list(ar = 0.8))
#ts.plot(fake_data, main = "Time Series Plot of Our Fake Data")


arma_fit <- fit_lssm(
  y = fake_data,
  model = "arma",
  transformation = "none",
  verbose = FALSE)

forecasts <- predict(
  arma_fit,
  horizon = 1,
  forecast_representation = "sample")
forecasts <- predict(
  arma_fit,
  horizon = 1,
  forecast_representation = "quantile")

plot_data <- data.frame(
  t = seq_len(200),
  y = as.vector(fake_data)
)

ggplot() +
  geom_line(data = plot_data, mapping = aes(x = t, y = y)) +
  geom_point(
    data = data.frame(t = 201, y = forecasts[1, 3]),
    mapping = aes(x = t, y = y),
    color = "purple"
  ) +
  geom_errorbar(
    data = data.frame(t = 201, ymin = forecasts[1, 2], ymax = forecasts[1, 4]),
    mapping = aes(x = t, ymin = ymin, ymax = ymax),
    color = "purple",
    width = 2
  ) +
  geom_errorbar(
    data = data.frame(t = 201, ymin = forecasts[1, 1], ymax = forecasts[1, 5]),
    mapping = aes(x = t, ymin = ymin, ymax = ymax),
    color = "purple",
    width = 3
  ) +
  theme_bw()



library(rstan)
rstan_options(auto_write = TRUE)
# AR1 fit
model <- stan_model("inst/stan_models/AR1.stan")

stan_data <- list(
  n = length(fake_data),
  p = 1,
  y = array(data=fake_data, dim = c(200,1))
)

estimate <- optimizing(model, stan_data)


# AR1 predict
model <- stan_model("inst/stan_models/AR1_predict.stan")

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

# ARp fit
fake_data_p <-arima.sim(n = 200, model = list(ar = c(.2, .5, .05)))
model <- stan_model("inst/stan_models/ARp.stan")

stan_data <- list(
  n = length(fake_data_p),
  p = 1,
  y = array(data=fake_data_p, dim = c(200,1)),
  p_ar =3
)

estimate_p <- optimizing(model, stan_data)

# ARp predict
model <- stan_model("inst/stan_models/ARp_predict.stan")

stan_data <- list(
  n = length(fake_data_p),
  p = 1,
  y = array(data=fake_data_p, dim = c(200,1)),
  steps_ahead = 1,
  p_ar = 3,
  # m = p_ar
  m = 3,
  phi = c(estimate_p$par["phi[1]"],estimate_p$par["phi[2]"],estimate_p$par["phi[3]"]),
  var_zeta =  estimate_p$par["var_zeta"],
  a1 = c(estimate_p$par["a1[1]"],estimate_p$par["a1[2]"],estimate_p$par["a1[3]"])
)
dim(stan_data$a1) = 3

prediction = optimizing(model, stan_data)

# ARMA fit 
fake_data_arma <- arima.sim(model = list(ar = c(.2, .5, .05), ma = .9 ), n = 200)
model <- stan_model("inst/stan_models/ARMA.stan")
stan_data <- list(
  n = length(fake_data_arma),
  p = 1,
  y = array(data=fake_data_arma, dim = c(200,1)),
  p_ar = 3,
  q_ma = 1
)

estimate_arma <- optimizing(model, stan_data)

# ARMA predict 
model <- stan_model("inst/stan_models/ARMA_predict.stan")

stan_data <- list(
  n = length(fake_data_arma),
  p = 1,
  y = array(data=fake_data_arma, dim = c(200,1)),
  steps_ahead = 1,
  p_ar = 3,
  q_ma = 1,
  # r = max(p_ar, q_ma+1)
  r = 3,
  # m = r
  m = 3,
  phi = c(estimate_arma$par["phi[1]"],estimate_arma$par["phi[2]"],estimate_arma$par["phi[3]"]),
  theta = estimate_arma$par["theta[1]"],
  var_zeta =  estimate_arma$par["var_zeta"],
  a1 = c(estimate_arma$par["a1[1]"],estimate_arma$par["a1[2]"],estimate_arma$par["a1[3]"] )
)
dim(stan_data$a1) = 3
dim(stan_data$phi) = 3
dim(stan_data$theta) = 1

prediction = optimizing(model, stan_data)

