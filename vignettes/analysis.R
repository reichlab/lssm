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
