stanfunctions <- lapply("functions", function(model_name) {
  # create C++ code for stan model
  stan_file <- file.path("inst", "stan", "lib", paste0(model_name, ".stan"))
  model_code <- paste(c("functions {", readLines(stan_file), "}"), collapse = "\n")
#  model_code <-  paste(readLines(stan_file), collapse = "\n")
  stanc_ret <- rstan::stanc(model_code = model_code, model_name = "Stan Functions",
                     allow_undefined = TRUE, obfuscate_model_name = FALSE)
  rstan::expose_stan_functions(stanc_ret, rebuild = TRUE, verbose = FALSE, env = environment())
})
