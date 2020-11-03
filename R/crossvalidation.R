#' Define grid of tuning parameters for the arma model
#' 
#' @param y 
#' @param len 
#' @param search 
#' @param transform 
#' @param transform_offset 
#' @param max_d 
#' @param max_D 
#' @param include_intercept
#' @param max_p_ar 
#' @param max_q_ma 
#' @param max_P_ar
#' @param max_Q_ma
#' @param min_order
#' @param max_order 
#' @param stationary
#' 
#' @return data frame with columns model, transform, transform_offset, d, D,
#' p_ar, q_ma, P_ar, Q_ma
#' 
#' @export
sarima_param_grid <- function(
  y,
  len = NULL,
  search = "grid",
  transformation = "box-cox",
  transform_offset = if (any(x <= 0)) { -1 * min(x) + 0.49} else { 0.0 },
  max_d = 2,
  max_D = 1,
  include_intercept = c(FALSE, TRUE),
  max_p_ar = 5,
  max_q_ma = 5,
  max_P_ar = 2,
  max_Q_ma = 2,
  min_order = 1,
  max_order = 5,
  stationary = 1L
) {
  library(lssm)
  
  # Only do seasonal differencing if y has a seasonal period
  if(frequency(y) < 2) {
    max_D <- 0
  }
  
  # For grid search, all combinations of d, D, p_ar, and q_ma satisfying
  # p_ar + q_ma <= max_order
  out <- expand.grid(
    model = "arma",
    transformation = transformation,
    transform_offset = transform_offset,
    d = int_seq(from = 0, to = max_d, max_len = len),
    D = int_seq(from = 0, to = max_D, max_len = len),
    include_intercept = include_intercept,
    p_ar = int_seq(from = 0, to = max_p_ar, max_len = len),
    q_ma = int_seq(from = 0, to = max_q_ma, max_len = len),
    P_ar = int_seq(from = 0, to = max_P_ar, max_len = len),
    Q_ma = int_seq(from = 0, to = max_Q_ma, max_len = len),
    stationary = stationary,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(
      p_ar + q_ma + P_ar + Q_ma >= min_order,
      p_ar + q_ma + P_ar + Q_ma <= max_order
    )
  
  # If data contains non-positive values, remove grid combinations where the
  # transformation is log or box-cox and the transform_offset is less than or
  # equal to the smallest data value
  min_data <- min(y)
  if(min_data <= 0) {
    out <- out %>%
      dplyr::filter(
        !(transformation %in% c("log", "box-cox") &
            transform_offset <= min_data)
      )
  }
  
  if (search == "random") {
    # select a random set of rows from the grid defined above
    out <- out %>%
      dplyr::sample_n(size = len)
  }
  
  return(out)
}


#' Obtain cross-validated performance for lssm models
#' 
#' @export
crossvalidate_lssm <- function(
  y,
  ts_frequency,
  initial_window = 150,
  crossval_start_horizon = 1,
  crossval_end_horizon = 1,
  fixed_window = FALSE,
  crossval_frequency = 1,
  tune_grid,
  verbose = TRUE,
  parallel = FALSE
) {
  if (parallel && foreach::getDoParWorkers() > 1) {
    `%do_op%` <- `%dopar%`
  } else {
    `%do_op%` <- `%do%`
  }
  
  # Define training and validation sets for time series cross-validation
  crossval_folds <- caret::createTimeSlices(
    y,
    initialWindow = initial_window,
    horizon = crossval_end_horizon,
    fixedWindow = fixed_window)
  crossval_horizon_inds <- seq(
    from = crossval_start_horizon, to = crossval_end_horizon)
  if (crossval_frequency != 1) {
    subset_inds <- seq(
      from = 1,
      by = crossval_frequency,
      to = length(crossval_folds$train))
    crossval_folds$train <- crossval_folds$train[subset_inds]
    crossval_folds$test <- lapply(
      crossval_folds$test[subset_inds],
      function(inds) {
        inds[crossval_horizon_inds]
      })
  }
  
  # Cross-validate each model specification in parallel
  crossval_results <- foreach(
    model_ind = seq_len(nrow(tune_grid)),
    .combine = rbind,
    .verbose = verbose) %dopar% {
    model_results <- tune_grid[model_ind, ] %>%
      dplyr::mutate(join_col = "a") %>%
      dplyr::right_join(
        data.frame(
          fold = seq_along(crossval_folds$train),
          log_score = NA_real_,
          run_time = NA_real_,
          join_col = "a",
          stringsAsFactors = FALSE
        ),
        by = "join_col"
      ) %>%
      dplyr::select(-join_col)
    
    # For first fit, initial parameter values are all 0;
    # later fits will use parameter estimates from the previous fit.
    init_par <- NULL
    
    # Iterate through folds and obtain log score for each.
    for (fold_ind in seq_along(crossval_folds$train)) {
      if (verbose) {
        print(paste0("Estimating fold ", fold_ind, " of ", length(crossval_folds$train),
                     " for model ", model_ind, " of ", nrow(tune_grid)))
      }
      
      tic <- Sys.time()
      # fit model to training set
      param <- as.list(tune_grid[model_ind, ])
      param$y <- y[crossval_folds$train[[fold_ind]]]
      param$ts_frequency <- ts_frequency
      param$verbose <- verbose
      if(!is.null(init_par)) {
        param$init_par <- init_par
      }
      
      model_fit <- do.call(fit_lssm, param)
      
      # save parameter estimates to use as initial values for estimation with
      # next fold
      init_par <- model_fit$par
      
      # validation set predictions
      named_dist_forecast <- predict(
        model_fit,
        horizon = crossval_end_horizon,
        forecast_representation = "named_dist")
      
      # log score for validation set predictions
      dfun <- paste0("d", named_dist_forecast$family[[1]])
      if(dfun == "dmvnorm") {
        dfun <- mvtnorm::dmvnorm
      }
      arg_names <- names(named_dist_forecast)[
        !(names(named_dist_forecast) %in% c("family", "h"))]
      call_args <- purrr::map(
        arg_names,
        function(pan) {
          named_dist_forecast[[pan]][[1]]
        }
      )
      names(call_args) <- arg_names
      call_args$mean <- call_args$mean[crossval_horizon_inds]
      call_args$sigma <- call_args$sigma[
        crossval_horizon_inds, crossval_horizon_inds, drop = FALSE]
      call_args$x <- y[crossval_folds$test[[fold_ind]]]
      call_args$log <- TRUE
      
      model_results$log_score[fold_ind] <- mean(do.call(dfun, call_args))
      
      toc <- Sys.time()
      model_results$run_time[fold_ind] <- toc - tic
    }
    
    return(model_results)
  }
}
