library(lssm)

context("initial transform")

################################################################################
## Tests for initial transformations
################################################################################

test_that("do_initial_transform/invert_initial_transform: box-cox", {
  ## data all positive
  y <- 1:10

  gamma <- 0.5
  bc_params <- car::powerTransform(y + gamma, family = "bcPower")
  y_trans <- do_initial_transform(
    y = y,
    transformation = "box-cox",
    transform_offset = gamma,
    bc_lambda = bc_params$lambda)

  y_trans2 <- car::bcPower(
    U = y + gamma,
   lambda = bc_params$lambda)

  y_out <- invert_initial_transform(
    y = y_trans,
    transformation = "box-cox",
    transform_offset = gamma,
    bc_lambda = bc_params$lambda)

  expect_equal(y_trans, y_trans2)
  expect_equal(y, y_out)



  ## data all non-negative
  y <- 0:10

  gamma <- 0.5
  bc_params <- car::powerTransform(y + gamma, family = "bcPower")
  y_trans <- do_initial_transform(
    y = y,
    transformation = "box-cox",
    transform_offset = gamma,
    bc_lambda = bc_params$lambda)

  y_trans2 <- car::bcPower(
    U = y + gamma,
    lambda = bc_params$lambda)

  y_out <- invert_initial_transform(
    y = y_trans,
    transformation = "box-cox",
    transform_offset = gamma,
    bc_lambda = bc_params$lambda)

  expect_equal(y_trans, y_trans2)
  expect_equal(y, y_out)
})

test_that("do_initial_transform/invert_initial_transform: log, no offset", {
  ## data all positive
  y <- 1:10

  y_trans <- do_initial_transform(
    y = y,
    transformation = "log",
    transform_offset = 0.0)

  y_trans2 <- log(y)

  y_out <- invert_initial_transform(
    y = y_trans,
    transformation = "log",
    transform_offset = 0.0)

  expect_equal(y_trans, y_trans2)
  expect_equal(y, y_out)
})

test_that("do_initial_transform/invert_initial_transform: log, offset", {
  ## data all positive
  y <- 1:10

  y_trans <- do_initial_transform(
    y = y,
    transformation = "log",
    transform_offset = 1.0)

  y_trans2 <- log(y + 1.0)

  y_out <- invert_initial_transform(
    y = y_trans,
    transformation = "log",
    transform_offset = 1.0)

  expect_equal(y_trans, y_trans2)
  expect_equal(y, y_out)
})

test_that("do_initial_transform/invert_initial_transform: transformation none",
  {
  ## data all positive
  y <- 1:10

  y_trans <- do_initial_transform(
    y = y,
    transformation = "none")

  y_trans2 <- y

  y_out <- invert_initial_transform(
    y = y_trans,
    transformation = "none")

  expect_equal(y_trans, y_trans2)
  expect_equal(y, y_out)
})
