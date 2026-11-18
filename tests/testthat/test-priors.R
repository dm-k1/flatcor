## The three priors as one family, R-squared ~ Beta(a, 1) with a uniform
## direction (R/priors.R)

test_that("each prior has its shape, and flat_partial_r is the default", {
  expect_equal(prior_shape("flat_partial_r", 6), 3)
  expect_equal(prior_shape("flat_r", 6), 1 / 2)
  expect_equal(prior_shape("flat_r_squared", 6), 1)
  expect_error(prior_shape(NULL, 3), "Choose a prior")
  expect_error(check_prior("flat"), "Unknown prior")
  ## without a prior, the default is used and announced once per session
  default_prior_announced$done <- NULL
  set.seed(1)
  expect_message(by_default <- draw_coefficients_from_prior(10, diag(2)), "default prior")
  set.seed(1)
  expect_equal(by_default, draw_coefficients_from_prior(10, diag(2), prior = "flat_partial_r"))
})

test_that("the tests' sample R-squared is that of the least-squares fit", {
  set.seed(1)
  predictors <- matrix(rnorm(60), 30)
  y <- drop(predictors %*% c(0.5, -0.3)) + rnorm(30)
  ## a perfect fit gives 1, a fit to unrelated noise gives a small value
  expect_equal(sample_r_squared(drop(predictors %*% c(1, 2)) + 3, predictors), 1)
  expect_lt(sample_r_squared(rnorm(30), predictors), 0.3)
  ## equals the squared correlation of y with its fitted values
  fitted <- y - qr.resid(qr(cbind(1, predictors)), y)
  expect_equal(sample_r_squared(y, predictors), cor(y, fitted)^2)
})
