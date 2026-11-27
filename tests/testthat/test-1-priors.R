## The four priors (R/priors.R): three are one family, R-squared ~ Beta(a, 1)
## with a uniform direction; independent_partial_r gives each partial
## correlation its own flat prior

test_that("each prior has its shape, and joint_partial_r is the default", {
  expect_equal(prior_shape("joint_partial_r", 6), 3)
  expect_equal(prior_shape("flat_r", 6), 1 / 2)
  expect_equal(prior_shape("flat_r_squared", 6), 1)
  expect_error(prior_shape(NULL, 3), "Choose a prior")
  expect_error(check_prior("flat"), "Unknown prior")
  ## without a prior, the default is used and announced once per session
  default_prior_announced$done <- NULL
  set.seed(1)
  expect_message(by_default <- draw_coefficients_from_prior(10, diag(2)), "default prior")
  set.seed(1)
  expect_equal(by_default, draw_coefficients_from_prior(10, diag(2), prior = "joint_partial_r"))
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

test_that("independent_partial_r is a prior; flat_partial_r, its old name, is not", {
  expect_identical(check_prior("independent_partial_r"), "independent_partial_r")
  expect_true(prior_is_independent("independent_partial_r"))
  expect_false(any(vapply(c("joint_partial_r", "flat_r", "flat_r_squared"),
                          prior_is_independent, logical(1))))
  ## it has no Beta shape
  expect_error(prior_shape("independent_partial_r", 3), "no R-squared")
  expect_error(check_prior("flat_partial_r"), "Unknown prior")
  expect_error(draw_coefficients_from_prior(10, diag(2), prior = "flat_partial_r"),
               "Unknown prior")
})

test_that("every prior has a label, in the order of the prior names", {
  labels <- prior_labels()
  expect_identical(names(labels), prior_names)
  expect_true(all(nzchar(labels)))
  expect_false(anyDuplicated(unname(labels)) > 0)
})
