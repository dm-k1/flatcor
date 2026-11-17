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
