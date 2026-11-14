## The three priors as one family, R-squared ~ Beta(a, 1) with a uniform
## direction (R/priors.R): fast versions of checks in objective-bayes,
## studies 2 and 4

test_that("each prior has its shape, and flat_partial_r is the default", {
  expect_equal(prior_shape("flat_partial_r", 6), 3)
  expect_equal(prior_shape("flat_r", 6), 1 / 2)
  expect_equal(prior_shape("flat_r_squared", 6), 1)
  expect_error(prior_shape(NULL, 3), "Choose a prior")
  expect_error(log_bayes_factor_vs_intercept_only(0.2, 50, 3, prior = "flat"), "Unknown prior")
  ## without a prior, the default is used and announced once per session
  default_prior_announced$done <- NULL
  expect_message(by_default <- log_bayes_factor_vs_intercept_only(0.2, 50, 3), "default prior")
  expect_equal(by_default, log_bayes_factor_vs_intercept_only(0.2, 50, 3, prior = "flat_partial_r"))
})

test_that("rows dropped for missing values are reported", {
  with_missing <- mtcars
  with_missing$hp[c(2, 5)] <- NA
  expect_warning(bayes_factor_models(lm(mpg ~ wt + hp, with_missing),
                                     lm(mpg ~ wt, with_missing[-c(2, 5), ]),
                                     prior = "flat_partial_r"), "2 rows with a missing value")
})

test_that("flat R-squared: the whole-model Bayes factor agrees with adaptive integration", {
  ## the block likelihood ratio integrated against the prior density of |f|,
  ## 2 f / (1 + f^2)^2 (R^2 = f^2 / (1 + f^2) flat), by integrate()
  by_integration <- function(r_squared, n, k) {
    integrand <- function(f) exp(log_likelihood_ratio_block(f, r_squared, n, k, k)) *
      2 * f / (1 + f^2)^2
    peak <- sqrt(r_squared / (1 - r_squared)); width <- 1 / sqrt(n)
    breaks <- sort(unique(c(0, max(0, peak - 20 * width), peak, peak + 20 * width, Inf)))
    log(sum(sapply(seq_len(length(breaks) - 1), function(i)
      integrate(integrand, breaks[i], breaks[i + 1], rel.tol = 1e-10, subdivisions = 2000)$value)))
  }
  for (case in list(c(0.1, 50, 1), c(0.3, 50, 3), c(0.001, 500, 3), c(0.05, 200, 5)))
    expect_equal(log_bayes_factor_vs_intercept_only(case[1], case[2], case[3], prior = "flat_r_squared"),
                 by_integration(case[1], case[2], case[3]), tolerance = 1e-4)
})

test_that("the one route holds at very strong evidence (large-argument Bessel function)", {
  for (case in list(c(0.99, 10000, 3), c(0.99, 10000, 30), c(0.999, 10000, 10)))
    expect_equal(log_bayes_factor_vs_intercept_only(case[1], case[2], case[3], prior = "flat_r"),
                 log_bayes_factor_flat_r_by_directions(case[1], case[2], case[3]), tolerance = 1e-6)
  expect_true(is.finite(log_bayes_factor_vs_intercept_only(0.99, 10000, 30, prior = "flat_r_squared")))
})
