## bayes_factor_models(): checks of the inputs (no Stan needed), and the
## Bayes factor by bridge sampling against independent references
## (helper-references.R) and the intercept-only model's closed form (Stan).

example_data <- local({
  set.seed(42)
  correlations <- matrix(0.4, 3, 3)
  diag(correlations) <- 1
  predictors <- matrix(rnorm(40 * 3), 40) %*% chol(correlations)
  data <- data.frame(x1 = predictors[, 1], x2 = predictors[, 2], x3 = predictors[, 3])
  data$y <- 0.4 * data$x1 + 0.2 * data$x2 + rnorm(40)
  data
})

## Bayes factors with cmdstanr take the log density from BridgeStan
skip_without_stan <- function() {
  skip_on_cran()
  skip_if_not(requireNamespace("cmdstanr", quietly = TRUE), "needs cmdstanr")
  skip_if_not(requireNamespace("bridgestan", quietly = TRUE), "needs bridgestan")
}

## a Bayes factor estimate agrees with a reference within its Monte Carlo
## error: 5 times bridgesampling's approximate error (which can be optimistic),
## and at least 0.02 on the log scale
expect_within_error <- function(estimate, reference, error) {
  expect_lt(abs(estimate - reference), max(5 * error, 0.02))
}

test_that("the inputs are checked before anything is fitted", {
  with_missing <- example_data
  with_missing$x2[c(2, 5)] <- NA
  expect_warning(rows <- complete_rows(with_missing, list(y ~ x1, y ~ x2)),
                 "2 of 40 rows have a missing value")
  expect_identical(nrow(rows), 38L)
  expect_identical(nrow(complete_rows(with_missing, list(y ~ x1))), 40L)
  expect_error(complete_rows(example_data, list(y ~ x9)), "Not in the data: x9")
  expect_error(predictors_and_outcome(y ~ 0 + x1, example_data), "intercept")
  expect_identical(ncol(predictors_and_outcome(y ~ 1, example_data)$predictors), 0L)
  expect_error(bayes_factor_models(y ~ x1, x2 ~ x1, example_data, prior = "flat_r"),
               "same outcome")
  expect_error(bayes_factor_models(y ~ x1, y ~ 1, example_data, prior = "flat"),
               "Unknown prior")
  expect_error(check_full_rank(cbind(example_data$x1, example_data$x2,
                                     example_data$x1 + example_data$x2)), "linearly dependent")
  expect_true(check_full_rank(matrix(0, 40, 0)))
})

test_that("the intercept-only model's marginal likelihood matches its closed form", {
  skip_without_stan()
  ## flat intercept, 1/sigma noise: n^(-1/2) (2 pi)^(-(n-1)/2) Gamma((n-1)/2) / 2 (SST/2)^(-(n-1)/2)
  y <- example_data$y
  n <- length(y)
  sst <- sum((y - mean(y))^2)
  closed_form <- -0.5 * log(n) - (n - 1) / 2 * log(2 * pi) + lgamma((n - 1) / 2) - log(2) -
    (n - 1) / 2 * log(sst / 2)
  sampled <- sample_posterior(matrix(0, n, 0), y, "flat_r", chains = 4, iter_warmup = 1000,
                              iter_sampling = 5000, adapt_delta = 0.95, seed = 1,
                              backend = "cmdstanr")
  bridge <- log_marginal_likelihood(sampled, "cmdstanr", flatcor_density(), "normal", seed = 1)
  expect_within_error(bridge$log_marginal_likelihood, closed_form, bridge$error)
})

test_that("every prior's Bayes factor against the intercept-only model matches its reference", {
  skip_without_stan()
  y <- example_data$y
  n <- length(y)
  k <- 2
  r_squared <- sample_r_squared(y, cbind(example_data$x1, example_data$x2))
  for (prior in c("flat_partial_r", "flat_r", "flat_r_squared")) {
    reference <- if (prior == "flat_r") log_bayes_factor_flat_r_by_directions(r_squared, n, k) else
      log_bayes_factor_by_g_mixture(r_squared, n, k, prior_shape(prior, k))
    result <- bayes_factor_models(y ~ x1 + x2, y ~ 1, example_data, prior = prior, seed = 1)
    expect_s3_class(result, "flatcor_bayes_factor")
    expect_equal(unname(result$divergent), c(0, 0))
    expect_within_error(result$log_bayes_factor, reference, result$error)
  }
  expect_output(print(result), "flat_r_squared prior")
})

test_that("rescaling the predictors and the outcome leaves the Bayes factor unchanged", {
  skip_without_stan()
  rescaled <- transform(example_data, x1 = 100 * x1 + 7, y = 1000 * y - 3)
  original <- bayes_factor_models(y ~ x1 + x2, y ~ x2, example_data, prior = "flat_r", seed = 2)
  changed  <- bayes_factor_models(y ~ x1 + x2, y ~ x2, rescaled, prior = "flat_r", seed = 2)
  expect_within_error(changed$log_bayes_factor, original$log_bayes_factor,
                      sqrt(original$error^2 + changed$error^2))
})

test_that("a seed, or set.seed() beforehand, makes the Bayes factor reproducible", {
  skip_without_stan()
  bayes_factor_once <- function(seed = NULL) {
    bayes_factor_models(y ~ x1 + x2, y ~ x2, example_data, prior = "flat_r",
                        iter_sampling = 2000, seed = seed)$log_bayes_factor
  }
  expect_identical(bayes_factor_once(seed = 3), bayes_factor_once(seed = 3))
  set.seed(7)
  first <- bayes_factor_once()
  set.seed(7)
  expect_identical(bayes_factor_once(), first)
})
