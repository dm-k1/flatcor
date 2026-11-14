## Fast versions of objective-bayes, study 2 (are the Bayes factors computed correctly?)

example_data <- local({
  set.seed(42)
  correlations <- matrix(0.4, 3, 3); diag(correlations) <- 1
  predictors <- matrix(rnorm(60 * 3), 60) %*% chol(correlations)
  data <- data.frame(x1 = predictors[, 1], x2 = predictors[, 2], x3 = predictors[, 3])
  data$y <- 0.3 * data$x1 + 0.2 * data$x2 + rnorm(60)
  data
})

test_that("the example's Bayes factors for each predictor are unchanged (model against model)", {
  full <- lm(y ~ x1 + x2 + x3, example_data)
  without <- function(term) lm(update(formula(full), paste(". ~ . -", term)), example_data)
  bayes_factors <- sapply(c("x1", "x2", "x3"), function(term)
    bayes_factor_models(full, without(term), prior = "flat_r"))
  expect_equal(unname(bayes_factors), c(18.696, 0.3827, 1.5872), tolerance = 1e-3)
})

test_that("the Bayes factor for one predictor depends on what the others explain", {
  ## same evidence about x1 (partial correlation 0.2, n = 100, k = 3); the
  ## other predictors explain 5% or 60% (objective-bayes/docs/2_how_it_calculates.md, section 10)
  joint_bf <- function(r2_without) {
    r2_with <- 1 - (1 - r2_without) * (1 - 0.2^2)
    exp(log_bayes_factor_flat_r_vs_intercept_only(r2_with, 100, 3) -
        log_bayes_factor_flat_r_vs_intercept_only(r2_without, 100, 2))
  }
  expect_gt(joint_bf(0.05), 2)     # leans toward an effect
  expect_lt(joint_bf(0.60), 1)     # leans toward no effect, with the same x1 evidence
})

test_that("the likelihood ratio of q columns is exp(-n f^2/2) 1F1((m+1)/2; q/2; n r^2 f^2/2)", {
  log_kummer <- function(a, c, z) {
    j <- 0:3000
    terms <- lgamma(a + j) - lgamma(a) - lgamma(c + j) + lgamma(c) - lgamma(j + 1) + j * log(z)
    max(terms) + log(sum(exp(terms - max(terms))))
  }
  for (case in list(c(40, 4, 2, 0.1, 0.3), c(200, 6, 3, 0.05, 0.2), c(100, 5, 5, 0.3, 0.5))) {
    n <- case[1]; k <- case[2]; q <- case[3]; r_squared <- case[4]; f <- case[5]
    m <- n - k + q - 2
    expect_equal(log_likelihood_ratio_block(f, r_squared, n, k, q),
                 -n * f^2 / 2 + log_kummer((m + 1) / 2, q / 2, n * r_squared * f^2 / 2),
                 tolerance = 1e-8)
  }
})

test_that("an aliased column stops with a clear message", {
  aliased <- transform(example_data, x4 = x1 + x2)
  expect_error(bayes_factor_models(lm(y ~ x1 + x2 + x4, aliased), lm(y ~ 1, aliased), prior = "flat_r"),
               "linearly dependent")
})

test_that("rescaling or recoding predictors leaves the Bayes factor unchanged", {
  original <- bayes_factor_models(lm(y ~ x1 + x2 + x3, example_data), lm(y ~ 1, example_data),
                                  prior = "flat_r")
  rescaled <- transform(example_data, x1 = 100 * x1 + 7, y = 1000 * y - 3)
  recoded  <- transform(example_data, sum_12 = x1 + x2, difference_12 = x1 - x2)
  expect_equal(bayes_factor_models(lm(y ~ x1 + x2 + x3, rescaled), lm(y ~ 1, rescaled),
                                   prior = "flat_r"), original)
  expect_equal(bayes_factor_models(lm(y ~ sum_12 + difference_12 + x3, recoded),
                                   lm(y ~ 1, recoded), prior = "flat_r"), original)
})

test_that("the grid agrees with adaptive integration (one predictor)", {
  ## one predictor: the g-prior Bayes factor averaged over g ~ Inverse-Gamma(1, n/2)
  r_squared <- 0.09; n <- 40
  by_integration <- integrate(function(g)
    (1 + g)^((n - 2) / 2) * (1 + g * (1 - r_squared))^(-(n - 1) / 2) *
      (n / 2) * g^(-2) * exp(-n / (2 * g)), 0, Inf)$value
  expect_equal(log_bayes_factor_flat_r_vs_intercept_only(r_squared, n, 1), log(by_integration),
               tolerance = 1e-4)
})

test_that("predictive matching: with n = k + 1 the Bayes factor is 1", {
  for (k in 1:5)
    expect_equal(exp(log_bayes_factor_flat_r_vs_intercept_only(1 - 1e-12, k + 1, k)), 1,
                 tolerance = 1e-3)
})

test_that("information consistency: evidence grows without limit as R-squared goes to 1", {
  evidence <- sapply(c(0.9, 0.99, 0.999, 0.9999), log_bayes_factor_flat_r_vs_intercept_only,
                     sample_size = 20, number_of_predictors = 3)
  expect_true(all(diff(evidence) > 0))
  expect_gt(evidence[4], 50)
})

test_that("models must share outcome and rows", {
  expect_error(bayes_factor_models(lm(y ~ x1, example_data), lm(y ~ x1, example_data[-1, ]),
                                   prior = "flat_r"), "same outcome")
})

## Flat partial r (R-squared ~ Beta(k/2, 1)): fast versions of objective-bayes/validation/02

test_that("flat partial r: the grid agrees with adaptive integration", {
  r_squared <- 0.2; n <- 50; k <- 3
  by_integration <- integrate(function(g)
    (1 + g)^((n - 1 - k) / 2) * (1 + g * (1 - r_squared))^(-(n - 1) / 2) *
      (n / 2) * g^(-2) * exp(-n / (2 * g)), 0, Inf)$value
  expect_equal(log_bayes_factor_flat_partial_r_vs_intercept_only(r_squared, n, k),
               log(by_integration), tolerance = 1e-4)
})

test_that("with one predictor flat partial r and flat R are the same", {
  expect_equal(log_bayes_factor_flat_partial_r_vs_intercept_only(0.09, 40, 1),
               log_bayes_factor_flat_r_vs_intercept_only(0.09, 40, 1), tolerance = 1e-4)
})

test_that("flat partial r: predictive matching and invariance", {
  for (k in 1:5)
    expect_equal(exp(log_bayes_factor_flat_partial_r_vs_intercept_only(1 - 1e-12, k + 1, k)), 1,
                 tolerance = 1e-3)
  original <- bayes_factor_models(lm(y ~ x1 + x2 + x3, example_data), lm(y ~ 1, example_data),
                                  prior = "flat_partial_r")
  recoded  <- transform(example_data, sum_12 = x1 + x2, difference_12 = x1 - x2)
  expect_equal(bayes_factor_models(lm(y ~ sum_12 + difference_12 + x3, recoded), lm(y ~ 1, recoded),
                                   prior = "flat_partial_r"), original)
})

test_that("a model without an intercept stops with a clear message", {
  expect_error(bayes_factor_models(lm(y ~ 0 + x1, example_data), lm(y ~ 1, example_data),
                                   prior = "flat_r"), "need an intercept")
})

test_that("the one route agrees with adaptive integration where the old grid did not", {
  ## flat partial r: the g-prior Bayes factor averaged over g ~ IG(1, n/2), by
  ## integrate() over log g. At k = 20 and R-squared = .001 the fixed grid of
  ## 2.0.0's first draft was off by 5e-4.
  by_integration <- function(r_squared, n, k) {
    log_integrand <- function(log_g) {
      g <- exp(log_g)
      (n - 1 - k) / 2 * log1p(g) - (n - 1) / 2 * log1p(g * (1 - r_squared)) +
        log(n / 2) - log_g - n / (2 * g)
    }
    peak <- optimize(log_integrand, c(-40, 60), maximum = TRUE)$maximum
    top <- log_integrand(peak)
    breaks <- peak + c(-80, -20, -5, -1, 0, 1, 5, 20, 80)
    top + log(sum(sapply(1:8, function(i)
      integrate(function(x) exp(log_integrand(x) - top), breaks[i], breaks[i + 1],
                rel.tol = 1e-12, subdivisions = 5000)$value)))
  }
  for (case in list(c(0.001, 1000, 20), c(0.3, 200, 8), c(0.5, 5, 3)))
    expect_equal(log_bayes_factor_vs_intercept_only(case[1], case[2], case[3], prior = "flat_partial_r"),
                 by_integration(case[1], case[2], case[3]), tolerance = 1e-6)
})

test_that("one residual degree of freedom works for every prior", {
  for (prior in c("flat_partial_r", "flat_r", "flat_r_squared"))
    expect_true(is.finite(log_bayes_factor_vs_intercept_only(0.5, 3, 1, prior = prior)))
})
