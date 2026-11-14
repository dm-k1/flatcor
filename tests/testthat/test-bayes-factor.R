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
  ## other predictors explain 5% or 60% (objective-bayes/docs/appendix/A_derivations.md, section 10)
  joint_bf <- function(r2_without) {
    r2_with <- 1 - (1 - r2_without) * (1 - 0.2^2)
    exp(log_bayes_factor_vs_intercept_only(r2_with, 100, 3, prior = "flat_r") -
        log_bayes_factor_vs_intercept_only(r2_without, 100, 2, prior = "flat_r"))
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

test_that("every prior agrees with its g-prior mixture, by integrate()", {
  for (prior in c("flat_partial_r", "flat_r", "flat_r_squared"))
    for (case in list(c(0.09, 40, 1), c(0.2, 50, 3), c(0.05, 200, 5), c(0.5, 100, 8),
                      c(0.001, 1000, 20), c(0.5, 5, 3))) {
      r_squared <- case[1]; n <- case[2]; k <- case[3]; a <- prior_shape(prior, k)
      if (a > k / 2) next   # flat R-squared with one predictor: no g-mixture (test-three-priors.R)
      reference <- if (prior == "flat_r") log_bayes_factor_flat_r_by_directions(r_squared, n, k) else
                   log_bayes_factor_by_g_mixture(r_squared, n, k, a)
      expect_equal(log_bayes_factor_vs_intercept_only(r_squared, n, k, prior = prior), reference,
                   tolerance = 1e-6)
    }
})

test_that("predictive matching: with n = k + 1 the Bayes factor is 1, for every prior", {
  for (prior in c("flat_partial_r", "flat_r", "flat_r_squared"))
    for (k in 1:5)
      expect_equal(log_bayes_factor_vs_intercept_only(1, k + 1, k, prior = prior), 0)
})

test_that("information consistency: evidence grows without limit as R-squared goes to 1", {
  for (prior in c("flat_partial_r", "flat_r", "flat_r_squared")) {
    evidence <- sapply(c(0.9, 0.99, 0.999, 0.9999), log_bayes_factor_vs_intercept_only,
                       sample_size = 20, number_of_predictors = 3, prior = prior)
    expect_true(all(diff(evidence) > 0))
    expect_gt(evidence[4], 50)
  }
})

test_that("models must share outcome and rows", {
  expect_error(bayes_factor_models(lm(y ~ x1, example_data), lm(y ~ x1, example_data[-1, ]),
                                   prior = "flat_r"), "same outcome")
})

test_that("with one predictor flat partial r and flat R are the same", {
  expect_equal(log_bayes_factor_vs_intercept_only(0.09, 40, 1, prior = "flat_partial_r"),
               log_bayes_factor_vs_intercept_only(0.09, 40, 1, prior = "flat_r"))
})

test_that("flat partial r: recoding the predictors leaves the Bayes factor unchanged", {
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

test_that("one residual degree of freedom works for every prior", {
  for (prior in c("flat_partial_r", "flat_r", "flat_r_squared"))
    expect_true(is.finite(log_bayes_factor_vs_intercept_only(0.5, 3, 1, prior = prior)))
})
