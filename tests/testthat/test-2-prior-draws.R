## The prior as a recipe gives R-squared ~ Beta(a, 1)
## (Objective-bayes/docs/2_the_four_priors.md)
test_that("flat R: implied R-squared is Beta(1/2, 1) for k = 1, 3, 10, correlated or not", {
  set.seed(1)
  correlated <- matrix(0.5, 3, 3)
  diag(correlated) <- 1
  for (covariance in list(matrix(1), correlated, diag(10))) {
    draws <- draw_coefficients_from_prior(50000, covariance, prior = "flat_r")
    f_squared <- rowSums((draws %*% covariance) * draws)
    r_squared <- f_squared / (1 + f_squared)
    distance <- suppressWarnings(ks.test(r_squared, "pbeta", 0.5, 1)$statistic)
    expect_lt(distance, 0.01)
  }
})

test_that("with one predictor the correlation is flat from -1 to 1", {
  set.seed(2)
  effect <- draw_coefficients_from_prior(50000, matrix(1), prior = "flat_r")[, 1]
  correlation <- effect / sqrt(1 + effect^2)
  expect_lt(suppressWarnings(ks.test(correlation, "punif", -1, 1)$statistic), 0.01)
})

test_that("joint flat partial r: implied R-squared is Beta(k/2, 1), correlated or not", {
  set.seed(3)
  correlated <- matrix(0.5, 3, 3)
  diag(correlated) <- 1
  for (covariance in list(matrix(1), correlated, diag(10))) {
    draws <- draw_coefficients_from_prior(50000, covariance, prior = "joint_partial_r")
    f_squared <- rowSums((draws %*% covariance) * draws)
    distance <- suppressWarnings(ks.test(f_squared / (1 + f_squared), "pbeta",
                                         ncol(covariance) / 2, 1)$statistic)
    expect_lt(distance, 0.01)
  }
})

test_that("joint flat partial r: every coefficient's partial correlation is flat, for any k", {
  set.seed(4)
  correlated <- matrix(0.5, 5, 5)
  diag(correlated) <- 1
  draws <- draw_coefficients_from_prior(50000, correlated, prior = "joint_partial_r")
  # effect / spread left after the others
  partial_f <- draws[, 1] * sqrt(1 / solve(correlated)[1, 1])
  partial_r <- partial_f / sqrt(1 + partial_f^2)
  expect_lt(suppressWarnings(ks.test(partial_r, "punif", -1, 1)$statistic), 0.01)
})

test_that("flat R-squared: implied R-squared is flat for every k, correlated or not", {
  set.seed(5)
  correlated <- matrix(0.5, 3, 3)
  diag(correlated) <- 1
  for (covariance in list(matrix(1), correlated, diag(10))) {
    draws <- draw_coefficients_from_prior(50000, covariance, prior = "flat_r_squared")
    f_squared <- rowSums((draws %*% covariance) * draws)
    expect_lt(suppressWarnings(ks.test(f_squared / (1 + f_squared), "punif")$statistic), 0.01)
  }
})

test_that("the default is joint_partial_r, and an unknown prior is refused", {
  set.seed(1)
  by_default <- suppressMessages(draw_coefficients_from_prior(10, diag(2)))
  set.seed(1)
  explicit <- draw_coefficients_from_prior(10, diag(2), prior = "joint_partial_r")
  expect_equal(by_default, explicit)
  expect_error(draw_coefficients_from_prior(10, diag(2), prior = "uniform_correlations"),
               "Unknown prior")
})

test_that("the prior draws carry what each set was drawn from", {
  set.seed(3)
  covariance <- matrix(c(1, 0.4, 0.4, 1), 2)
  draws <- draw_coefficients_from_prior(1000, covariance, prior = "flat_r_squared")
  r_squared <- attr(draws, "R_squared")
  f_squared <- attr(draws, "f_squared")
  expect_length(r_squared, 1000)
  expect_equal(f_squared, r_squared / (1 - r_squared))
  ## each set's f^2 is the drawn amount, and its direction the drawn one
  expect_equal(unname(rowSums((draws %*% covariance) * draws)), f_squared)
  direction <- attr(draws, "direction")
  coefficients <- matrix(c(draws), nrow(draws)) # without the attributes
  expect_equal(coefficients / sqrt(rowSums(coefficients^2)),
               direction / sqrt(rowSums(direction^2)))
})

## k = 3 predictors correlated 0.6, and 1e5 draws from each prior
correlated_three <- matrix(0.6, 3, 3)
diag(correlated_three) <- 1

test_that("each partial correlation is flat, for independent and joint partial r", {
  for (prior in c("independent_partial_r", "joint_partial_r")) {
    set.seed(11)
    draws <- draw_coefficients_from_prior(1e5, correlated_three, prior = prior)
    partial_r <- attr(draws, "partial_r")
    expect_equal(dim(partial_r), c(1e5, 3))
    for (column in 1:3) {
      expect_gt(suppressWarnings(ks.test(partial_r[, column], "punif", -1, 1)$p.value), 0.01)
    }
  }
})

test_that("the partial correlations are unlinked for independent, linked for joint", {
  set.seed(12)
  independent <- attr(
    draw_coefficients_from_prior(1e5, correlated_three, prior = "independent_partial_r"),
    "partial_r"
  )
  joint <- attr(
    draw_coefficients_from_prior(1e5, correlated_three, prior = "joint_partial_r"),
    "partial_r"
  )
  ## the size of one effect says nothing about the size of another, or it does
  expect_lt(abs(cor(abs(independent[, 1]), abs(independent[, 2]))), 0.02)
  expect_gt(cor(abs(joint[, 1]), abs(joint[, 2])), 0.05)
})

test_that("independent_partial_r carries what it was drawn from", {
  set.seed(13)
  draws <- draw_coefficients_from_prior(1000, correlated_three, prior = "independent_partial_r")
  f_squared <- attr(draws, "f_squared")
  expect_equal(unname(rowSums((draws %*% correlated_three) * draws)), f_squared)
  expect_equal(attr(draws, "R_squared"), f_squared / (1 + f_squared))
  expect_equal(attr(draws, "direction"), matrix(c(draws), 1000))
})

test_that("with one predictor, independent and joint partial r are the same prior", {
  set.seed(14)
  independent <- draw_coefficients_from_prior(1e5, matrix(1), prior = "independent_partial_r")
  joint <- draw_coefficients_from_prior(1e5, matrix(1), prior = "joint_partial_r")
  expect_gt(
    suppressWarnings(ks.test(attr(independent, "R_squared"), attr(joint, "R_squared"))$p.value),
    0.01
  )
})

test_that("Stan, sampling the prior only, gives the draws of the recipe", {
  skip_without_cmdstan()
  skip_if_cores_limited()
  set.seed(15)
  predictors <- matrix(rnorm(200 * 2), 200) %*% chol(matrix(c(1, 0.6, 0.6, 1), 2))
  for (prior in c("independent_partial_r", "joint_partial_r")) {
    stan_data <- stan_data_of(predictors, rep(0, 200), prior)
    stan_data$prior_only <- 1L
    stan_data$proper_stand_ins <- 1L
    fit <- flatcor_model()$sample(
      data = stan_data,
      chains = 4,
      parallel_chains = 4,
      iter_warmup = 1000,
      iter_sampling = 2500,
      seed = 1,
      refresh = 0,
      show_messages = FALSE
    )
    ## the coefficients in units of the noise, as the recipe draws them
    draws <- posterior::as_draws_df(fit$draws(c("b", "sigma", "R_squared")))
    scaled_b <- cbind(draws[["b[1]"]], draws[["b[2]"]]) / draws$sigma
    centered <- scale(predictors, scale = FALSE)
    recipe <- draw_coefficients_from_prior(1e5, crossprod(centered) / 200, prior = prior)
    ## the sampler's draws are autocorrelated, so the distance is allowed some room
    for (column in 1:2) {
      distance <- suppressWarnings(ks.test(scaled_b[, column], recipe[, column])$statistic)
      expect_lt(distance, 0.05)
    }
    distance <- suppressWarnings(
      ks.test(draws$R_squared, attr(recipe, "R_squared"))$statistic
    )
    expect_lt(distance, 0.05)
  }
})

test_that("partial_r is each coefficient's partial correlation, by the textbook definition", {
  ## The partial correlation of y and x_j given the other predictors, read
  ## from the inverse of the covariance matrix of (y, x1 ... xk):
  ## -P[y, j] / sqrt(P[y, y] * P[j, j]). This checks the conversion
  ## b_j = rho_j / (s_j sqrt(1 - rho_j^2)) against a definition that does not use it.
  set.seed(16)
  for (prior in prior_names) {
    draws <- draw_coefficients_from_prior(5, correlated_three, prior = prior)
    for (draw in 1:5) {
      b <- draws[draw, ]
      ## covariance of (y, x1 ... xk) with sigma = 1: Var(y) = b' Sigma b + 1
      covariance <- rbind(
        c(sum(b * (correlated_three %*% b)) + 1, b %*% correlated_three),
        cbind(correlated_three %*% b, correlated_three)
      )
      precision <- solve(covariance)
      textbook <- -precision[1, -1] / sqrt(precision[1, 1] * diag(precision)[-1])
      expect_equal(attr(draws, "partial_r")[draw, ], textbook, tolerance = 1e-8)
    }
  }
})
