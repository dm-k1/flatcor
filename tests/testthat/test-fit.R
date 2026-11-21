test_that("fit_regression() runs and works with print, loo and pp_check", {
  skip_on_cran()
  has_stan <- requireNamespace("cmdstanr", quietly = TRUE) ||
    requireNamespace("rstan", quietly = TRUE)
  skip_if_not(has_stan, "needs cmdstanr or rstan")
  set.seed(1)
  data <- data.frame(x1 = rnorm(50), x2 = rnorm(50))
  data$y <- 0.4 * data$x1 + rnorm(50)
  fit <- fit_regression(y ~ x1 + x2, data, chains = 2, iter_warmup = 500,
                        iter_sampling = 500, seed = 1, prior = "flat_r")
  expect_s3_class(fit, "flatcor_fit")
  expect_setequal(posterior::variables(fit$draws), c("Intercept", "x1", "x2", "sigma", "R_squared"))
  expect_output(print(fit), "flat_r prior")
  ## what the backend computed comes out too
  expect_identical(fit$seed, 1)
  expect_identical(nrow(fit$sampler), 2L)
  expect_identical(sum(fit$sampler$divergent), fit$divergent)
  expect_gt(fit$seconds, 0)
  expect_identical(fit$stan_data$N, 50L)
  stan_variables <- posterior::variables(stan_draws(fit$fit, fit$backend))
  expect_true(all(c("alpha", "R", "f", "direction[1]", "b[1]", "lp__") %in% stan_variables))
  expect_s3_class(loo(fit), "loo")
  expect_s3_class(pp_check(fit), "ggplot")
  expect_error(fit_regression(y ~ 0 + x1, data, prior = "flat_r"), "intercept")
  for (prior in c("flat_partial_r", "flat_r_squared")) {
    other <- fit_regression(y ~ x1 + x2, data, chains = 2, iter_warmup = 500,
                            iter_sampling = 500, seed = 1, prior = prior)
    expect_output(print(other), paste(prior, "prior"))
  }
  expect_error(fit_regression(y ~ x1 + x2, data, prior = "flat"), "Unknown prior")
})

test_that("cmdstanr fits leave no files in R's temporary folder, and keep working", {
  skip_without_cmdstan()
  set.seed(1)
  data <- data.frame(x1 = rnorm(50), x2 = rnorm(50))
  data$y <- 0.4 * data$x1 + rnorm(50)
  model <- flatcor_model("cmdstanr")
  ## per fit, cmdstanr writes CSV files (to flatcor's flatcor-fit-* folder)
  ## and the data (standata-*.json); a long-running session filled its disk
  left_over <- function() {
    list.files(tempdir(), pattern = "^(flatcor-fit-|standata-)", all.files = TRUE)
  }
  before <- left_over()
  fit_with_seed <- function(seed) {
    fit_regression(y ~ x1 + x2, data, prior = "flat_r", chains = 2, iter_warmup = 300,
                   iter_sampling = 300, seed = seed, backend = "cmdstanr", model = model)
  }
  fits <- lapply(c(1, 2, 1), fit_with_seed)
  expect_identical(left_over(), before)
  ## the draws and diagnostics are in memory; the same seed gives the same draws
  expect_false(any(file.exists(fits[[1]]$fit$output_files())))
  expect_identical(fits[[3]]$fit$draws(), fits[[1]]$fit$draws())
  expect_identical(fits[[3]]$draws, fits[[1]]$draws)
  expect_identical(fits[[3]]$sampler$divergent, fits[[1]]$sampler$divergent)
  expect_false(identical(fits[[2]]$draws, fits[[1]]$draws))
  expect_s3_class(fits[[1]]$fit$sampler_diagnostics(), "draws_array")
  expect_equal(fits[[1]]$fit$metadata()$max_treedepth, 10)
  expect_s3_class(loo(fits[[1]]), "loo")
})
