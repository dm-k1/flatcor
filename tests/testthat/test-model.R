## flatcor_model(), flatcor_density(), flatcor_stan_file(), flatcor_stan_data():
## what they give (no Stan needed), that the functions give the same result
## whether they build the model themselves or are given it, and that the Bayes
## factor from BridgeStan's log density is identical to the one from
## cmdstanr's model methods, which flatcor used before (version 2.0).

model_data <- local({
  set.seed(3)
  data <- data.frame(x1 = rnorm(50), x2 = rnorm(50), x3 = rnorm(50))
  data$y <- 0.4 * data$x1 + rnorm(50)
  data
})

test_that("the Stan program and its data are what the functions give Stan", {
  expect_true(file.exists(flatcor_stan_file()))
  stan_data <- flatcor_stan_data(y ~ x1 + x2, model_data, prior = "flat_r_squared")
  parts <- predictors_and_outcome(y ~ x1 + x2, model_data)
  expect_identical(stan_data, stan_data_of(parts$predictors, parts$outcome, "flat_r_squared"))
  expect_identical(stan_data$K, 2L)
  expect_identical(stan_data$r_squared_shape, 1)
  expect_identical(flatcor_stan_data(y ~ x1, model_data, prior = "flat_r", TRUE)$per_observation, 1L)
  expect_error(flatcor_stan_data(y ~ x1, model_data, prior = "flat"), "Unknown prior")
})

test_that("BridgeStan's parameter names become cmdstanr's", {
  expect_identical(
    stan_variable_names(c("alpha", "R_size.1", "z.12", "a.1.2")),
    c("alpha", "R_size[1]", "z[12]", "a[1,2]")
  )
})

test_that("a model of the wrong kind is refused", {
  expect_error(check_model("not a model", "cmdstanr"), "must be a cmdstanr model")
  expect_error(check_model(list(), "rstan"), "must be a rstan model")
})

test_that("given the model, or building it themselves, the functions give the same result", {
  skip_on_cran()
  skip_if_not(requireNamespace("cmdstanr", quietly = TRUE), "needs cmdstanr")
  model <- flatcor_model("cmdstanr")
  expect_s3_class(model, "CmdStanModel")
  expect_identical(flatcor_model("cmdstanr"), model)
  given <- fit_regression(y ~ x1 + x2, model_data, prior = "flat_r", chains = 2,
                          iter_warmup = 500, iter_sampling = 500, seed = 4, model = model)
  built <- fit_regression(y ~ x1 + x2, model_data, prior = "flat_r", chains = 2,
                          iter_warmup = 500, iter_sampling = 500, seed = 4)
  expect_identical(given$draws, built$draws)

  skip_if_not(requireNamespace("bridgestan", quietly = TRUE), "needs bridgestan")
  density <- flatcor_density()
  expect_s3_class(density, "flatcor_density")
  expect_true(file.exists(density$library))
  given <- bayes_factor_models(y ~ x1 + x2, y ~ x2, model_data, prior = "flat_r",
                               iter_sampling = 1000, seed = 5, model = model, density = density)
  built <- bayes_factor_models(y ~ x1 + x2, y ~ x2, model_data, prior = "flat_r",
                               iter_sampling = 1000, seed = 5)
  expect_identical(given$log_bayes_factor, built$log_bayes_factor)
})

test_that("the cache folder can be chosen, and the builds are found there", {
  skip_on_cran()
  skip_if_not(requireNamespace("cmdstanr", quietly = TRUE), "needs cmdstanr")
  cache <- file.path(tempdir(), "flatcor_cache_test")
  model <- flatcor_model("cmdstanr", cache = cache)
  expect_true(startsWith(normalizePath(model$exe_file()), normalizePath(cache)))
  skip_if_not(requireNamespace("bridgestan", quietly = TRUE), "needs bridgestan")
  density <- flatcor_density(cache = cache)
  expect_true(startsWith(normalizePath(density$library), normalizePath(cache)))
})

test_that("the Bayes factor from BridgeStan is identical to the one from model methods", {
  skip_on_cran()
  skip_if_not(requireNamespace("cmdstanr", quietly = TRUE), "needs cmdstanr")
  skip_if_not(requireNamespace("bridgestan", quietly = TRUE), "needs bridgestan")
  ## flatcor 2.0's route: the model compiled with cmdstanr's model methods,
  ## and bridgesampling's method for cmdstanr fits, which uses them
  methods_folder <- file.path(tempdir(), "flatcor_methods")
  dir.create(methods_folder, showWarnings = FALSE)
  with_methods <- cmdstanr::cmdstan_model(
    flatcor_stan_file(),
    dir = methods_folder,
    compile_model_methods = TRUE,
    force_recompile = TRUE,
    quiet = TRUE
  )
  pairs <- list(
    list(formulas = list(y ~ x1 + x2, y ~ x2), prior = "flat_partial_r"),
    list(formulas = list(y ~ x1 + x2, y ~ 1), prior = "flat_r"),
    list(formulas = list(y ~ x1 + x2 + x3, y ~ x2 + x3), prior = "flat_r_squared")
  )
  for (pair in pairs) {
    rows <- complete_rows(model_data, pair$formulas)
    by_methods <- vapply(c(0, 1), function(offset) {
      parts <- predictors_and_outcome(pair$formulas[[offset + 1]], rows)
      sampled <- sample_posterior(parts$predictors, parts$outcome, pair$prior, chains = 4,
                                  iter_warmup = 1000, iter_sampling = 1000, adapt_delta = 0.95,
                                  seed = 6 + offset, backend = "cmdstanr", model = with_methods)
      with_seed(6 + offset, bridgesampling::bridge_sampler(sampled$fit, silent = TRUE))$logml
    }, numeric(1))
    by_bridgestan <- bayes_factor_models(pair$formulas[[1]], pair$formulas[[2]], model_data,
                                         prior = pair$prior, iter_sampling = 1000, seed = 6)
    expect_identical(by_bridgestan$log_bayes_factor, by_methods[1] - by_methods[2])
  }
})
