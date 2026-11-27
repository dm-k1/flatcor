## The rstan backend, checked mechanically: it runs, and its Bayes factor
## matches the independent reference (helper-references.R). The validation
## (flatcor-validation) uses cmdstanr only.
##
## rstan runs in its own R process: once rstan is loaded, cmdstanr can no
## longer load a model's log density (both bring their own copy of the TBB
## library), so the two backends cannot share an R session for Bayes factors.

run_with_rstan <- function(code) {
  package_folder <- normalizePath(test_path("..", ".."), mustWork = FALSE)
  callr::r(
    function(package_folder, code) {
      if (file.exists(file.path(package_folder, "DESCRIPTION"))) {
        pkgload::load_all(package_folder, quiet = TRUE)
      } else {
        library(flatcor)
      }
      eval(code)
    },
    args = list(package_folder = package_folder, code = code)
  )
}

test_that("fit_regression() and bayes_factor_models() run with rstan", {
  skip_on_cran()
  ## installed is enough: rstan is only loaded in the separate process
  skip_if_not(nzchar(system.file(package = "rstan")), "needs rstan")
  skip_if_not(requireNamespace("callr", quietly = TRUE), "needs callr")
  skip_if_cores_limited()
  set.seed(2)
  data <- data.frame(x1 = rnorm(40), x2 = rnorm(40))
  data$y <- 0.4 * data$x1 + rnorm(40)
  result <- run_with_rstan(bquote({
    data <- .(data)
    fit <- fit_regression(y ~ x1 + x2, data, prior = "flat_r", chains = 2,
                          iter_warmup = 500, iter_sampling = 1000, seed = 1, backend = "rstan")
    bayes_factor <- bayes_factor_models(y ~ x1 + x2, y ~ 1, data, prior = "joint_partial_r",
                                        iter_sampling = 4000, seed = 5, backend = "rstan")
    list(variables = posterior::variables(fit$draws), sampler = fit$sampler,
         stan_variables = posterior::variables(posterior::as_draws(as.array(fit$fit))),
         bayes_factor = bayes_factor)
  }))
  expect_setequal(result$variables, c("Intercept", "x1", "x2", "sigma", "R_squared"))
  ## the sampler's diagnostics and every Stan variable come out with rstan too
  expect_identical(nrow(result$sampler), 2L)
  expect_false(anyNA(result$sampler))
  expect_true(all(c("alpha", "R", "f", "direction[1]") %in% result$stan_variables))
  model_a <- result$bayes_factor$models$a
  expect_identical(nrow(model_a$sampler), 4L)
  expect_true(all(c("R", "rhat", "ess_tail") %in% c(model_a$summary$variable,
                                                       names(model_a$summary))))
  expect_gt(model_a$bridge_sampling$iterations, 0)
  r_squared <- sample_r_squared(data$y, cbind(data$x1, data$x2))
  reference <- log_bayes_factor_by_g_mixture(r_squared, 40, 2, prior_shape("joint_partial_r", 2))
  difference <- abs(result$bayes_factor$log_bayes_factor - reference)
  expect_lt(difference, max(5 * result$bayes_factor$error, 0.02))
})
