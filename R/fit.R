## =====================================================================
##  The Stan model: estimation (fit_regression()) and the fits behind
##  every Bayes factor (bayes_factor_models())
## =====================================================================
## The Stan model samples the prior in its own coordinates: R and a
## direction (uniform), with the coefficients computed from them. Sampling
## the coefficients directly misses the prior's peak at zero for larger k
## (flatcor-validation, studies 1a and 2a). With no predictors it is the
## intercept-only model.

compiled_models <- new.env()   # compile once per session

## The Stan model ships in inst/stan. For estimation, cmdstanr writes its
## executable to a user cache folder, one per flatcor version, because the
## installed package folder may not be writable and an executable from another
## version must never be reused. Bayes factors also need the model methods
## (the log density on the unconstrained scale, for bridgesampling), and
## cmdstanr can only provide those in the R session that compiled the model:
## that version is compiled once per session into the session's temporary
## folder, so that parallel R processes never compile into the same file.
compile_stan_model <- function(backend, model_methods = FALSE) {
  key <- paste(backend, if (model_methods && backend == "cmdstanr") "methods" else "plain")
  if (is.null(compiled_models[[key]])) {
    stan_file <- system.file("stan", "flatcor.stan", package = "flatcor", mustWork = TRUE)
    compiled_models[[key]] <- if (backend == "cmdstanr") {
      executable_folder <- if (model_methods) file.path(tempdir(), "flatcor_stan") else
        file.path(tools::R_user_dir("flatcor", which = "cache"), as.character(utils::packageVersion("flatcor")))
      dir.create(executable_folder, recursive = TRUE, showWarnings = FALSE)
      cmdstanr::cmdstan_model(stan_file, dir = executable_folder, compile_model_methods = model_methods,
                              quiet = TRUE)
    } else {
      rstan::stan_model(stan_file)
    }
  }
  compiled_models[[key]]
}

default_backend <- function()
  if (requireNamespace("cmdstanr", quietly = TRUE)) "cmdstanr" else "rstan"

## Sample the posterior of one model: predictors (a matrix, possibly with no
## columns) and outcome. Returns the backend's fit and its number of
## divergent transitions.
sample_posterior <- function(predictors, outcome, prior, chains, iter_warmup, iter_sampling,
                             adapt_delta, seed, backend, model_methods = FALSE) {
  check_full_rank(predictors)
  stan_data <- list(N = length(outcome), K = ncol(predictors), X = predictors, y = outcome,
                    prior_only = 0L, proper_stand_ins = 0L,
                    r_squared_shape = prior_shape(prior, ncol(predictors)))
  model <- compile_stan_model(backend, model_methods)
  if (backend == "cmdstanr") {
    fit <- model$sample(data = stan_data, chains = chains, parallel_chains = chains,
                        iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                        adapt_delta = adapt_delta, seed = seed, refresh = 0,
                        show_messages = FALSE)
    divergent <- sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)
  } else {
    fit <- rstan::sampling(model, data = stan_data, chains = chains, cores = chains,
                           iter = iter_warmup + iter_sampling, warmup = iter_warmup,
                           control = list(adapt_delta = adapt_delta), seed = seed, refresh = 0)
    divergent <- rstan::get_num_divergent(fit)
  }
  if (divergent > 0) warning(divergent, " divergent transitions; check the fit.", call. = FALSE)
  list(fit = fit, divergent = divergent)
}

## A column that is a combination of the others (aliased, or constant) would
## make the predictors' covariance singular: stop with a clear message instead.
check_full_rank <- function(predictors) {
  if (ncol(predictors) == 0) return(invisible(TRUE))
  centered <- scale(predictors, scale = FALSE)
  if (qr(centered)$rank < ncol(predictors))
    stop("The predictor columns are linearly dependent (one is a combination of ",
         "the others, or constant). Remove the redundant columns.", call. = FALSE)
  invisible(TRUE)
}

## lm() and model.frame() drop every row with a missing value in a variable of
## the formula. The sample size n enters every Bayes factor, so say so.
warn_dropped_rows <- function(data, rows_used) {
  rows_dropped <- nrow(data) - rows_used
  if (rows_dropped > 0)
    warning(sprintf(paste("%d of %d rows have a missing value in a variable of the formula",
                          "and were left out; the analysis uses n = %d."),
                    rows_dropped, nrow(data), rows_used), call. = FALSE)
}

#' Fit a linear regression with one of the three priors
#'
#' Samples the prior in its own coordinates, the overall correlation R
#' and a direction (uniform), with the coefficients computed from them
#' (\code{inst/stan/flatcor.stan}). Intercept flat, noise level \eqn{1/\sigma}.
#'
#' @param formula A model formula with an intercept and at least one predictor.
#' @param data A data frame.
#' @param chains,iter_warmup,iter_sampling,adapt_delta,seed Sampler settings.
#' @param backend \code{"cmdstanr"} (default when installed) or \code{"rstan"}.
#' @param prior \code{"flat_partial_r"} (the default; R-squared ~ Beta(k/2, 1)),
#'   \code{"flat_r"} (R flat) or \code{"flat_r_squared"} (R-squared flat).
#' @return A \code{flatcor_fit}: \code{print()} shows estimates and diagnostics,
#'   \code{loo()} and \code{pp_check()} work as usual, and \code{$draws}
#'   holds all posterior draws as a \code{posterior::draws_df}.
#' @examples
#' \dontrun{
#' fit <- fit_regression(mpg ~ wt + hp, mtcars, prior = "flat_partial_r")
#' fit
#' loo(fit)
#' }
#' @export
fit_regression <- function(formula, data, prior = "flat_partial_r", chains = 4, iter_warmup = 1000,
                           iter_sampling = 1000, adapt_delta = 0.95, seed = NULL,
                           backend = default_backend()) {
  prior <- check_prior(prior, defaulted = missing(prior))
  model_frame <- model.frame(formula, data)
  warn_dropped_rows(data, nrow(model_frame))
  design      <- model.matrix(formula, model_frame)
  if (!"(Intercept)" %in% colnames(design))
    stop("fit_regression() needs a model with an intercept.")
  predictors <- design[, colnames(design) != "(Intercept)", drop = FALSE]
  if (ncol(predictors) == 0)
    stop("fit_regression() needs at least one predictor.")
  outcome <- as.numeric(model.response(model_frame))
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)

  sampled <- sample_posterior(predictors, outcome, prior, chains, iter_warmup, iter_sampling,
                              adapt_delta, seed, backend)
  all_draws <- if (backend == "cmdstanr") posterior::as_draws_df(sampled$fit$draws()) else
                 posterior::as_draws_df(as.array(sampled$fit))

  ## Coefficients under their own names; the intercept back on the data's scale
  ## (the model's alpha is the mean of y, at the predictors' means)
  coefficients <- unclass(posterior::as_draws_matrix(posterior::subset_draws(all_draws, variable = "b")))
  colnames(coefficients) <- colnames(predictors)
  intercept <- all_draws$alpha - drop(coefficients %*% colMeans(predictors))
  draws <- posterior::as_draws_df(cbind(
    data.frame(Intercept = intercept, check.names = FALSE),
    as.data.frame(coefficients, check.names = FALSE),
    data.frame(sigma = all_draws$sigma, R_squared = all_draws$R_squared)))

  structure(list(draws        = draws,
                 log_lik      = unclass(posterior::as_draws_matrix(posterior::subset_draws(all_draws, variable = "log_lik"))),
                 y_rep        = unclass(posterior::as_draws_matrix(posterior::subset_draws(all_draws, variable = "y_rep"))),
                 chain_ids    = all_draws$.chain,
                 outcome      = outcome,
                 predictors   = predictors,
                 formula      = formula,
                 prior        = prior,
                 divergent    = sampled$divergent),
            class = "flatcor_fit")
}

#' @export
print.flatcor_fit <- function(x, digits = 3, ...) {
  cat("Linear regression with the", x$prior, "prior\n")
  cat("Formula:", deparse(x$formula), "\n")
  cat("Observations:", length(x$outcome), "  Predictors:", ncol(x$predictors),
      "  Divergent transitions:", x$divergent, "\n\n")
  estimates <- posterior::summarise_draws(
    x$draws, mean, sd,
    ~ quantile(.x, probs = c(0.025, 0.975)),
    rhat = posterior::rhat, ess_bulk = posterior::ess_bulk)
  print(as.data.frame(estimates), digits = digits, row.names = FALSE)
  cat("\nR_squared is the population explained share; its prior is",
      sprintf("Beta(%g, 1).\n", prior_shape(x$prior, ncol(x$predictors))))
  invisible(x)
}

#' Leave-one-out cross-validation for a flatcor fit
#'
#' Runs on one core by default: fast enough here, and forking the R session is
#' not allowed in Positron or RStudio.
#' @param x A \code{flatcor_fit}.
#' @param cores Number of cores.
#' @param ... Passed to \code{loo::loo()}.
#' @exportS3Method loo::loo
loo.flatcor_fit <- function(x, cores = 1, ...) {
  relative_efficiency <- loo::relative_eff(exp(x$log_lik), chain_id = x$chain_ids, cores = cores)
  loo::loo(x$log_lik, r_eff = relative_efficiency, cores = cores, ...)
}
#' Posterior predictive check for a flatcor fit
#' @param object A \code{flatcor_fit}.
#' @param ndraws Number of predicted datasets to draw.
#' @param ... Passed to \code{bayesplot::ppc_dens_overlay()}.
#' @exportS3Method bayesplot::pp_check
pp_check.flatcor_fit <- function(object, ndraws = 50, ...) {
  shown <- sample(nrow(object$y_rep), min(ndraws, nrow(object$y_rep)))
  bayesplot::ppc_dens_overlay(object$outcome, object$y_rep[shown, , drop = FALSE], ...)
}

#' @importFrom loo loo
#' @export
loo::loo

#' @importFrom bayesplot pp_check
#' @export
bayesplot::pp_check
