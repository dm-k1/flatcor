## =====================================================================
##  The Bayes factor between two models, for the three priors
## =====================================================================
## Each model gets the Stan model (R/fit.R, inst/stan/flatcor.stan) with the
## chosen prior; its marginal likelihood is estimated from the posterior
## draws by bridge sampling (bridgesampling; Gronau, Singmann & Wagenmakers,
## 2020), and the Bayes factor is their ratio. Nothing is computed by a
## formula of ours.
##
## Following the bridgesampling authors:
##   - testing needs about ten times the draws of estimation, so the default
##     is 4 chains of 10,000 draws;
##   - the error reported is bridgesampling's approximate error (the
##     coefficient of variation of each marginal likelihood, so the SE of its
##     log), combined over the two models; it is available for method =
##     "normal" only. Their "gold standard" is to repeat the whole Bayes
##     factor with fresh posterior draws and look at the spread;
##   - method = "warp3" is at least as precise and takes about twice as long;
##     it helps when the posterior is skewed.
## The intercept and the noise have improper priors (flat and 1/sigma), the
## same in both models, so their arbitrary constant cancels in the ratio.
##
## Bridge sampling evaluates each model's log density at many points on
## Stan's unconstrained scale. With cmdstanr that log density comes from
## BridgeStan's library (flatcor_density()), built once and loaded with each
## model's data; with rstan, from rstan's fit itself. Both are Stan's log
## density up to a constant, with the Jacobian (what cmdstanr's and rstan's
## log_prob give), which is enough here because flatcor.stan writes every
## density with all its normalising constants.

## The log marginal likelihood of one sampled model (sample_posterior()) and
## its approximate error, by bridge sampling, with bridge sampling's whole
## result ($bridge) and the seconds it took. seed: of R's random numbers,
## which bridge sampling uses (put back afterwards).
log_marginal_likelihood <- function(sampled, backend, density, method, seed) {
  started <- proc.time()[["elapsed"]]
  bridge <- if (backend == "rstan") {
    with_seed(seed, bridgesampling::bridge_sampler(sampled$fit, method = method, silent = TRUE))
  } else {
    if (!inherits(density, "flatcor_density")) {
      stop("density = must be flatcor_density().", call. = FALSE)
    }
    # the log density: BridgeStan's library, loaded with this model's data
    data_file <- tempfile(fileext = ".json")
    on.exit(unlink(data_file))
    cmdstanr::write_stan_json(sampled$stan_data, data_file)
    log_density <- bridgestan::StanModel$new(density$library, data_file, seed = 1, warn = FALSE)

    # the posterior draws on Stan's unconstrained scale
    draws <- posterior::as_draws_matrix(
      sampled$fit$draws(stan_variable_names(log_density$param_names()))
    )
    unconstrained <- t(apply(unclass(draws), 1, log_density$param_unconstrain))
    colnames(unconstrained) <- log_density$param_unc_names()
    bound <- stats::setNames(rep(Inf, ncol(unconstrained)), colnames(unconstrained))

    with_seed(seed, bridgesampling::bridge_sampler(
      unconstrained,
      log_posterior = log_posterior_bridgestan,
      data = log_density,
      lb = -bound,
      ub = bound,
      method = method,
      silent = TRUE
    ))
  }
  seconds <- proc.time()[["elapsed"]] - started
  error_measures <- if (method == "normal") bridgesampling::error_measures(bridge)
  list(
    log_marginal_likelihood = bridge$logml,
    error = if (method == "normal") error_measures$cv else NA_real_,
    error_measures = error_measures,
    bridge = bridge,
    seconds = seconds
  )
}

## Everything about one model of a Bayes factor: its marginal likelihood and
## bridge sampling, the sampler's diagnostics, a summary of every Stan
## variable, the timings and the seed. keep_fits: also the Stan fit (read
## into memory), bridge sampling's object and the data Stan got.
model_details <- function(sampled, bridged, backend, keep_fits) {
  details <- list(
    log_marginal_likelihood = bridged$log_marginal_likelihood,
    error = bridged$error,
    bridge_sampling = list(
      method = bridged$bridge$method,
      iterations = bridged$bridge$niter,
      repetitions = bridged$bridge$repetitions,
      error_measures = bridged$error_measures
    ),
    divergent = sampled$divergent,
    sampler = sampler_by_chain(sampled$fit, backend),
    summary = summarise_stan_draws(stan_draws(sampled$fit, backend)),
    seconds = c(sampling = sampled$seconds, bridge_sampling = bridged$seconds),
    seed = sampled$seed,
    predictors = sampled$stan_data$K,
    n = sampled$stan_data$N
  )
  if (keep_fits) {
    details$fit <- fit_in_memory(sampled$fit, backend)
    details$bridge_sampler <- bridged$bridge
    details$stan_data <- sampled$stan_data
  }
  details
}

## The log density at one draw (unconstrained), for bridgesampling, which
## calls it with these argument names
log_posterior_bridgestan <- function(s.row, data) { # nolint: object_name_linter.
  tryCatch(
    data$log_density(as.numeric(s.row), propto = TRUE, jacobian = TRUE),
    error = function(e) -Inf
  )
}

## BridgeStan names the elements of a parameter z.1 (and a.1.2); cmdstanr
## z[1] (and a[1,2])
stan_variable_names <- function(bridgestan_names) {
  vapply(strsplit(bridgestan_names, ".", fixed = TRUE), function(parts) {
    if (length(parts) == 1) parts else
      sprintf("%s[%s]", parts[1], paste(parts[-1], collapse = ","))
  }, character(1))
}

#' Bayes factor between any two linear models
#'
#' Nested or not, as long as both have the same outcome. Each model is fitted
#' with the chosen prior (\code{inst/stan/flatcor.stan}), its marginal
#' likelihood is estimated by bridge sampling
#' (\code{bridgesampling::bridge_sampler()}), and the Bayes factor is their
#' ratio. A model with no predictors (\code{y ~ 1}) is the intercept-only
#' model. Both models use the same rows: those with no missing value in any
#' variable of either formula.
#'
#' The estimate has a Monte Carlo error. The reported error is
#' bridgesampling's approximate error of the log Bayes factor (method
#' \code{"normal"} only); it can be optimistic. To assess the error directly,
#' run the function again with other seeds and look at the spread, which the
#' bridgesampling authors recommend where it is affordable. Their largest
#' worked example has 96 parameters (Gronau et al., 2017); a model here has
#' k + 3 (k predictor columns: the intercept, the noise, R and a direction of
#' k numbers), so up to about 90 predictors stays within what has been shown.
#' @param formula_a,formula_b The two models, as formulas with an intercept
#'   and the same outcome, such as \code{y ~ x1 + x2} and \code{y ~ x2}.
#' @inheritParams fit_regression
#' @param chains,iter_warmup,iter_sampling,adapt_delta Sampler settings for
#'   each model (the default draws are ten times those for estimation).
#' @param seed Makes the result reproducible: Stan's seed for the first model
#'   (the second uses \code{seed + 1}), and the seed of R's random numbers that
#'   bridge sampling uses (R's own random state is put back afterwards).
#'   Without it the seed is drawn from R's random numbers, so
#'   \code{set.seed()} beforehand does the same.
#' @param method \code{"normal"} (the default) or \code{"warp3"}, passed to
#'   \code{bridgesampling::bridge_sampler()}.
#' @param model The model to sample, for both models:
#'   \code{flatcor_model(backend)}.
#' @param density With cmdstanr, the log density bridge sampling evaluates:
#'   \code{flatcor_density()} (needs the bridgestan package). Not used with
#'   rstan. See \code{flatcor_model()} for parallel processes.
#' @param keep_fits Also return each model's Stan fit (with all draws, read
#'   into memory), bridge sampling's object and the data Stan got. Off by
#'   default: the fits hold every draw, too much to keep for many Bayes
#'   factors; the summaries below come out either way.
#' @return A \code{flatcor_bayes_factor}:
#'   \describe{
#'     \item{\code{$bayes_factor}, \code{$log_bayes_factor}}{Model A against
#'       model B.}
#'     \item{\code{$error}}{Approximate SE of the log Bayes factor, the two
#'       models' errors combined (\code{NA} for \code{"warp3"}).}
#'     \item{\code{$log_marginal_likelihoods}, \code{$divergent}}{Per model.}
#'     \item{\code{$models}}{Everything per model (\code{$a} and \code{$b}):
#'       the log marginal likelihood and its own error; bridge sampling's
#'       method, iterations and error measures; \code{$sampler}, the
#'       sampler's diagnostics per chain (divergent transitions, hits of the
#'       maximum tree depth, E-BFMI, step size, mean leapfrog steps, seconds
#'       of warmup and sampling); \code{$summary}, every Stan variable's
#'       mean, median, sd, mad, quantiles, R-hat and bulk and tail ESS;
#'       \code{$seconds} of sampling and of bridge sampling; the seed, the
#'       number of predictors and n. With \code{keep_fits = TRUE} also
#'       \code{$fit}, \code{$bridge_sampler} and \code{$stan_data}.}
#'     \item{\code{$seed}, \code{$prior}, \code{$method}, \code{$draws},
#'       \code{$predictors}, \code{$backend}}{The seed used (drawn when none
#'       was given), and the settings.}
#'   }
#' @examples
#' # needs CmdStan (cmdstanr::install_cmdstan())
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
#'   bayes_factor_models(mpg ~ wt + hp, mpg ~ wt, mtcars, prior = "flat_partial_r", seed = 1)
#'   bayes_factor_models(mpg ~ wt + hp, mpg ~ 1, mtcars, prior = "flat_r", seed = 1)
#' }
#' @export
bayes_factor_models <- function(formula_a, formula_b, data, prior = "flat_partial_r", chains = 4,
                                iter_warmup = 1000, iter_sampling = 10000, adapt_delta = 0.95,
                                seed = NULL, method = "normal", backend = default_backend(),
                                model = flatcor_model(backend),
                                density = if (backend == "cmdstanr") flatcor_density(),
                                keep_fits = FALSE) {
  prior  <- check_prior(prior, defaulted = missing(prior))
  method <- match.arg(method, c("normal", "warp3"))
  rows <- complete_rows(data, list(formula_a, formula_b))
  parts <- list(
    a = predictors_and_outcome(formula_a, rows),
    b = predictors_and_outcome(formula_b, rows)
  )
  if (!identical(parts$a$outcome, parts$b$outcome)) {
    stop("Both models need the same outcome.", call. = FALSE)
  }
  if (is.null(seed)) seed <- random_seed()

  ## each model with its own seed: Stan's, and that of R's random numbers for
  ## bridge sampling. (Only bridge sampling: cmdstanr names its output files
  ## with R's random numbers, and seeding those would give two fits with the
  ## same seed the same file names.)
  per_model <- lapply(c(a = 0, b = 1), function(offset) {
    model_parts <- parts[[offset + 1]]
    sampled <- sample_posterior(model_parts$predictors, model_parts$outcome, prior, chains,
                                iter_warmup, iter_sampling, adapt_delta, seed + offset, backend,
                                model)
    bridged <- log_marginal_likelihood(sampled, backend, density, method, seed + offset)
    model_details(sampled, bridged, backend, keep_fits)
  })
  log_bayes_factor <- per_model$a$log_marginal_likelihood - per_model$b$log_marginal_likelihood
  structure(list(bayes_factor     = exp(log_bayes_factor),
                 log_bayes_factor = log_bayes_factor,
                 error            = sqrt(per_model$a$error^2 + per_model$b$error^2),
                 log_marginal_likelihoods = c(a = per_model$a$log_marginal_likelihood,
                                              b = per_model$b$log_marginal_likelihood),
                 divergent        = c(a = per_model$a$divergent, b = per_model$b$divergent),
                 prior = prior, method = method, draws = chains * iter_sampling,
                 predictors = c(a = ncol(parts$a$predictors), b = ncol(parts$b$predictors)),
                 seed = seed, backend = backend, models = per_model),
            class = "flatcor_bayes_factor")
}

#' @export
print.flatcor_bayes_factor <- function(x, digits = 4, ...) {
  cat(sprintf("Bayes factor, model A (%d predictors) against model B (%d), %s prior\n",
              x$predictors[["a"]], x$predictors[["b"]], x$prior))
  cat(sprintf("  BF = %s   log BF = %s", format(x$bayes_factor, digits = digits),
              format(x$log_bayes_factor, digits = digits)))
  if (!is.na(x$error)) cat(sprintf("   (approximate SE of the log: %.4f)", x$error))
  cat(sprintf("\n  bridge sampling (%s) on %s posterior draws per model;",
              x$method, format(x$draws, big.mark = ",")))
  cat(sprintf(" divergent transitions: %d and %d\n", x$divergent[["a"]], x$divergent[["b"]]))
  invisible(x)
}
