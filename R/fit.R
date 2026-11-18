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
      executable_folder <- if (model_methods) {
        file.path(tempdir(), "flatcor_stan")
      } else {
        version <- as.character(utils::packageVersion("flatcor"))
        file.path(tools::R_user_dir("flatcor", which = "cache"), version)
      }
      dir.create(executable_folder, recursive = TRUE, showWarnings = FALSE)
      cmdstanr::cmdstan_model(
        stan_file,
        dir = executable_folder,
        compile_model_methods = model_methods,
        quiet = TRUE
      )
    } else {
      rstan::stan_model(stan_file)
    }
  }
  compiled_models[[key]]
}

## A random seed when none is given, drawn from R's random numbers (so
## set.seed() beforehand makes the result reproducible, as in rstanarm). One
## below the largest integer, so that bayes_factor_models() can give its
## second model seed + 1.
random_seed <- function() sample.int(.Machine$integer.max - 1, 1)

## Evaluate `code` with R's random numbers seeded by `seed`, and put the
## caller's random state back afterwards. For bridge sampling, which uses
## R's random numbers instead of Stan's.
with_seed <- function(seed, code) {
  caller_seed <- get0(".Random.seed", envir = globalenv())
  on.exit({
    if (is.null(caller_seed)) {
      rm(".Random.seed", envir = globalenv())
    } else {
      assign(".Random.seed", caller_seed, envir = globalenv())
    }
  })
  set.seed(seed)
  code
}

## cmdstanr when CmdStan itself is installed, otherwise rstan: having only the
## cmdstanr R package, without CmdStan, is not enough to fit anything
default_backend <- function() {
  has_cmdstan <- requireNamespace("cmdstanr", quietly = TRUE) &&
    !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))
  if (has_cmdstan) "cmdstanr" else "rstan"
}

## Sample the posterior of one model: predictors (a matrix, possibly with no
## columns) and outcome. Returns the backend's fit and its number of
## divergent transitions. per_observation: also draw log_lik and y_rep (for
## loo and pp_check; a Bayes factor does not need them).
sample_posterior <- function(predictors, outcome, prior, chains, iter_warmup, iter_sampling,
                             adapt_delta, seed, backend, model_methods = FALSE,
                             per_observation = FALSE) {
  check_full_rank(predictors)
  stan_data <- list(N = length(outcome), K = ncol(predictors), X = predictors, y = outcome,
                    prior_only = 0L, proper_stand_ins = 0L,
                    per_observation = as.integer(per_observation),
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

## The rows of data with no missing value in any variable of the formulas
## (all columns for a formula with "."). The sample size n enters every
## Bayes factor, so a dropped row is reported.
complete_rows <- function(data, formulas) {
  variables <- unique(unlist(lapply(formulas, all.vars)))
  if ("." %in% variables) variables <- names(data)
  unknown <- setdiff(variables, names(data))
  if (length(unknown) > 0)
    stop("Not in the data: ", paste(unknown, collapse = ", "), ".", call. = FALSE)
  complete <- stats::complete.cases(data[, variables, drop = FALSE])
  if (!all(complete)) {
    warning(sprintf(paste("%d of %d rows have a missing value in a variable of the formula",
                          "and were left out; the analysis uses n = %d."),
                    sum(!complete), nrow(data), sum(complete)), call. = FALSE)
  }
  data[complete, , drop = FALSE]
}

## The predictor columns (without the intercept) and the outcome of a formula
predictors_and_outcome <- function(formula, rows) {
  model_frame <- model.frame(formula, rows)
  design <- model.matrix(formula, model_frame)
  if (!"(Intercept)" %in% colnames(design)) {
    stop("The model needs an intercept; ", deparse(formula), " has none.", call. = FALSE)
  }
  list(
    predictors = design[, colnames(design) != "(Intercept)", drop = FALSE],
    outcome    = as.numeric(model.response(model_frame))
  )
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
#' @param backend \code{"cmdstanr"} (the default when CmdStan is installed)
#'   or \code{"rstan"} (the default otherwise).
#'   Use one backend per R session: once rstan is loaded, cmdstanr can no
#'   longer compute Bayes factors in that session. The validation
#'   (flatcor-validation) used cmdstanr.
#' @param prior \code{"flat_partial_r"} (every correlation matrix equally
#'   likely, R-squared ~ Beta(k/2, 1)), \code{"flat_r"} (R-squared ~
#'   Beta(1/2, 1) for every k) or \code{"flat_r_squared"} (R-squared ~
#'   Beta(1, 1) for every k). The default is \code{"flat_partial_r"}.
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
  parts <- predictors_and_outcome(formula, complete_rows(data, list(formula)))
  predictors <- parts$predictors
  outcome <- parts$outcome
  if (ncol(predictors) == 0) stop("fit_regression() needs at least one predictor.", call. = FALSE)
  if (is.null(seed)) seed <- random_seed()

  sampled <- sample_posterior(predictors, outcome, prior, chains, iter_warmup, iter_sampling,
                              adapt_delta, seed, backend, per_observation = TRUE)
  all_draws <- if (backend == "cmdstanr") {
    posterior::as_draws_df(sampled$fit$draws())
  } else {
    posterior::as_draws_df(as.array(sampled$fit))
  }
  draws_of <- function(variable) {
    unclass(posterior::as_draws_matrix(posterior::subset_draws(all_draws, variable = variable)))
  }

  ## Coefficients under their own names; the intercept back on the data's scale
  ## (the model's alpha is the mean of y, at the predictors' means)
  coefficients <- draws_of("b")
  colnames(coefficients) <- colnames(predictors)
  intercept <- all_draws$alpha - drop(coefficients %*% colMeans(predictors))
  draws <- posterior::as_draws_df(cbind(
    data.frame(Intercept = intercept, check.names = FALSE),
    as.data.frame(coefficients, check.names = FALSE),
    data.frame(sigma = all_draws$sigma, R_squared = all_draws$R_squared)
  ))

  structure(list(draws        = draws,
                 log_lik      = draws_of("log_lik"),
                 y_rep        = draws_of("y_rep"),
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
    rhat = posterior::rhat, ess_bulk = posterior::ess_bulk
  )
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
#' @return A \code{loo} object (see \code{loo::loo()}): the expected log
#'   predictive density, its SE, and the Pareto k diagnostics per observation.
#' @examples
#' \dontrun{
#' fit <- fit_regression(mpg ~ wt + hp, mtcars, prior = "flat_partial_r", seed = 1)
#' loo(fit)
#' }
#' @exportS3Method loo::loo
loo.flatcor_fit <- function(x, cores = 1, ...) {
  relative_efficiency <- loo::relative_eff(exp(x$log_lik), chain_id = x$chain_ids, cores = cores)
  loo::loo(x$log_lik, r_eff = relative_efficiency, cores = cores, ...)
}
#' Posterior predictive check for a flatcor fit
#' @param object A \code{flatcor_fit}.
#' @param ndraws Number of predicted datasets to draw.
#' @param ... Passed to \code{bayesplot::ppc_dens_overlay()}.
#' @return A ggplot: the density of the outcome against the densities of
#'   \code{ndraws} datasets predicted from the posterior.
#' @examples
#' \dontrun{
#' fit <- fit_regression(mpg ~ wt + hp, mtcars, prior = "flat_partial_r", seed = 1)
#' pp_check(fit)
#' }
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
