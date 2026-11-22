## =====================================================================
##  Sampling the posterior: estimation (fit_regression()) and the fits
##  behind every Bayes factor (bayes_factor_models())
## =====================================================================
## The Stan model samples the prior in its own coordinates: R and a
## direction (uniform), with the coefficients computed from them. Sampling
## the coefficients directly misses the prior's peak at zero for larger k
## (flatcor-simulation-studies, studies 1a and 2a). With no predictors it is the
## intercept-only model.

## The model itself is built by flatcor_model() and its log density by
## flatcor_density() (R/model.R). compile_stan_model() is kept for code that
## called it before (flatcor-simulation-studies): the model for a backend. The model
## methods it used to add are no longer used; the Bayes factors take the log
## density from flatcor_density().
compile_stan_model <- function(backend, model_methods = FALSE) flatcor_model(backend)

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
## columns) and outcome, with `model` (flatcor_model()). Returns the backend's
## fit, its number of divergent transitions, the data Stan got, the seed and
## the seconds sampling took (wall clock, all chains).
## per_observation: also draw log_lik and y_rep (for loo and pp_check; a
## Bayes factor does not need them).
sample_posterior <- function(predictors, outcome, prior, chains, iter_warmup, iter_sampling,
                             adapt_delta, seed, backend, model = flatcor_model(backend),
                             per_observation = FALSE) {
  check_full_rank(predictors)
  check_model(model, backend)
  stan_data <- stan_data_of(predictors, outcome, prior, per_observation)
  started <- proc.time()[["elapsed"]]
  if (backend == "cmdstanr") {
    ## Each fit writes its CSV files to a folder of its own. cmdstanr names
    ## the files from the minute and R's random numbers, so two fits after the
    ## same set.seed() in one minute would share files, and cmdstanr deletes a
    ## fit's files once that fit is garbage: the other fit's draws would be
    ## gone. tempfile() does not use R's random numbers, so seeds and results
    ## are unchanged.
    ## Once the fit is read into memory (fit_in_memory()) nothing reads its
    ## files again, so the folder is deleted when this function ends, also
    ## when sampling fails. cmdstanr leaves files in a folder it was given,
    ## and a long-running R session (a validation worker) would otherwise
    ## fill the disk with them.
    output_folder <- tempfile("flatcor-fit-")
    dir.create(output_folder)
    on.exit(unlink(output_folder, recursive = TRUE), add = TRUE)
    fit <- model$sample(data = stan_data, chains = chains, parallel_chains = chains,
                        iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                        adapt_delta = adapt_delta, seed = seed, refresh = 0,
                        show_messages = FALSE, output_dir = output_folder)
    ## cmdstanr also wrote the data to a file of its own in R's temporary
    ## folder (standata-*.json); flatcor keeps the data itself ($stan_data)
    on.exit(unlink(fit$data_file()), add = TRUE)
    fit <- fit_in_memory(fit, backend)
    divergent <- as.integer(sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent))
  } else {
    fit <- rstan::sampling(model, data = stan_data, chains = chains, cores = chains,
                           iter = iter_warmup + iter_sampling, warmup = iter_warmup,
                           control = list(adapt_delta = adapt_delta), seed = seed, refresh = 0)
    divergent <- as.integer(rstan::get_num_divergent(fit))
  }
  seconds <- proc.time()[["elapsed"]] - started
  if (divergent > 0) warning(divergent, " divergent transitions; check the fit.", call. = FALSE)
  list(fit = fit, divergent = divergent, stan_data = stan_data, seed = seed, seconds = seconds)
}

## ---- What a fit holds besides its draws ---------------------------------------------
## Everything the backend computed comes out (the package's rule: nothing
## that would be there without flatcor goes missing because of it).

## All draws of every Stan variable, as a posterior::draws_array
stan_draws <- function(fit, backend) {
  if (backend == "cmdstanr") fit$draws() else posterior::as_draws_array(as.array(fit))
}

## The sampler's diagnostics per chain: divergent transitions, transitions
## that hit the maximum tree depth, E-BFMI, the adapted step size, the mean
## number of leapfrog steps, and the seconds of warmup and sampling Stan
## reports
sampler_by_chain <- function(fit, backend) {
  if (backend == "cmdstanr") {
    summary <- fit$diagnostic_summary(quiet = TRUE)
    per_draw <- posterior::as_draws_df(fit$sampler_diagnostics())
    times <- fit$time()$chains
    by_chain <- split(as.data.frame(per_draw), per_draw$.chain)
    data.frame(
      chain = seq_along(by_chain),
      divergent = as.integer(summary$num_divergent),
      max_treedepth_hits = as.integer(summary$num_max_treedepth),
      max_treedepth = fit$metadata()$max_treedepth,
      ebfmi = summary$ebfmi,
      step_size = vapply(by_chain, function(chain) chain$stepsize__[1], numeric(1)),
      mean_leapfrog_steps = vapply(by_chain, function(chain) mean(chain$n_leapfrog__), numeric(1)),
      warmup_seconds = times$warmup,
      sampling_seconds = times$sampling,
      row.names = NULL
    )
  } else {
    per_chain <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    max_treedepth <- fit@stan_args[[1]]$control$max_treedepth
    if (is.null(max_treedepth)) max_treedepth <- 10
    times <- rstan::get_elapsed_time(fit)
    data.frame(
      chain = seq_along(per_chain),
      divergent = vapply(per_chain, function(chain) {
        as.integer(sum(chain[, "divergent__"]))
      }, integer(1)),
      max_treedepth_hits = vapply(per_chain, function(chain) {
        sum(chain[, "treedepth__"] >= max_treedepth)
      }, integer(1)),
      max_treedepth = max_treedepth,
      ebfmi = unname(rstan::get_bfmi(fit)),
      step_size = vapply(per_chain, function(chain) chain[1, "stepsize__"], numeric(1)),
      mean_leapfrog_steps = vapply(per_chain, function(chain) mean(chain[, "n_leapfrog__"]),
                                   numeric(1)),
      warmup_seconds = unname(times[, "warmup"]),
      sampling_seconds = unname(times[, "sample"]),
      row.names = NULL
    )
  }
}

## A summary of every Stan variable: posterior's default summaries (mean,
## median, sd, mad, 5% and 95% quantiles), the 2.5% and 97.5% quantiles,
## R-hat and the bulk and tail ESS
summarise_stan_draws <- function(draws) {
  as.data.frame(posterior::summarise_draws(
    draws,
    "mean", "median", "sd", "mad",
    ~ posterior::quantile2(.x, probs = c(0.025, 0.05, 0.95, 0.975)),
    "rhat", "ess_bulk", "ess_tail"
  ))
}

## A cmdstanr fit reads its draws from CSV files, which sample_posterior()
## deletes right after this. Read everything into the object, so the fit
## keeps working, also when saved (what cmdstanr's $save_object() does).
## What needs the files themselves no longer works: $output_files() names
## files that are gone, and $cmdstan_diagnose() and $cmdstan_summary() fail.
fit_in_memory <- function(fit, backend) {
  if (backend == "cmdstanr") {
    fit$draws()
    try(fit$sampler_diagnostics(), silent = TRUE)
    try(fit$init(), silent = TRUE)
  }
  fit
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
#'   (flatcor-simulation-studies) used cmdstanr.
#' @param model The model to sample: \code{flatcor_model(backend)}, built once
#'   and reused (see \code{flatcor_model()} for parallel processes).
#' @param prior \code{"flat_partial_r"} (every correlation matrix equally
#'   likely, R-squared ~ Beta(k/2, 1)), \code{"flat_r"} (R-squared ~
#'   Beta(1/2, 1) for every k) or \code{"flat_r_squared"} (R-squared ~
#'   Beta(1, 1) for every k). The default is \code{"flat_partial_r"}.
#' @return A \code{flatcor_fit}: \code{print()} shows estimates and diagnostics,
#'   \code{loo()} and \code{pp_check()} work as usual. It holds:
#'   \describe{
#'     \item{\code{$draws}}{The posterior draws of the coefficients (under
#'       their own names, the intercept on the data's scale), sigma and
#'       R_squared, as a \code{posterior::draws_df}.}
#'     \item{\code{$log_lik}, \code{$y_rep}}{Per draw and observation, for
#'       \code{loo()} and \code{pp_check()}.}
#'     \item{\code{$fit}}{The backend's own fit (\code{CmdStanMCMC} or
#'       \code{stanfit}) with every Stan variable, such as R, Cohen's f and
#'       the direction; its draws are read into memory, so a saved fit keeps
#'       working.}
#'     \item{\code{$sampler}}{The sampler's diagnostics per chain: divergent
#'       transitions, hits of the maximum tree depth, E-BFMI, step size, mean
#'       leapfrog steps, seconds of warmup and sampling.}
#'     \item{\code{$seconds}, \code{$seed}, \code{$backend},
#'       \code{$stan_data}}{Wall-clock seconds of sampling, the seed used
#'       (drawn when none was given), the backend and the data Stan got.}
#'     \item{\code{$outcome}, \code{$predictors}, \code{$formula},
#'       \code{$prior}, \code{$divergent}, \code{$chain_ids}}{What was
#'       fitted, and the divergent transitions in all chains.}
#'   }
#' @examples
#' # needs CmdStan (cmdstanr::install_cmdstan())
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
#'   fit <- fit_regression(mpg ~ wt + hp, mtcars, prior = "flat_partial_r")
#'   fit
#'   loo(fit)
#' }
#' @export
fit_regression <- function(formula, data, prior = "flat_partial_r", chains = 4, iter_warmup = 1000,
                           iter_sampling = 1000, adapt_delta = 0.95, seed = NULL,
                           backend = default_backend(), model = flatcor_model(backend)) {
  prior <- check_prior(prior, defaulted = missing(prior))
  parts <- predictors_and_outcome(formula, complete_rows(data, list(formula)))
  predictors <- parts$predictors
  outcome <- parts$outcome
  if (ncol(predictors) == 0) stop("fit_regression() needs at least one predictor.", call. = FALSE)
  if (is.null(seed)) seed <- random_seed()

  sampled <- sample_posterior(predictors, outcome, prior, chains, iter_warmup, iter_sampling,
                              adapt_delta, seed, backend, model, per_observation = TRUE)
  all_draws <- posterior::as_draws_df(stan_draws(sampled$fit, backend))
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
                 divergent    = sampled$divergent,
                 sampler      = sampler_by_chain(sampled$fit, backend),
                 seconds      = sampled$seconds,
                 seed         = seed,
                 backend      = backend,
                 stan_data    = sampled$stan_data,
                 fit          = sampled$fit),
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
#' # needs CmdStan (cmdstanr::install_cmdstan())
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
#'   fit <- fit_regression(mpg ~ wt + hp, mtcars, prior = "flat_partial_r", seed = 1)
#'   loo(fit)
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
#' # needs CmdStan (cmdstanr::install_cmdstan())
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
#'   fit <- fit_regression(mpg ~ wt + hp, mtcars, prior = "flat_partial_r", seed = 1)
#'   pp_check(fit)
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
