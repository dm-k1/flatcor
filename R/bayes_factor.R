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
#' @return A \code{flatcor_bayes_factor}: \code{$bayes_factor} (model A
#'   against model B), \code{$log_bayes_factor}, \code{$error} (approximate
#'   SE of the log Bayes factor; \code{NA} for \code{"warp3"}), the two log
#'   marginal likelihoods and the divergent transitions of each fit.
#' @examples
#' \dontrun{
#' bayes_factor_models(mpg ~ wt + hp, mpg ~ wt, mtcars, prior = "flat_partial_r", seed = 1)
#' bayes_factor_models(mpg ~ wt + hp, mpg ~ 1, mtcars, prior = "flat_r", seed = 1)
#' }
#' @export
bayes_factor_models <- function(formula_a, formula_b, data, prior = "flat_partial_r", chains = 4,
                                iter_warmup = 1000, iter_sampling = 10000, adapt_delta = 0.95,
                                seed = NULL, method = "normal", backend = default_backend()) {
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
                                model_methods = TRUE)
    bridge <- with_seed(
      seed + offset,
      bridgesampling::bridge_sampler(sampled$fit, method = method, silent = TRUE)
    )
    list(log_marginal_likelihood = bridge$logml,
         error = if (method == "normal") bridgesampling::error_measures(bridge)$cv else NA_real_,
         divergent = sampled$divergent)
  })
  log_bayes_factor <- per_model$a$log_marginal_likelihood - per_model$b$log_marginal_likelihood
  structure(list(bayes_factor     = exp(log_bayes_factor),
                 log_bayes_factor = log_bayes_factor,
                 error            = sqrt(per_model$a$error^2 + per_model$b$error^2),
                 log_marginal_likelihoods = c(a = per_model$a$log_marginal_likelihood,
                                              b = per_model$b$log_marginal_likelihood),
                 divergent        = c(a = per_model$a$divergent, b = per_model$b$divergent),
                 prior = prior, method = method, draws = chains * iter_sampling,
                 predictors = c(a = ncol(parts$a$predictors), b = ncol(parts$b$predictors))),
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
