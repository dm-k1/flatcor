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

## Pull predictors and outcome out of an lm() fit or a flatcor_fit
extract_predictors_and_outcome <- function(model) {
  if (inherits(model, "flatcor_fit"))
    return(list(predictors = model$predictors, outcome = model$outcome))
  ## without an intercept the first column below would be a predictor
  if (attr(terms(model), "intercept") != 1)
    stop("Both models need an intercept; ", deparse(formula(model)), " has none.",
         call. = FALSE)
  if (!is.null(model$na.action))
    warning(sprintf("%d rows with a missing value were left out of %s; it uses n = %d.",
                    length(model$na.action), deparse(formula(model)), nobs(model)),
            call. = FALSE)
  list(predictors = model.matrix(model)[, -1, drop = FALSE],
       outcome    = as.numeric(model.response(model.frame(model))))
}

#' Bayes factor between any two linear models
#'
#' Nested or not, as long as both use the same outcome and data rows. Each
#' model is fitted with the chosen prior (\code{inst/stan/flatcor.stan}), its
#' marginal likelihood is estimated by bridge sampling
#' (\code{bridgesampling::bridge_sampler()}), and the Bayes factor is their
#' ratio. A model with no predictors (\code{lm(y ~ 1)}) is the
#' intercept-only model.
#'
#' The estimate has a Monte Carlo error. The reported error is
#' bridgesampling's approximate error of the log Bayes factor (method
#' \code{"normal"} only); it can be optimistic. To assess the error directly,
#' run the function again with other seeds and look at the spread, which the
#' bridgesampling authors recommend where it is affordable. Their largest
#' worked example has 96 parameters (Gronau et al., 2017); a model here has
#' k + 3 (k predictor columns: the intercept, the noise, R and a direction of
#' k numbers), so up to about 90 predictors stays within what has been shown.
#' @param model_a,model_b Fits from \code{lm()} or \code{fit_regression()}.
#' @param prior \code{"flat_partial_r"} (every correlation matrix equally
#'   likely, R-squared ~ Beta(k/2, 1)), \code{"flat_r"} (R-squared ~
#'   Beta(1/2, 1) for every k) or \code{"flat_r_squared"} (R-squared ~
#'   Beta(1, 1) for every k). The default is \code{"flat_partial_r"}.
#' @param chains,iter_warmup,iter_sampling,adapt_delta Sampler settings for
#'   each model (the default draws are ten times those for estimation).
#' @param seed A seed for the first model; the second uses \code{seed + 1}.
#' @param method \code{"normal"} (the default) or \code{"warp3"}, passed to
#'   \code{bridgesampling::bridge_sampler()}.
#' @param backend \code{"cmdstanr"} (default when installed) or \code{"rstan"}.
#' @return A \code{flatcor_bayes_factor}: \code{$bayes_factor} (model A
#'   against model B), \code{$log_bayes_factor}, \code{$error} (approximate
#'   SE of the log Bayes factor; \code{NA} for \code{"warp3"}), the two log
#'   marginal likelihoods and the divergent transitions of each fit.
#' @examples
#' \dontrun{
#' bayes_factor_models(lm(mpg ~ wt + hp, mtcars), lm(mpg ~ wt, mtcars),
#'                     prior = "flat_partial_r", seed = 1)
#' bayes_factor_models(lm(mpg ~ wt + hp, mtcars), lm(mpg ~ 1, mtcars),
#'                     prior = "flat_r", seed = 1)
#' }
#' @export
bayes_factor_models <- function(model_a, model_b, prior = "flat_partial_r", chains = 4,
                                iter_warmup = 1000, iter_sampling = 10000, adapt_delta = 0.95,
                                seed = NULL, method = "normal", backend = default_backend()) {
  prior  <- check_prior(prior, defaulted = missing(prior))
  method <- match.arg(method, c("normal", "warp3"))
  parts <- list(a = extract_predictors_and_outcome(model_a), b = extract_predictors_and_outcome(model_b))
  if (!isTRUE(all.equal(parts$a$outcome, parts$b$outcome)))
    stop("Both models need the same outcome and the same data rows.", call. = FALSE)
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max - 1, 1)

  per_model <- lapply(c(a = 0, b = 1), function(offset) {
    model_parts <- parts[[offset + 1]]
    sampled <- sample_posterior(model_parts$predictors, model_parts$outcome, prior, chains,
                                iter_warmup, iter_sampling, adapt_delta, seed + offset, backend,
                                model_methods = TRUE)
    bridge <- bridgesampling::bridge_sampler(sampled$fit, method = method, silent = TRUE)
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
  cat(sprintf("\n  bridge sampling (%s) on %s posterior draws per model; divergent transitions: %d and %d\n",
              x$method, format(x$draws, big.mark = ","), x$divergent[["a"]], x$divergent[["b"]]))
  invisible(x)
}
