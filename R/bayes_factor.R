## =====================================================================
##  Bayes factors for whole models, for the three priors
## =====================================================================
## The three priors are R^2 ~ Beta(a, 1) with a uniform direction (R/priors.R):
##   flat_partial_r   a = k/2: every correlation matrix of (y, x1 ... xk)
##                    equally likely (stan_lm with eta = 1)
##   flat_r           a = 1/2: R flat for every k
##   flat_r_squared   a = 1:   R^2 flat for every k
## Every model gets an exact Bayes factor against the intercept-only model; any
## two models then compare by dividing. The data enter only through the sample
## R-squared, n and k.
##
## One route for every prior: the likelihood ratio of the model's k columns
## against none (R/likelihood.R), averaged over the prior on the length f of
## the effect in noise units by integrate(). R^2 = f^2 / (1 + f^2), so
## R^2 ~ Beta(a, 1) fixes the prior on f, and only a differs between the priors.
## (Closed forms through g-prior mixtures exist for each prior; they are kept
## in the tests, tests/testthat/helper-closed-forms.R, as references.)

log_sum_exp <- function(log_values) {
  largest <- max(log_values)
  largest + log(sum(exp(log_values - largest)))
}

## Sample R-squared (always with intercept). Pure arithmetic, no fitting.
## A column that is a combination of the others (aliased) adds nothing to
## R-squared but would still count in k, so it stops instead.
sample_r_squared <- function(predictor_columns, outcome) {
  if (is.null(predictor_columns) || NCOL(predictor_columns) == 0) return(0)
  decomposition <- qr(cbind(1, predictor_columns))
  if (decomposition$rank < ncol(decomposition$qr))
    stop("The predictor columns are linearly dependent (one is a combination of ",
         "the others, or constant). Remove the redundant columns.")
  residuals_left <- qr.resid(decomposition, outcome)
  1 - sum(residuals_left^2) / sum((outcome - mean(outcome))^2)
}

## Any a: the likelihood ratio of the k columns against none, averaged over the
## prior on f, the length of the effect in noise units. R^2 = f^2 / (1 + f^2),
## so R^2 ~ Beta(a, 1) gives |f| the density 2a f^(2a-1) / (1 + f^2)^(a+1).
## The average is one call to integrate() (adaptive quadrature) over log f,
## centred on the peak of the integrand, as BayesFactor computes its regression
## Bayes factor by integrate() over log g. Only a differs between the priors.
## At n = k + 1 every model fits the data exactly (sample R^2 = 1) and the data
## say nothing about the effect: the Bayes factor is 1 for any proper prior on
## it (predictive matching, Bayarri et al., 2012).
log_bayes_factor_shape_vs_intercept_only <- function(explained_share, sample_size,
                                                     number_of_predictors, shape) {
  if (number_of_predictors == 0) return(0)
  if (sample_size < number_of_predictors + 1)
    stop("Need at least one more observation than predictors.")
  if (sample_size == number_of_predictors + 1) return(0)
  r_squared <- pmin(explained_share, 1 - 1e-12)
  ## log of likelihood ratio x prior density of f x f, as a function of log f
  ## (beyond f = 1e8 or below 1e-8 the integrand is negligible)
  log_integrand <- function(log_f) {
    f <- exp(log_f)
    value <- log_likelihood_ratio_block(f, r_squared, sample_size, number_of_predictors,
                                        number_of_predictors) +
      log(2 * shape) + 2 * shape * log_f - (shape + 1) * log1p(f^2)
    ifelse(is.finite(value), value, -Inf)
  }
  ## the peak: on a coarse grid, then refined; its width from the curvature
  limits <- log(c(1e-8, 1e8))
  grid   <- seq(limits[1], limits[2], by = 1)
  start  <- grid[which.max(log_integrand(grid))]
  peak   <- stats::optimize(log_integrand, start + c(-1, 1), maximum = TRUE,
                            tol = 1e-10)$maximum
  top    <- log_integrand(peak)
  step   <- 1e-4
  curvature <- (log_integrand(peak + step) - 2 * top + log_integrand(peak - step)) / step^2
  width  <- if (is.finite(curvature) && curvature < 0) 1 / sqrt(-curvature) else 1
  breaks <- sort(unique(pmin(limits[2], pmax(limits[1],
    c(limits, peak + width * c(-30, -10, -4, -1, 0, 1, 4, 10, 30))))))
  pieces <- vapply(seq_len(length(breaks) - 1), function(i)
    stats::integrate(function(x) exp(log_integrand(x) - top), breaks[i], breaks[i + 1],
                     rel.tol = 1e-10, subdivisions = 1000L)$value, numeric(1))
  top + log(sum(pieces))
}

#' Log Bayes factor of a model with k predictors against the intercept-only model
#'
#' The data enter only through the sample R-squared, n and k. The three priors
#' are R-squared ~ Beta(a, 1) with a uniform direction: \code{"flat_partial_r"}
#' (a = k/2, every correlation matrix of the outcome and predictors equally
#' likely), \code{"flat_r"} (a = 1/2, R flat) and \code{"flat_r_squared"}
#' (a = 1, R-squared flat).
#' @param explained_share Sample R-squared of the model.
#' @param sample_size Number of observations, n.
#' @param number_of_predictors Number of predictor columns, k.
#' @param prior \code{"flat_partial_r"} (the default), \code{"flat_r"} or
#'   \code{"flat_r_squared"}.
#' @return The log Bayes factor.
#' @details One route for every prior: the likelihood ratio of the model's
#'   columns against none, with the noise and the direction of the effect
#'   integrated out, averaged over R-squared ~ Beta(a, 1) by
#'   \code{stats::integrate()}. Checked on 174 settings (n = 3 to 10^6, k = 1
#'   to 50, R-squared from 10^-6 to .999999): for flat partial r within 2e-7
#'   of adaptive integration of its g-prior form; for flat R and flat
#'   R-squared within the error of the closed forms used up to 1.0.2, which
#'   were the less accurate side (refining them, or a third route, gives this
#'   route's value). The validation checks it against independent routes
#'   (objective-bayes, study 2).
#' @examples
#' log_bayes_factor_vs_intercept_only(0.2, 50, 3, prior = "flat_partial_r")
#' log_bayes_factor_vs_intercept_only(0.2, 50, 3, prior = "flat_r_squared")
#' @export
log_bayes_factor_vs_intercept_only <- function(explained_share, sample_size, number_of_predictors,
                                               prior = "flat_partial_r") {
  prior <- check_prior(prior, defaulted = missing(prior))
  log_bayes_factor_shape_vs_intercept_only(explained_share, sample_size, number_of_predictors,
                                           shape = prior_shape(prior, number_of_predictors))
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

## Pull predictors and outcome out of a flatcor_fit or lm fit
extract_predictors_and_outcome <- function(model) {
  if (inherits(model, "flatcor_fit")) {
    predictors <- model$predictors
    outcome    <- model$outcome
  } else {
    ## without an intercept the first column below would be a predictor
    if (attr(terms(model), "intercept") != 1)
      stop("Both models need an intercept; ", deparse(formula(model)), " has none.",
           call. = FALSE)
    if (!is.null(model$na.action))
      warning(sprintf("%d rows with a missing value were left out of %s; it uses n = %d.",
                      length(model$na.action), deparse(formula(model)), nobs(model)),
              call. = FALSE)
    predictors <- model.matrix(model)[, -1, drop = FALSE]
    outcome    <- as.numeric(model.response(model.frame(model)))
  }
  list(predictors = predictors, outcome = outcome)
}

#' Bayes factor between any two linear models
#'
#' Nested or not, as long as both use the same outcome and data rows.
#' @param model_a,model_b Fits from \code{fit_regression()} or \code{lm()}.
#' @param prior \code{"flat_partial_r"} (every correlation matrix equally
#'   likely, R-squared ~ Beta(k/2, 1)), \code{"flat_r"} (R-squared ~
#'   Beta(1/2, 1) for every k) or \code{"flat_r_squared"} (R-squared ~
#'   Beta(1, 1) for every k). The default is \code{"flat_partial_r"}, the
#'   prior the validation recommends
#'   (\url{https://github.com/dm-k1/objective-bayes}).
#' @return The Bayes factor of model A against model B.
#' @examples
#' bayes_factor_models(lm(mpg ~ wt + hp, mtcars), lm(mpg ~ wt, mtcars),
#'                     prior = "flat_partial_r")
#' bayes_factor_models(lm(mpg ~ wt + hp, mtcars), lm(mpg ~ wt, mtcars),
#'                     prior = "flat_r")
#' @export
bayes_factor_models <- function(model_a, model_b, prior = "flat_partial_r") {
  prior <- check_prior(prior, defaulted = missing(prior))
  parts_a <- extract_predictors_and_outcome(model_a)
  parts_b <- extract_predictors_and_outcome(model_b)
  if (!isTRUE(all.equal(parts_a$outcome, parts_b$outcome)))
    stop("Both models need the same outcome and the same data rows.")
  sample_size <- length(parts_a$outcome)
  log_bayes_factor_of <- function(parts)
    log_bayes_factor_vs_intercept_only(
      sample_r_squared(parts$predictors, parts$outcome), sample_size, ncol(parts$predictors), prior)
  exp(log_bayes_factor_of(parts_a) - log_bayes_factor_of(parts_b))
}
