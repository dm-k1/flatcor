#' flatcor: objective Bayes for linear regression with flat correlations
#'
#' Define what's possible, not what you expect. With one predictor, Jeffreys'
#' flat correlation is the only choice; with several, not everything can stay
#' flat. Three priors, built and validated side by side, all R-squared ~
#' Beta(a, 1) with a uniform direction:
#' \itemize{
#'   \item \code{"flat_partial_r"}: every correlation matrix of the outcome and
#'     the predictors equally likely; each partial correlation flat on (-1, 1)
#'     whatever the number of predictors (a = k/2).
#'   \item \code{"flat_r"}: the multiple correlation R flat for every k (a = 1/2).
#'   \item \code{"flat_r_squared"}: R-squared flat for every k (a = 1).
#' }
#' The default is \code{"flat_partial_r"}; a function that uses it because no
#' prior was given says so once per session. The validation of the three
#' priors is at \url{https://github.com/dm-k1/flatcor-simulation-studies}. All three divide
#' the total over the predictors with no preference:
#' P(total) x P(this combination | total) x unit conversion.
#'
#' One Stan model (\code{inst/stan/flatcor.stan}) does everything, for all
#' three priors: the posterior (\code{fit_regression()}) and, through bridge
#' sampling on it, the Bayes factor between two models
#' (\code{bayes_factor_models()}). \code{draw_coefficients_from_prior()}
#' gives the prior as a recipe.
#'
#' @importFrom stats model.frame model.matrix model.response terms formula nobs
#'   runif quantile sd
#' @keywords internal
"_PACKAGE"
