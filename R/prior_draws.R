## =====================================================================
##  The prior as a recipe
## =====================================================================

#' Draw coefficient sets from one of the three priors
#'
#' The prior as a recipe, in units of the noise (sigma = 1):
#' 1. P(total): draw the overall correlation R; Cohen's f^2 = R^2 / (1 - R^2).
#'    R^2 ~ Beta(a, 1), i.e. R = U^(1/(2a)) with U flat: a = k/2 for
#'    \code{"flat_partial_r"}, 1/2 for \code{"flat_r"} (R flat) and 1 for
#'    \code{"flat_r_squared"} (R^2 flat).
#' 2. P(combination | total): draw a direction with no preference in explained
#'    variance (normal draws with covariance solve(Sigma)).
#' 3. Scale each direction so its f^2 is exactly the drawn amount.
#' @param number_of_draws Number of coefficient sets.
#' @param predictor_covariance Covariance matrix of the predictors.
#' @param prior \code{"flat_partial_r"} (the default), \code{"flat_r"} or
#'   \code{"flat_r_squared"}.
#' @return A matrix with one coefficient set per row, with what each set was
#'   drawn from as attributes: \code{"R_squared"} and \code{"f_squared"}
#'   (one per row) and \code{"direction"} (the unscaled direction, a matrix
#'   like the result).
#' @examples
#' draws <- draw_coefficients_from_prior(1000, diag(3), prior = "flat_partial_r")
#' r_squared <- rowSums(draws^2) / (1 + rowSums(draws^2))   # ~ Beta(3/2, 1)
#' @export
draw_coefficients_from_prior <- function(number_of_draws, predictor_covariance,
                                         prior = "flat_partial_r") {
  prior <- check_prior(prior, defaulted = missing(prior))
  number_of_predictors <- ncol(predictor_covariance)
  shape <- prior_shape(prior, number_of_predictors)
  overall_correlation <- runif(number_of_draws, min = 0, max = 1)^(1 / (2 * shape))
  f_squared  <- overall_correlation^2 / (1 - overall_correlation^2)
  directions <- MASS::mvrnorm(number_of_draws, mu = rep(0, number_of_predictors),
                              Sigma = solve(predictor_covariance))
  # keep a matrix when k = 1
  directions <- matrix(directions, number_of_draws, number_of_predictors)
  f_squared_of_direction <- rowSums((directions %*% predictor_covariance) * directions)
  coefficients <- sqrt(f_squared / f_squared_of_direction) * directions
  structure(coefficients, R_squared = overall_correlation^2, f_squared = f_squared,
            direction = directions)
}
