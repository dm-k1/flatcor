## =====================================================================
##  The prior as a recipe
## =====================================================================

#' Draw coefficient sets from one of the four priors
#'
#' The prior as a recipe, in units of the noise (sigma = 1). The three
#' R-squared ~ Beta(a, 1) priors (\code{"joint_partial_r"}, \code{"flat_r"} and
#' \code{"flat_r_squared"}) share one recipe:
#' 1. P(total): draw the overall correlation R; Cohen's f^2 = R^2 / (1 - R^2).
#'    R^2 ~ Beta(a, 1), i.e. R = U^(1/(2a)) with U flat: a = k/2 for
#'    \code{"joint_partial_r"}, 1/2 for \code{"flat_r"} (R flat) and 1 for
#'    \code{"flat_r_squared"} (R^2 flat).
#' 2. P(combination | total): draw a direction with no preference in explained
#'    variance (normal draws with covariance solve(Sigma)).
#' 3. Scale each direction so its f^2 is exactly the drawn amount.
#'
#' \code{"independent_partial_r"} has its own recipe, which does not go
#' through the total:
#' 1. Draw each effect's partial correlation (with the outcome, given the
#'    other predictors) flat on its own: rho_j ~ U(-1, 1).
#' 2. What each predictor has that the others don't: its residual SD given
#'    the others, s_j = 1 / sqrt(diag(solve(Sigma))).
#' 3. Turn each rho_j into its coefficient: b_j = rho_j / (s_j sqrt(1 - rho_j^2))
#'    (sigma = 1).
#' 4. What this implies for the total: f^2 = b' Sigma b, R^2 = f^2 / (1 + f^2).
#'    For k > 1 this R^2 has no closed form and depends on Sigma; these
#'    draws are how to see it.
#' @param number_of_draws Number of coefficient sets.
#' @param predictor_covariance Covariance matrix of the predictors.
#' @param prior \code{"joint_partial_r"} (the default),
#'   \code{"independent_partial_r"}, \code{"flat_r"} or \code{"flat_r_squared"}.
#' @return A matrix with one coefficient set per row, with what each set was
#'   drawn from as attributes: \code{"R_squared"} and \code{"f_squared"}
#'   (one per row), \code{"direction"} (the unscaled direction, a matrix
#'   like the result; for \code{"independent_partial_r"} the coefficients
#'   themselves) and \code{"partial_r"} (each coefficient's partial
#'   correlation with the outcome, a matrix like the result).
#' @examples
#' draws <- draw_coefficients_from_prior(1000, diag(3), prior = "joint_partial_r")
#' r_squared <- rowSums(draws^2) / (1 + rowSums(draws^2))   # ~ Beta(3/2, 1)
#' @export
draw_coefficients_from_prior <- function(number_of_draws, predictor_covariance,
                                         prior = "joint_partial_r") {
  prior <- check_prior(prior, defaulted = missing(prior))
  draws <- if (prior_is_independent(prior)) {
    draw_from_independent_partial_r(number_of_draws, predictor_covariance)
  } else {
    draw_from_joint_family(number_of_draws, predictor_covariance, prior)
  }
  # each coefficient's partial correlation: rho_j = b_j s_j / sqrt(1 + (b_j s_j)^2)
  residual_sd <- 1 / sqrt(diag(solve(predictor_covariance)))
  scaled_effects <- sweep(draws, 2, residual_sd, "*")
  attr(draws, "partial_r") <- scaled_effects / sqrt(1 + scaled_effects^2)
  draws
}

## The recipe of the three R-squared ~ Beta(a, 1) priors (R-squared ~ Beta(a, 1), uniform
## direction), steps 1 to 3 above.
draw_from_joint_family <- function(number_of_draws, predictor_covariance, prior) {
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

## The recipe of independent_partial_r, steps 1 to 4 above.
draw_from_independent_partial_r <- function(number_of_draws, predictor_covariance) {
  number_of_predictors <- ncol(predictor_covariance)

  # step 1: each effect's partial correlation, flat on its own
  partial_correlations <- matrix(
    runif(number_of_draws * number_of_predictors, min = -1, max = 1),
    number_of_draws,
    number_of_predictors
  )

  # step 2: what each predictor has that the others don't
  residual_sd <- 1 / sqrt(diag(solve(predictor_covariance)))

  # step 3: each partial correlation to its coefficient (sigma = 1)
  coefficients <- sweep(
    partial_correlations / sqrt(1 - partial_correlations^2),
    2,
    residual_sd,
    "/"
  )

  # step 4: what this implies for the total
  f_squared <- rowSums((coefficients %*% predictor_covariance) * coefficients)
  structure(
    coefficients,
    R_squared = f_squared / (1 + f_squared),
    f_squared = f_squared,
    direction = coefficients
  )
}
