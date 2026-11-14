## Closed forms of the whole-model Bayes factor, one per prior, through
## g-prior mixtures (objective-bayes/docs/appendix/A_derivations.md, section 8).
## flatcor computes every prior by one route (R/bayes_factor.R); these were
## its formulas up to 1.0.2 and are kept here as independent references. They
## have grid errors of their own: up to about 3e-5 for flat R (its direction
## grid) and 5e-4 for flat R-squared at n = 10^6 (its table of the mixing
## density), so the tests compare with tolerances of 1e-4 to 2e-3.

log_sum_exp_per_row <- function(log_matrix) {
  largest_per_row <- apply(log_matrix, 1, max)
  largest_per_row + log(rowSums(exp(log_matrix - largest_per_row)))
}

## flat_r, closed form: the one-predictor Bayes factor averaged over every
## direction the effect could take.
log_bayes_factor_flat_r_vs_intercept_only <- function(explained_share, sample_size,
                                                      number_of_predictors,
                                                      grid_points_g = 800,
                                                      grid_points_direction = 200) {
  if (number_of_predictors == 0) return(0)
  if (sample_size < number_of_predictors + 1)
    stop("Need at least one more observation than predictors.")

  ## 1. Average over g: the one-predictor prior written as a mixture of
  ##    normal priors with variance g (g ~ Inverse-Gamma(1, n/2)), on a log-g grid.
  g_grid   <- inverse_gamma_g_grid(sample_size, grid_points_g)
  g_values <- g_grid$g_values

  ## 2. Average over directions: only the alignment matters, the cosine of the
  ##    angle between a direction and the direction the data point in
  if (number_of_predictors == 1) {
    alignment_squared     <- 1
    log_direction_weights <- 0
  } else {
    ## Angle between a direction and the data's direction. Only (0, pi/2) is
    ## needed: angles phi and pi - phi give the same alignment. The grid packs
    ## its points near angle 0, where the integrand becomes a narrow peak when
    ## the evidence is strong (packing_power = 1 would be an even grid).
    grid_position <- (seq_len(grid_points_direction) - 0.5) / grid_points_direction
    packing_power <- 3
    angles        <- (pi / 2) * grid_position^packing_power
    angle_steps   <- (pi / 2) * packing_power * grid_position^(packing_power - 1) /
                     grid_points_direction
    ## exact normalizing constant: integral of sin^(k-2) over (0, pi/2)
    log_total_weight <- 0.5 * log(pi) - log(2) +
      lgamma((number_of_predictors - 1) / 2) - lgamma(number_of_predictors / 2)
    alignment_squared     <- cos(angles)^2
    log_direction_weights <- (number_of_predictors - 2) * log(sin(angles)) +
                             log(angle_steps) - log_total_weight
  }

  shrink_factor  <- g_values / (1 + g_values)
  log_fit_matrix <- -(sample_size - 1) / 2 *
    log1p(-pmin(explained_share, 1 - 1e-12) * outer(shrink_factor, alignment_squared)) +
    matrix(log_direction_weights, grid_points_g, length(alignment_squared), byrow = TRUE)

  ## 3. The one-predictor Bayes factor for each (g, alignment), averaged
  log_sum_exp(g_grid$log_weights - 0.5 * log1p(g_values) +
              log_sum_exp_per_row(log_fit_matrix))
}

## g ~ Inverse-Gamma(1, n/2) on a log-g grid: the scale of both priors
## (objective-bayes/docs/appendix/A_derivations.md, section 8: g/n = share x W, W ~ IG(1, 1/2)).
## Each weight is the IG density times g times the step (the step is in log g).
inverse_gamma_g_grid <- function(sample_size, grid_points_g) {
  log_g_values <- seq(log(1e-6), log(1e10), length.out = grid_points_g)
  g_values     <- exp(log_g_values)
  list(g_values    = g_values,
       log_weights = log(sample_size / 2) - log_g_values - sample_size / (2 * g_values) +
                     log(log_g_values[2] - log_g_values[1]))
}

## flat_partial_r, closed form: R-squared ~ Beta(k/2, 1) with a uniform
## direction (rstanarm stan_lm with eta = 1). As a mixture of g-priors its
## scale is the same for every k, g ~ Inverse-Gamma(1, n/2), so the Bayes factor
## is the standard g-prior Bayes factor averaged over g: one integral, on a
## log-g grid.
log_bayes_factor_flat_partial_r_vs_intercept_only <- function(explained_share, sample_size,
                                                                    number_of_predictors,
                                                                    grid_points_g = 800) {
  if (number_of_predictors == 0) return(0)
  if (sample_size < number_of_predictors + 1)
    stop("Need at least one more observation than predictors.")
  g_grid   <- inverse_gamma_g_grid(sample_size, grid_points_g)
  g_values <- g_grid$g_values
  ## the g-prior Bayes factor (Liang et al., 2008) for each g, averaged
  log_g_prior_bayes_factor <- (sample_size - 1 - number_of_predictors) / 2 * log1p(g_values) -
    (sample_size - 1) / 2 * log1p(g_values * (1 - pmin(explained_share, 1 - 1e-12)))
  log_sum_exp(g_grid$log_weights + log_g_prior_bayes_factor)
}

## flat_r_squared for k >= 3, closed form: R^2 ~ Beta(1, 1) with a uniform
## direction is the g-prior mixture g/n = t = B W, B ~ Beta(1, k/2 - 1),
## W ~ Inverse-Gamma(1, 1/2) (objective-bayes/docs/appendix/A_derivations.md, section 8). The
## density of t, p(t) = E[B exp(-B / (2t))] / (2 t^2), depends only on k: it is
## computed once per k (integrate() over B at each point of a log-t grid) and
## kept, so each Bayes factor is one sum over the grid, as for flat_partial_r.
## (k = 2: B = 1, the flat_partial_r prior; k = 1: the block route.)
flat_r_squared_hyperprior <- new.env()
log_t_grid_flat_r_squared <- function(number_of_predictors, points_t = 1600) {
  key <- as.character(number_of_predictors)
  if (is.null(flat_r_squared_hyperprior[[key]])) {
    log_t <- seq(log(1e-9), log(1e9), length.out = points_t)
    other_shape <- number_of_predictors / 2 - 1
    log_density <- vapply(exp(log_t), function(t) {
      weighted_beta <- function(b) b * stats::dbeta(b, 1, other_shape) * exp(-b / (2 * t))
      split <- min(1, 60 * t)   # the integrand sits near b = 0 when t is small
      value <- stats::integrate(weighted_beta, 0, split, rel.tol = 1e-10)$value +
        (if (split < 1) stats::integrate(weighted_beta, split, 1, rel.tol = 1e-10)$value else 0)
      log(value) - log(2) - 2 * log(t)
    }, numeric(1))
    flat_r_squared_hyperprior[[key]] <- list(
      log_t = log_t, log_weights = log_density + log_t + log(log_t[2] - log_t[1]))
  }
  flat_r_squared_hyperprior[[key]]
}
log_bayes_factor_flat_r_squared_vs_intercept_only <- function(explained_share, sample_size,
                                                              number_of_predictors) {
  if (number_of_predictors == 0) return(0)
  if (number_of_predictors <= 2)
    return(if (number_of_predictors == 2)
      log_bayes_factor_flat_partial_r_vs_intercept_only(explained_share, sample_size, 2) else
      log_bayes_factor_shape_vs_intercept_only(explained_share, sample_size, 1, shape = 1))
  if (sample_size < number_of_predictors + 1)
    stop("Need at least one more observation than predictors.")
  grid <- log_t_grid_flat_r_squared(number_of_predictors)
  g_values <- sample_size * exp(grid$log_t)
  log_g_prior_bayes_factor <- (sample_size - 1 - number_of_predictors) / 2 * log1p(g_values) -
    (sample_size - 1) / 2 * log1p(g_values * (1 - pmin(explained_share, 1 - 1e-12)))
  keep <- is.finite(grid$log_weights)
  log_sum_exp(grid$log_weights[keep] + log_g_prior_bayes_factor[keep])
}

