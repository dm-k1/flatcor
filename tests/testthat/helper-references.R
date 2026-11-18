## Independent references for the whole-model Bayes factor against the
## intercept-only model: the prior as a mixture of g-priors
## (Objective-bayes/docs/appendix/A_derivations.md, sections 7 and 8),
## evaluated by integrate(). They share nothing with the package, which
## estimates each marginal likelihood from its Stan model by bridge sampling.
## A g-mixture exists for a <= k/2; for flat R the tests use the direction
## form at the end of this file, which is more precise than the B-mixture
## when Beta(1/2, (k - 1)/2) is sharp (large k).

## The sample R-squared of outcome y on the columns of predictors (with an
## intercept): 1 - residual sum of squares / total sum of squares, from R's
## QR decomposition. The references below take it as their input.
sample_r_squared <- function(y, predictors) {
  residuals <- qr.resid(qr(cbind(1, predictors)), y)
  1 - sum(residuals^2) / sum((y - mean(y))^2)
}

## the g-prior Bayes factor for k predictors, on the log scale
log_g_prior_bayes_factor <- function(g, r_squared, n, k) {
  (n - 1 - k) / 2 * log1p(g) - (n - 1) / 2 * log1p(g * (1 - r_squared))
}

## integrate exp(log_integrand) over x, centred on its peak in (lower, upper)
log_integral_around_peak <- function(log_integrand, lower, upper) {
  peak   <- optimize(log_integrand, c(lower, upper), maximum = TRUE)$maximum
  height <- log_integrand(peak)
  breaks <- sort(unique(pmin(upper, pmax(lower, peak + c(-Inf, -20, -5, -1, 0, 1, 5, 20, Inf)))))
  piece <- function(i) {
    integrate(
      function(x) exp(log_integrand(x) - height),
      breaks[i], breaks[i + 1],
      rel.tol = 1e-10, subdivisions = 2000L
    )$value
  }
  height + log(sum(vapply(seq_len(length(breaks) - 1), piece, numeric(1))))
}

## R-squared ~ Beta(a, 1) as a g-mixture: g/n = t = B W, W ~ Inverse-Gamma(1, 1/2),
## B ~ Beta(a, k/2 - a), or B = 1 when a = k/2 (flat partial r); so
## p(t) = E[B exp(-B / (2t))] / (2 t^2)
log_density_g_over_n <- function(t, k, a) {
  if (isTRUE(all.equal(a, k / 2))) return(-log(2) - 2 * log(t) - 1 / (2 * t))
  ## over u = b^a, which removes the singularity of Beta(a, .) at 0 (a < 1)
  weighted <- function(u) {
    b <- u^(1 / a)
    b * exp(-b / (2 * t)) * (1 - b)^(k / 2 - a - 1) / (a * beta(a, k / 2 - a))
  }
  split <- min(1, (60 * t)^a)
  below <- integrate(weighted, 0, split, rel.tol = 1e-10)$value
  above <- if (split < 1) integrate(weighted, split, 1, rel.tol = 1e-10)$value else 0
  log(below + above) - log(2) - 2 * log(t)
}
log_bayes_factor_by_g_mixture <- function(r_squared, n, k, a) {
  stopifnot(a <= k / 2)
  log_integrand <- function(log_t) {
    vapply(log_t, function(lt) {
      value <- log_g_prior_bayes_factor(n * exp(lt), r_squared, n, k) +
        log_density_g_over_n(exp(lt), k, a) + lt
      if (is.finite(value)) value else -Inf
    }, numeric(1))
  }
  log_integral_around_peak(log_integrand, -25, 25)
}

## flat R in the form of appendix A, section 7, steps 1 to 3: the one-predictor
## g-prior Bayes factor, g ~ Inverse-Gamma(1, n/2), averaged over the angle phi
## between the effect's direction and the data's (weight sin(phi)^(k - 2)).
## More precise than the B-mixture above when Beta(1/2, (k - 1)/2) is sharp.
log_direction_average_by_integration <- function(scaled_share, n, k) {
  log_at_zero <- -(n - 1) / 2 * log1p(-scaled_share)
  if (k == 1) return(log_at_zero)
  log_total_weight <- 0.5 * log(pi) - log(2) + lgamma((k - 1) / 2) - lgamma(k / 2)
  integrand <- function(angle) {
    log_value <- -(n - 1) / 2 * log1p(-scaled_share * cos(angle)^2) - log_at_zero +
      (k - 2) * log(sin(angle)) - log_total_weight
    exp(log_value)
  }
  width  <- if (scaled_share > 0) sqrt((1 - scaled_share) / ((n - 1) * scaled_share)) else Inf
  inner  <- c(3, 10, 40) * width
  breaks <- c(0, inner[inner < pi / 2], pi / 2)
  piece <- function(i) {
    integrate(integrand, breaks[i], breaks[i + 1], rel.tol = 1e-10, subdivisions = 1000L)$value
  }
  log_at_zero + log(sum(vapply(seq_len(length(breaks) - 1), piece, numeric(1))))
}
log_bayes_factor_flat_r_by_directions <- function(r_squared, n, k) {
  log_integrand <- function(log_g) {
    vapply(log_g, function(lg) {
      g <- exp(lg)
      log(n / 2) - lg - n / (2 * g) - 0.5 * log1p(g) +
        log_direction_average_by_integration(r_squared * g / (1 + g), n, k)
    }, numeric(1))
  }
  log_integral_around_peak(log_integrand, -20, 40)
}
