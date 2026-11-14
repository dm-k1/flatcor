## =====================================================================
##  The likelihood ratio of an effect's size against zero
## =====================================================================
## No prior enters here. For q columns of a model with k predictor columns
## (q = k: the whole model; q = 1: one coefficient), f is the length of their
## effect in noise units, with its direction uniform within the q columns. The
## other coefficients and the intercept get flat priors and the noise 1/sigma.
## Integrating them out leaves the likelihood ratio
##   p(y | f) / p(y | 0) = exp(-n f^2/2) E_c[I(sqrt(n) r f c)] / I(0),
##   I(b) = integral_0^inf v^m exp(-v^2/2 + b v) dv,   m = n - k + q - 2,
## with r^2 the columns' sample partial R^2 given the other columns and c the
## cosine between the effect and the data's direction (q = 1: c = +-1, the two
## signs). It depends on the data only through r^2, n, k and q.
##
## Two uses: the whole model's Bayes factor against the intercept-only model is
## this ratio with q = k averaged over the prior on f (R/bayes_factor.R), and
## the support curves show it as it is (R/support_curve.R).
## objective-bayes/docs/appendix/A_derivations.md.

## log I(b) for a vector of b: a grid around the peak of the integrand
log_peak_integral <- function(b, m, half_width = 12, points = 241) {
  peak    <- (b + sqrt(b^2 + 4 * m)) / 2
  width   <- 1 / sqrt(1 + m / peak^2)
  offsets <- seq(-half_width, half_width, length.out = points)
  v <- matrix(peak, length(b), points) + outer(width, offsets)
  v[v <= 0] <- NA
  log_integrand <- m * log(v) - v^2 / 2 + b * v
  largest <- apply(log_integrand, 1, max, na.rm = TRUE)
  largest + log(rowSums(exp(log_integrand - largest), na.rm = TRUE)) +
    log(width * (offsets[2] - offsets[1]))
}

## Block of q columns: the effect is a vector gamma in noise units (whitened by
## the block's covariance given the other predictors), |gamma| = f, with a
## uniform direction. Integrating out sigma as for one column leaves
## exp(-n f^2/2) E_c[I(sqrt(n) r f c)] / I(0), with r^2 the block's sample
## partial R^2, c the cosine between gamma and the data's direction, and
## m = n - k + q - 2. The average over directions is the 0F1 below:
##   E_c[exp(t c)] = 0F1(; q/2; t^2/4) = Gamma(q/2) (t/2)^(1 - q/2) I_{q/2-1}(t)
## (Bessel I; cosh(t) when q = 1, the two signs). Its log, for t >= 0:
## Above t = 1e4, besselI() loses its accuracy and returns 0 beyond about 1e5
## (very strong evidence), so there the large-argument expansion
##   I_v(t) e^(-t) = (2 pi t)^(-1/2) (1 - (4v^2 - 1) / (8t) + (4v^2 - 1)(4v^2 - 9) / (2 (8t)^2) - ...)
## is used; at t = 1e4 its error is below 1e-12.
log_scaled_bessel_i <- function(t, order) {
  mu <- 4 * order^2
  ## besselI() only where its value is used: at large arguments it is slow
  ## (NA, for grid points left out, stays NA)
  result <- rep(NA_real_, length(t))
  large  <- !is.na(t) & t > 1e4
  small  <- !is.na(t) & !large
  result[small] <- log(besselI(t[small], order, expon.scaled = TRUE))
  t_large <- t[large]
  result[large] <- -0.5 * log(2 * pi * t_large) +
    log1p(-(mu - 1) / (8 * t_large) + (mu - 1) * (mu - 9) / (2 * (8 * t_large)^2) -
            (mu - 1) * (mu - 9) * (mu - 25) / (6 * (8 * t_large)^3))
  result
}
log_direction_average <- function(t, number_of_columns) {
  if (number_of_columns == 1) return(t + log1p(exp(-2 * t)) - log(2))
  order <- number_of_columns / 2 - 1
  small <- t < 1e-6
  t_large <- ifelse(small, 1, t)
  ifelse(small, t^2 / (2 * number_of_columns),
         lgamma(number_of_columns / 2) - order * log(t_large / 2) +
           log_scaled_bessel_i(t_large, order) + t_large)
}

## log of integral_0^inf v^m exp(-v^2/2) E_c[exp(b v c)] dv for a vector of
## b >= 0: a grid around the peak of the integrand, as log_peak_integral().
## For large b v the direction average grows like exp(b v) v^(-(q-1)/2), so the
## peak is that of I(b) with m lowered by (q - 1)/2.
log_peak_integral_block <- function(b, m, number_of_columns, half_width = 12, points = 241) {
  m_peak  <- pmax(1, m - (number_of_columns - 1) / 2 * (b > 0))
  peak    <- (b + sqrt(b^2 + 4 * m_peak)) / 2
  width   <- 1 / sqrt(1 + m_peak / peak^2)
  ## the grid runs from max(0, peak - half_width widths) to peak + half_width
  ## widths, so it starts exactly at 0 when the peak is near 0 (few residual
  ## degrees of freedom): the result then changes smoothly with b, which
  ## integrate() needs. Simpson's rule (points is odd); the integrand is 0 at
  ## v = 0 for m > 0.
  lower <- pmax(0, peak - half_width * width)
  upper <- peak + half_width * width
  step  <- (upper - lower) / (points - 1)
  v <- matrix(lower, length(b), points) + outer(step, seq(0, points - 1))
  log_integrand <- m * log(v) - v^2 / 2 +
    matrix(log_direction_average(as.vector(b * v), number_of_columns), length(b), points)
  log_integrand[!is.finite(log_integrand)] <- -Inf
  simpson <- c(1, rep(c(4, 2), (points - 3) / 2), 4, 1) / 3
  log_weights <- matrix(log(simpson), length(b), points, byrow = TRUE)
  largest <- apply(log_integrand, 1, max)
  largest + log(rowSums(exp(log_integrand + log_weights - largest))) + log(step)
}

## Log likelihood ratio p(y | f) / p(y | 0) for a block of q columns, for a
## vector of f: depends only on the block's sample partial R^2, n, k and q.
## (In closed form exp(-n f^2/2) 1F1((m+1)/2; q/2; n r^2 f^2 / 2), which is hard
## to evaluate for large n; the integral over v is not.)
log_likelihood_ratio_block <- function(f_values, partial_r_squared, sample_size,
                                       number_of_predictors, number_of_columns) {
  m <- sample_size - number_of_predictors + number_of_columns - 2
  b <- sqrt(sample_size * partial_r_squared) * f_values
  -sample_size * f_values^2 / 2 + log_peak_integral_block(b, m, number_of_columns) -
    log_peak_integral_block(0, m, number_of_columns)
}
