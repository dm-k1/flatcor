// One of the three priors, R^2 ~ Beta(a, 1) with a uniform direction, sampled
// in its own coordinates:
//   R          the overall correlation (factor 1): R^2 ~ Beta(a, 1), with
//              a = K/2 (flat_partial_r: every correlation matrix of
//              (y, x1 ... xK) equally likely), 1/2 (flat_r) or 1 (flat_r_squared)
//   direction  how the total is divided, uniform on the sphere (factor 2)
// and the coefficients follow from them (factor 3 is the conversion):
//   b = sigma * f * (direction expressed in predictor units),  f = R / sqrt(1 - R^2)
// In these coordinates the prior has no peak at zero, so the sampler
// reaches small effects for any number of predictors.
//
// The sampler works with the size of R, R_magnitude in (0, 1). With several
// predictors the sign of each coefficient lives in the direction. With one
// predictor the sign of R is + or - with probability 1/2 each: it is
// marginalised in the model block (the likelihood averaged over both signs)
// and drawn in generated quantities. Sampling a signed R on (-1, 1) instead
// would put a zero of the density at R = 0 whenever a > 1/2 (flat_r_squared:
// density |R|), a wall the sampler cannot cross without divergent transitions.
data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;                 // predictors (centered here)
  vector[N] y;
  int<lower=0, upper=1> prior_only;      // 1: ignore y (sample the prior)
  int<lower=0, upper=1> proper_stand_ins; // 1: proper priors on alpha and sigma (for checks)
  real<lower=0> r_squared_shape;         // a in R^2 ~ Beta(a, 1)
}
transformed data {
  matrix[N, K] Xc;
  for (j in 1:K) Xc[, j] = X[, j] - mean(X[, j]);
  // Any square root of the predictor covariance works; Cholesky is the cheapest.
  // It turns a direction in "explained variance" units into predictor units.
  matrix[K, K] root = cholesky_decompose(crossprod(Xc) / N);
}
parameters {
  real alpha;
  real<lower=0> sigma;
  real<lower=0, upper=1> R_magnitude;    // |R|
  // the direction is z / |z|. Any density for z that depends only on its length
  // makes the direction uniform; this one keeps |z| near 1, away from 0, where
  // the direction would be undefined (Stan's unit_vector lets |z| reach 0,
  // which gives divergent transitions)
  vector[K > 1 ? K : 0] z;
}
transformed parameters {
  real f_magnitude = R_magnitude / sqrt(1 - square(R_magnitude));   // |Cohen's f|
  vector[K] direction;
  if (K == 1) direction = rep_vector(1, 1);
  else        direction = z / norm2(z);
  // b' Sigma b = sigma^2 f^2; with one predictor this is |b|, its sign follows
  vector[K] b_magnitude = sigma * f_magnitude * mdivide_right_tri_low(direction', root)';
}
model {
  // factor 1: R^2 ~ Beta(a, 1), which is density 2a R^(2a - 1) on R in (0, 1)
  // (flat for a = 1/2; with K = 1, flat_partial_r and flat_r both have a = 1/2)
  if (r_squared_shape != 0.5) target += (2 * r_squared_shape - 1) * log(R_magnitude);
  // factor 2: direction uniform
  if (K > 1) target += -0.5 * square((norm2(z) - 1) / 0.1);
  if (proper_stand_ins) {
    // Checks that sample the prior, or simulate data from it, need proper priors
    // on alpha and sigma. R, f and the direction don't depend on these.
    alpha ~ normal(0, 10);
    sigma ~ lognormal(0, 1);
  } else {
    target += -log(sigma);                   // Jeffreys prior on the noise; alpha flat
  }
  if (!prior_only) {
    if (K == 1)                              // the sign of R, + or - with probability 1/2
      target += log_mix(0.5, normal_lpdf(y | alpha + Xc * b_magnitude, sigma),
                             normal_lpdf(y | alpha - Xc * b_magnitude, sigma));
    else
      y ~ normal(alpha + Xc * b_magnitude, sigma);
  }
}
generated quantities {
  // the sign of R (K = 1 only), drawn with the probability its likelihood gives it
  real sign_of_R = 1;
  if (K == 1) {
    real log_odds_plus = 0;
    if (!prior_only)
      log_odds_plus = normal_lpdf(y | alpha + Xc * b_magnitude, sigma)
                      - normal_lpdf(y | alpha - Xc * b_magnitude, sigma);
    sign_of_R = bernoulli_logit_rng(log_odds_plus) ? 1 : -1;
  }
  real R = sign_of_R * R_magnitude;          // the overall correlation (signed when K = 1)
  real f = sign_of_R * f_magnitude;          // Cohen's f (signed when K = 1)
  vector[K] b = sign_of_R * b_magnitude;     // the coefficients
  // per observation: log-likelihood (for loo) and a predicted value (for pp_check)
  vector[N] log_lik;
  array[N] real y_rep;
  real R_squared = square(R_magnitude);      // population explained share
  {
    vector[N] mu = alpha + Xc * b;
    for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | mu[n], sigma);
    y_rep = normal_rng(mu, sigma);
  }
}
