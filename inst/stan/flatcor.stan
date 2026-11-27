// One of four priors, sampled in its own coordinates. The first,
// independent_partial_r, does not go through R at all: each
// predictor's partial correlation with y (given the other predictors) gets
// its own flat prior on (-1, 1), not linked to the others (data
// independent_effects = 1). The sampler works with those partial
// correlations directly, and the coefficients follow from them:
//   b_j = sigma * rho_j / (s_j * sqrt(1 - rho_j^2)),  s_j = residual SD of x_j
//                                                     given the other predictors
// R and f are then what these coefficients imply.
//
// The other three are R^2 ~ Beta(a, 1) with a uniform direction:
//   R          the overall correlation (factor 1): R^2 ~ Beta(a, 1), with
//              a = K/2 (joint_partial_r: one prior over every valid correlation
//              matrix of (y, x1 ... xK); each partial r is flat on its own,
//              but they are linked), 1/2 (flat_r) or 1 (flat_r_squared)
//   direction  how the total is divided, uniform on the sphere (factor 2)
// and the coefficients follow from them (factor 3 is the conversion):
//   b = sigma * f * (direction expressed in predictor units),  f = R / sqrt(1 - R^2)
// In these coordinates the prior has no peak at zero, so the sampler
// reaches small effects for any number of predictors.
//
// For the three R-squared ~ Beta(a, 1) priors the sampler works with the size of R,
// R_magnitude in (0, 1). With several predictors the sign of each coefficient lives in the direction. With one
// predictor the sign of R is + or - with probability 1/2 each: it is
// marginalised in the model block (the likelihood averaged over both signs)
// and drawn in generated quantities. Sampling a signed R on (-1, 1) instead
// would put a zero of the density at R = 0 whenever a > 1/2 (flat_r_squared:
// density |R|), a wall the sampler cannot cross without divergent transitions.
// independent_partial_r needs no such mixture: the sign is in partial_r.
//
// K = 0 is the intercept-only model: only alpha and sigma.
//
// Every density is written with target += and all its normalising constants,
// because bayes_factor_models() estimates each model's marginal likelihood
// from this program with bridgesampling, which needs them (Gronau, Singmann &
// Wagenmakers, 2020). The intercept (flat) and the noise (1/sigma) have
// improper priors; they are the same in every model, so their arbitrary
// constant cancels in a Bayes factor.
functions {
  // the radial part of the direction's auxiliary density (see parameters),
  // r^(K - 1) exp(-((r - 1) / 0.1)^2 / 2), for its normalising constant
  real radial_density(real r, real r_complement, array[] real theta,
                      array[] real x_r, array[] int x_i) {
    return exp((x_i[1] - 1) * log(r) - 0.5 * square((r - 1) / 0.1));
  }
}
data {
  int<lower=1> N;
  int<lower=0> K;
  matrix[N, K] X;                 // predictors (centered here)
  vector[N] y;
  int<lower=0, upper=1> prior_only;      // 1: ignore y (sample the prior)
  int<lower=0, upper=1> proper_stand_ins; // 1: proper priors on alpha and sigma (for checks)
  int<lower=0, upper=1> per_observation; // 1: log_lik and y_rep per observation (loo,
                                         // pp_check); 0: skip them (Bayes factors, where
                                         // they would be N numbers per draw for nothing)
  real<lower=0> r_squared_shape;         // a in R^2 ~ Beta(a, 1) (unused when K = 0 or
                                         // independent_effects = 1)
  int<lower=0, upper=1> independent_effects; // 1: independent_partial_r, 0: the three
                                         // R-based priors
}
transformed data {
  matrix[N, K] Xc;
  for (j in 1:K) Xc[, j] = X[, j] - mean(X[, j]);
  // Any square root of the predictor covariance works; Cholesky is the cheapest.
  // It turns a direction in "explained variance" units into predictor units.
  matrix[K, K] predictor_covariance;
  matrix[K, K] root;
  // what each predictor has that the others don't: the SD of x_j left after
  // the other predictors, 1 / sqrt(diag(Sigma^-1))_j (independent_partial_r)
  vector[K] residual_sd;
  if (K > 0) {
    predictor_covariance = crossprod(Xc) / N;
    root = cholesky_decompose(predictor_covariance);
    residual_sd = 1 ./ sqrt(diagonal(inverse_spd(predictor_covariance)));
  }
  // log of the normalising constant of z's density on R^K (K > 1): the
  // surface area of the unit sphere, 2 pi^(K/2) / Gamma(K/2), times the
  // integral of the radial part
  real log_direction_constant = 0;
  if (K > 1)
    log_direction_constant = log2() + K / 2.0 * log(pi()) - lgamma(K / 2.0) +
      log(integrate_1d(radial_density, 0, positive_infinity(), {0.0}, {0.0}, {K}, 1e-10));
}
parameters {
  real alpha;
  real<lower=0> sigma;
  vector<lower=0, upper=1>[K > 0 && !independent_effects ? 1 : 0] R_size;   // |R|, absent when
                                                    // K = 0 or independent_effects
  // the direction is z / |z|. Any density for z that depends only on its length
  // makes the direction uniform; this one keeps |z| near 1, away from 0, where
  // the direction would be undefined (Stan's unit_vector lets |z| reach 0,
  // which gives divergent transitions)
  vector[K > 1 && !independent_effects ? K : 0] z;
  // independent_partial_r: each predictor's partial correlation with y, flat on its own
  vector<lower=-1, upper=1>[independent_effects ? K : 0] partial_r;
}
transformed parameters {
  real R_magnitude = 0;
  real f_magnitude = 0;                          // |Cohen's f|
  vector[K] direction = rep_vector(0, K);
  // b' Sigma b = sigma^2 f^2. For the R-based priors this is |b| with one
  // predictor, and its sign follows; with independent_partial_r, b_magnitude
  // is already signed (the sign is in partial_r)
  vector[K] b_magnitude = rep_vector(0, K);
  if (K > 0 && !independent_effects) {
    R_magnitude = R_size[1];
    f_magnitude = R_magnitude / sqrt(1 - square(R_magnitude));
    if (K == 1) direction = rep_vector(1, 1);
    if (K > 1)  direction = z / norm2(z);
    b_magnitude = sigma * f_magnitude * mdivide_right_tri_low(direction', root)';
  }
  if (K > 0 && independent_effects) {
    // step 1: each partial correlation is a parameter, partial_r
    // step 2: residual_sd (transformed data)
    // step 3: each rho_j to its coefficient
    b_magnitude = sigma * partial_r ./ sqrt(1 - square(partial_r)) ./ residual_sd;
    // step 4: what this implies for the total
    f_magnitude = sqrt(quad_form(predictor_covariance, b_magnitude)) / sigma;
    R_magnitude = f_magnitude / sqrt(1 + square(f_magnitude));
    // the direction, in the same units as for the other priors
    direction = (root' * b_magnitude) / norm2(root' * b_magnitude);
  }
}
model {
  // factor 1: R^2 ~ Beta(a, 1), which is density 2a R^(2a - 1) on R in (0, 1)
  if (K > 0 && !independent_effects)
    target += log(2 * r_squared_shape) + (2 * r_squared_shape - 1) * log(R_magnitude);
  // factor 2: direction uniform (z's density, normalised)
  if (K > 1 && !independent_effects)
    target += -0.5 * square((norm2(z) - 1) / 0.1) - log_direction_constant;
  // independent_partial_r: each rho_j ~ U(-1, 1) has density 1/2, so 2^-K in all
  // (the full constant, which bridge sampling needs)
  if (independent_effects) target += -K * log2();

  // A long trajectory sometimes overshoots far into a region of practically
  // zero probability, where R rounds to 1 (an infinite f), sigma rounds to 0
  // or to infinity, or the mean of y is not finite. The density there is
  // zero: such a proposal is rejected by giving it log density -infinity,
  // instead of letting a density function stop with an error, which Stan
  // would print and then reject the proposal just the same.
  vector[N] mu = alpha + Xc * b_magnitude;
  real mu_sum = sum(mu);
  if (sigma == 0 || is_inf(sigma)
      || (!prior_only && (is_inf(mu_sum) || is_nan(mu_sum)))) {
    target += negative_infinity();
  } else {
    if (proper_stand_ins) {
      // Checks that sample the prior, or simulate data from it, need proper
      // priors on alpha and sigma. R, f and the direction don't depend on these.
      target += normal_lpdf(alpha | 0, 10) + lognormal_lpdf(sigma | 0, 1);
    } else {
      target += -log(sigma);                 // Jeffreys prior on the noise; alpha flat
    }
    if (!prior_only) {
      if (K == 1 && !independent_effects)    // the sign of R, + or - with probability 1/2
        target += log_mix(0.5, normal_lpdf(y | mu, sigma),
                               normal_lpdf(y | alpha - Xc * b_magnitude, sigma));
      else
        target += normal_lpdf(y | mu, sigma);
    }
  }
}
generated quantities {
  // the sign of R (K = 1 only), drawn with the probability its likelihood gives it
  real sign_of_R = 1;
  if (K == 1 && independent_effects) sign_of_R = partial_r[1] >= 0 ? 1 : -1;
  if (K == 1 && !independent_effects) {
    real log_odds_plus = 0;
    if (!prior_only)
      log_odds_plus = normal_lpdf(y | alpha + Xc * b_magnitude, sigma)
                      - normal_lpdf(y | alpha - Xc * b_magnitude, sigma);
    sign_of_R = bernoulli_logit_rng(log_odds_plus) ? 1 : -1;
  }
  real R = sign_of_R * R_magnitude;          // the overall correlation (signed when K = 1)
  real f = sign_of_R * f_magnitude;          // Cohen's f (signed when K = 1)
  // the coefficients (with independent_partial_r b_magnitude is already signed)
  vector[K] b = independent_effects ? b_magnitude : sign_of_R * b_magnitude;
  // per observation: log-likelihood (for loo) and a predicted value (for
  // pp_check); empty when per_observation = 0
  vector[per_observation ? N : 0] log_lik;
  array[per_observation ? N : 0] real y_rep;
  real R_squared = square(R_magnitude);      // population explained share
  if (per_observation) {
    vector[N] mu = alpha + Xc * b;
    for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | mu[n], sigma);
    y_rep = normal_rng(mu, sigma);
  }
}
