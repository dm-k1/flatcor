## =====================================================================
##  Example: flatcor on one dataset, start to finish
##  Estimation, then Bayes factors for model comparisons, with the three
##  priors side by side: flat_partial_r, flat_r and flat_r_squared. The
##  default is flat_partial_r; prior = selects the others (the validation:
##  https://github.com/dm-k1/objective-bayes). Everything runs on Stan.
## =====================================================================
library(flatcor)   # install once: remotes::install_github("dm-k1/flatcor")

## ---- 1. Data ---------------------------------------------------------
## Three correlated predictors (x1 and x2 matter, x3 doesn't) and one
## unrelated extra predictor x4 for a non-nested comparison.
set.seed(42)
sample_size <- 60
predictor_correlations <- matrix(0.4, nrow = 3, ncol = 3); diag(predictor_correlations) <- 1
predictor_values <- matrix(rnorm(sample_size * 3), sample_size) %*% chol(predictor_correlations)
my_data <- data.frame(x1 = predictor_values[, 1],
                      x2 = predictor_values[, 2],
                      x3 = predictor_values[, 3])
my_data$y  <- 0.3 * my_data$x1 + 0.2 * my_data$x2 + rnorm(sample_size)
my_data$x4 <- rnorm(sample_size)

## ---- 2. Estimation ---------------------------------------------------
full_model <- fit_regression(y ~ x1 + x2 + x3, data = my_data, prior = "flat_partial_r",
                             chains = 4, iter_warmup = 1000, iter_sampling = 1000, seed = 1)
full_model                       # estimates, intervals, Rhat, ESS, divergences
loo(full_model)                  # leave-one-out cross-validation
pp_check(full_model)             # posterior predictive check
bayesplot::mcmc_trace(full_model$draws)   # trace plots (any bayesplot plot works on $draws)

## The three priors next to the ordinary estimates: at n = 60 with three
## predictors each prior pulls the estimates toward zero, by a different amount
priors <- c("flat_partial_r", "flat_r", "flat_r_squared")
other_fits <- lapply(c(flat_r = "flat_r", flat_r_squared = "flat_r_squared"), function(prior)
  fit_regression(y ~ x1 + x2 + x3, data = my_data, prior = prior, chains = 4,
                 iter_warmup = 1000, iter_sampling = 1000, seed = 1))
posterior_means <- function(fit) colMeans(as.data.frame(fit$draws)[, c("x1", "x2", "x3")])
cbind(flat_partial_r = posterior_means(full_model),
      flat_r         = posterior_means(other_fits$flat_r),
      flat_r_squared = posterior_means(other_fits$flat_r_squared),
      ordinary       = coef(lm(y ~ x1 + x2 + x3, data = my_data))[c("x1", "x2", "x3")])

## ---- 3. Comparing models ------------------------------------------
## Takes lm() fits (or fit_regression() fits); each model is fitted with the
## prior and its marginal likelihood estimated by bridge sampling, so every
## Bayes factor comes with a Monte Carlo error
full_lm <- lm(y ~ x1 + x2 + x3, data = my_data)
bayes_factor_models(full_lm, lm(y ~ x2 + x3, data = my_data), prior = "flat_partial_r", seed = 1)  # evidence for x1
bayes_factor_models(full_lm, lm(y ~ x4, data = my_data), prior = "flat_partial_r", seed = 1)       # non-nested
bayes_factor_models(full_lm, lm(y ~ 1, data = my_data), prior = "flat_partial_r", seed = 1)        # against the intercept-only model

## Evidence for each predictor under each prior: the model with it against
## the model without it
without <- function(predictor) update(full_lm, paste(". ~ . -", predictor))
priors_by_predictor <- sapply(priors, function(prior)
  sapply(c("x1", "x2", "x3"), function(predictor)
    bayes_factor_models(full_lm, without(predictor), prior = prior, seed = 1)$bayes_factor))
priors_by_predictor   # BF10 per predictor (rows) and prior (columns)

## The Monte Carlo error directly: the same Bayes factor with other seeds
sapply(1:5, function(seed)
  bayes_factor_models(full_lm, without("x1"), prior = "flat_partial_r", seed = seed)$log_bayes_factor)

## ---- 4. Comparison with the default (JZS / JASP) --------------------
## Optional; needs the BayesFactor package
if (requireNamespace("BayesFactor", quietly = TRUE)) {
  default_bf_x1 <- BayesFactor::lmBF(y ~ x1 + x2 + x3, data = my_data) /
                   BayesFactor::lmBF(y ~ x2 + x3, data = my_data)
  print(default_bf_x1)
}
