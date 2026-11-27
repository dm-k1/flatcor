## =====================================================================
##  Example: flatcor on one dataset, start to finish
##  Estimation, then Bayes factors for model comparisons, with the four
##  priors side by side: joint_partial_r, independent_partial_r, flat_r and
##  flat_r_squared. The
##  default is joint_partial_r; prior = selects the others (the validation:
##  https://github.com/dm-k1/flatcor-validation). Everything runs on Stan.
## =====================================================================
library(flatcor)   # install once: remotes::install_github("dm-k1/flatcor")

## ---- 1. Data ---------------------------------------------------------
## Three correlated predictors (x1 and x2 matter, x3 doesn't) and one
## unrelated extra predictor x4 for a non-nested comparison.
set.seed(42)
sample_size <- 60
predictor_correlations <- matrix(0.4, nrow = 3, ncol = 3)
diag(predictor_correlations) <- 1
predictor_values <- matrix(rnorm(sample_size * 3), sample_size) %*% chol(predictor_correlations)
my_data <- data.frame(x1 = predictor_values[, 1],
                      x2 = predictor_values[, 2],
                      x3 = predictor_values[, 3])
my_data$y  <- 0.3 * my_data$x1 + 0.2 * my_data$x2 + rnorm(sample_size)
my_data$x4 <- rnorm(sample_size)

## ---- 2. Estimation ---------------------------------------------------
full_model <- fit_regression(y ~ x1 + x2 + x3, data = my_data, prior = "joint_partial_r",
                             chains = 4, iter_warmup = 1000, iter_sampling = 1000, seed = 1)
full_model                       # estimates, intervals, Rhat, ESS, divergences
loo(full_model)                  # leave-one-out cross-validation
pp_check(full_model)             # posterior predictive check
bayesplot::mcmc_trace(full_model$draws)   # trace plots (any bayesplot plot works on $draws)

## The four priors side by side: at n = 60 with three predictors each prior
## pulls the estimates toward zero, by a different amount
priors <- c("joint_partial_r", "independent_partial_r", "flat_r", "flat_r_squared")
other_fits <- lapply(c(independent_partial_r = "independent_partial_r",
                       flat_r = "flat_r",
                       flat_r_squared = "flat_r_squared"), function(prior) {
  fit_regression(y ~ x1 + x2 + x3, data = my_data, prior = prior, chains = 4,
                 iter_warmup = 1000, iter_sampling = 1000, seed = 1)
})
posterior_means <- function(fit) colMeans(as.data.frame(fit$draws)[, c("x1", "x2", "x3")])
cbind(joint_partial_r       = posterior_means(full_model),
      independent_partial_r = posterior_means(other_fits$independent_partial_r),
      flat_r                = posterior_means(other_fits$flat_r),
      flat_r_squared        = posterior_means(other_fits$flat_r_squared))

## ---- 3. Comparing models ------------------------------------------
## Takes two formulas and the data; each model is fitted with the prior and
## its marginal likelihood estimated by bridge sampling, so every Bayes
## factor comes with a Monte Carlo error
full <- y ~ x1 + x2 + x3
# evidence for x1
bayes_factor_models(full, y ~ x2 + x3, my_data, prior = "joint_partial_r", seed = 1)
# non-nested
bayes_factor_models(full, y ~ x4, my_data, prior = "joint_partial_r", seed = 1)
# against the intercept-only model
bayes_factor_models(full, y ~ 1, my_data, prior = "joint_partial_r", seed = 1)

## Evidence for each predictor under each prior: the model with it against
## the model without it
without <- function(predictor) update(full, paste(". ~ . -", predictor))
priors_by_predictor <- sapply(priors, function(prior) {
  sapply(c("x1", "x2", "x3"), function(predictor) {
    bayes_factor_models(full, without(predictor), my_data, prior = prior, seed = 1)$bayes_factor
  })
})
priors_by_predictor   # BF10 per predictor (rows) and prior (columns)

## The Monte Carlo error directly: the same Bayes factor with other seeds
sapply(1:5, function(seed) {
  result <- bayes_factor_models(full, without("x1"), my_data, prior = "joint_partial_r",
                                seed = seed)
  result$log_bayes_factor
})

## ---- 4. Comparison with the default (JZS / JASP) --------------------
## Optional; needs the BayesFactor package
if (requireNamespace("BayesFactor", quietly = TRUE)) {
  default_bf_x1 <- BayesFactor::lmBF(y ~ x1 + x2 + x3, data = my_data) /
    BayesFactor::lmBF(y ~ x2 + x3, data = my_data)
  print(default_bf_x1)
}
