# flatcor 1.0.2

- `log_bayes_factor_conditional()` with `prior = "flat_r"` or
  `"flat_r_squared"` stopped with an error when the partial correlation was
  so near 1 that every effect on its grid exceeded f = 100 (for example
  r = .99999 with n = 1000). It now works.
- With few residual degrees of freedom (n - k < 10) and a partial correlation
  near 1, the conditional test's grid missed the broad tail of the likelihood
  ratio, an error of up to 0.02 on the log scale. It now adds a wide grid,
  as the whole-model test already did; the error is below about 3e-4.
- `bayes_factor_models()` and `bayes_factor_per_term()` stop with a clear
  message for a model without an intercept (before, `bayes_factor_models()`
  dropped the first predictor and returned a wrong Bayes factor).
- `bayes_factor_per_term(test = "joint")` works with a tibble.
- `pp_check()` no longer fails when `ndraws` exceeds the number of draws.

# flatcor 1.0.1

- With one predictor, the Stan model now samples the size of R on (0, 1) and
  marginalises its sign (the likelihood averaged over both signs; the sign is
  drawn in generated quantities). Before, it sampled a signed R on (-1, 1),
  and under `flat_r_squared` the density |R| vanishes at R = 0: the sampler
  then had to cross a zero of the density whenever the sign changed, which
  gave divergent transitions. The prior and the posterior are unchanged (the
  same distributions); `R`, `f` and `b` are now generated quantities, and the
  sampled parameters are `R_magnitude` and, with several predictors, `z`.
- `log_bayes_factor_vs_intercept_only()` documents its precision at
  astronomical evidence (a relative error of about 5e-7 when the log Bayes
  factor is about 5.8 million).
- Comments and help pages point to the reorganised validation repository
  (https://github.com/dm-k1/objective-bayes).

# flatcor 1.0.0

- First release: the three priors (flat partial r, flat R, flat R-squared),
  Bayes factors between linear models, the conditional test per term, and
  estimation with Stan.
