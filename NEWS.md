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
