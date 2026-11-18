# flatcor (development version, to be 3.0.0)

One Stan model now does everything, for all three priors: the posterior, and
the Bayes factor between two models. No Bayes factor is computed by a formula
of the package's own any more.

- `bayes_factor_models(formula_a, formula_b, data)` takes the two models as
  formulas, like `fit_regression()`, instead of `lm()` fits: nothing was ever
  taken from an `lm()` fit but its predictors and outcome. Both models use
  the rows with no missing value in any variable of either formula.
- `seed` now makes a Bayes factor fully reproducible: it seeds Stan and also
  R's random numbers that bridge sampling uses. Without `seed`, the seed is
  drawn from R's random numbers, so `set.seed()` beforehand does the same.
- The Stan model draws `log_lik` and `y_rep` (N numbers per draw each) only
  for `fit_regression()`, which needs them for `loo()` and `pp_check()`; a
  Bayes factor fit skips them (data `per_observation`). At n = 5000 with the
  default draws they took about 4.5 GB on disk and 3.2 GB in memory per
  model, enough to crash a worker.
- No more "proposal rejected" messages from Stan. A long trajectory now and
  then overshoots into a region of practically zero probability, where R
  rounded to 1 (an infinite f) or sigma to 0 or infinity, and a density
  function stopped with an error that Stan printed before rejecting the
  proposal. Such proposals now get log density -infinity instead, which
  Stan rejects the same way without a message: the draws and divergent
  transitions are identical to before, bit for bit.
- A Bayes factor no longer fails after flatcor is reloaded in the same R
  session ("Model methods cannot be used with a pre-compiled Stan
  executable"): the version of the model with methods is compiled afresh
  once per loaded package.
- The README's examples live in `inst/examples/quick_start.R`, which
  `README.Rmd` reads in; the help pages' examples run whenever CmdStan is
  installed (no longer `\dontrun`). The checks run both on every push.
- The default backend is cmdstanr only when CmdStan itself is installed;
  with just the cmdstanr R package it is rstan (before, fitting stopped with
  "CmdStan path has not been set yet").
- rstan and cmdstanr cannot both be used for Bayes factors in one R session;
  the help says so, and the tests run rstan in its own R process.
- `bayes_factor_models()` fits each model with the chosen prior and estimates
  its marginal likelihood by bridge sampling (bridgesampling, now imported);
  the Bayes factor is their ratio. It returns a `flatcor_bayes_factor` with
  the Bayes factor, its log, bridgesampling's approximate error of the log,
  both log marginal likelihoods and each fit's divergent transitions. By
  default each model gets 4 chains of 10,000 draws (ten times those for
  estimation, as the bridgesampling authors advise for testing); `seed`,
  `method = "warp3"` and the sampler settings can be set. It needs Stan
  (cmdstanr or rstan).
- A model with no predictors (`lm(y ~ 1)`) is the intercept-only model: the
  Stan model now accepts K = 0.
- The Stan model keeps every normalising constant (`target +=` throughout, the
  2a of R's density and the constant of the direction's density), as bridge
  sampling needs.
- Removed `log_bayes_factor_vs_intercept_only()` and the likelihood ratio code
  behind it: it went wrong with many predictors (about 90 or more), at
  near-perfect fits with few observations and at very strong evidence.
- Removed the support curves (`support_per_partial_correlation()`,
  `support_per_partial_r_squared()`, `support_regions()`): they use no prior,
  and ran on the same likelihood ratio code.
- The compiled Stan model is kept in a folder per flatcor version, so an
  executable from another version is never reused. Bayes factors need the
  model's log density as well, which cmdstanr can only provide in the R
  session that compiled the model: that version is compiled once per session,
  into the session's temporary folder.
- The help pages and the package description no longer say that the
  validation recommends flat partial r: the validation is being redone, and
  the default stays flat partial r.
- `stats::dt()` and other unused functions are no longer imported; ggplot2 is
  no longer suggested.

# flatcor 2.0.0

flatcor now does two things, the same way for all three priors: the Bayes
factor between two models, and the posterior.

- Removed `log_bayes_factor_conditional()`, the conditional test of one
  coefficient or term. It is a Bayes factor of one model only under flat
  partial r (it equals the Savage–Dickey ratio of the full model); under flat
  R and flat R-squared it is not, so it cannot be the test for all three.
- Removed `bayes_factor_per_term()`. The evidence for a term is the model with
  it against the model without it: `bayes_factor_models(full, reduced)`.
- The support curves (`support_per_partial_correlation()`,
  `support_per_partial_r_squared()`, `support_regions()`) stay. They use no
  prior, so they are the same under all three priors; their help pages no
  longer refer to the conditional test.
- One route computes the whole-model Bayes factor for all three priors: the
  likelihood ratio of the model's columns against none, averaged over
  R-squared ~ Beta(a, 1) by `integrate()` over log f, centred on the peak (as
  BayesFactor computes its regression Bayes factor over log g). Only a
  differs between the priors. It replaces three prior-specific closed forms,
  which had grid errors of up to 3e-5 (flat R) and 5e-4 (flat R-squared at
  n = 10^6); against adaptive integration of flat partial r's g-prior form
  the new route is within 2e-7 on 174 settings. It is slower: about 0.2 s
  per Bayes factor with more than one predictor (10 ms with one).
- The likelihood ratio behind the support curves and the Bayes factors moved
  to `R/likelihood.R`. Its inner integral now starts at zero and uses
  Simpson's rule, which keeps it smooth with one residual degree of freedom
  (within 6e-8 of the exact hypergeometric series), and `besselI()` is only
  called where its value is used (it is slow at large arguments).

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
  (https://github.com/dm-k1/flatcor-validation).

# flatcor 1.0.0

- First release: the three priors (flat partial r, flat R, flat R-squared),
  Bayes factors between linear models, the conditional test per term, and
  estimation with Stan.
