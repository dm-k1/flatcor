# flatcor

Objective Bayes factors and estimation for linear regression, with priors
that extend Jeffreys' flat correlation to several predictors.

With one predictor, the natural objective prior is flat on the correlation.
With several, not every correlation-based quantity can be flat at once, so
the package offers three priors from one family, *R*² ~ Beta(*a*, 1) with a
uniform direction:

| `prior =` | *R*² prior | Kept flat as predictors are added |
|---|---|---|
| `"flat_partial_r"` (default) | Beta(*k*/2, 1) | each partial correlation (LKJ(1); `stan_lm` with η = 1) |
| `"flat_r"` | Beta(½, 1) | the multiple correlation *R* |
| `"flat_r_squared"` | Beta(1, 1) | *R*² |

The default, flat partial r, is the prior the validation recommends: its test
of one coefficient confirms true nulls at the same rate whatever the number
of other predictors, and it is the Bayes factor of one joint model. The
validation, the comparison of the three priors and the arguments behind them
are in a separate repository,
[dm-k1/objective-bayes](https://github.com/dm-k1/objective-bayes).

## Install

```r
remotes::install_github("dm-k1/flatcor")
```

Bayes factors need nothing beyond the package. Estimation needs Stan, through
[cmdstanr](https://mc-stan.org/cmdstanr/) (recommended) or rstan.

## Use

```r
library(flatcor)

## Bayes factor per term (the conditional test for each coefficient or factor)
bayes_factor_per_term(mpg ~ wt + hp + qsec, mtcars)

## any two linear models, nested or not
bayes_factor_models(lm(mpg ~ wt + hp, mtcars), lm(mpg ~ wt, mtcars))

## estimation with Stan
fit <- fit_regression(mpg ~ wt + hp, mtcars)
fit; loo(fit); pp_check(fit)

## another prior
bayes_factor_per_term(mpg ~ wt + hp + qsec, mtcars, prior = "flat_r")
```

Without `prior =` the functions use flat partial r and say so once per
session. Rows with a missing value in a variable of the formula are left out,
with a warning that gives the sample size used. A worked example on one
dataset, start to finish: `system.file("examples", "example.R", package = "flatcor")`.

## Tests

```r
devtools::test()
```

## License

MIT, D. M. Kloppenburg.
