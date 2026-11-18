## The examples of flatcor's README (README.Rmd reads this file in), run by
## the package's checks (.github/workflows/check.yml) on every push.
library(flatcor)

## the Bayes factor between any two linear models, nested or not
bayes_factor_models(mpg ~ wt + hp, mpg ~ wt, mtcars, seed = 1)

## a model against the intercept-only model
bayes_factor_models(mpg ~ wt + hp, mpg ~ 1, mtcars, seed = 1)

## estimation
fit <- fit_regression(mpg ~ wt + hp, mtcars)
fit
loo(fit)
pp_check(fit)

## another prior
bayes_factor_models(mpg ~ wt + hp, mpg ~ wt, mtcars, prior = "flat_r", seed = 1)
