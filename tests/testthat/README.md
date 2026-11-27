# flatcor's tests

The test files are numbered in the order to read them. Each one builds on
the ones before it: first the priors, then their draws, then the Stan
model, then fitting, then Bayes factors, and last the second backend
(rstan).

The four priors are `joint_partial_r` (the default), `independent_partial_r`,
`flat_r` and `flat_r_squared`. The tables show which tests cover
`independent_partial_r` specifically. The others are covered throughout.

## Where the tests run

| GitHub CI job | What it runs |
|---|---|
| `check` | R CMD check. The Stan tests skip themselves when CmdStan is missing (`helper-skips.R`). |
| `tests-cmdstanr` | Every test with cmdstanr, CmdStan and bridgestan, including all the Bayes factors. |
| `tests-rstan` | Files 4 and 6 (`filter = "rstan|fit"`) with rstan. |
| `lint` | lintr |

Tests marked *bridgestan* skip on a machine without the bridgestan
package. Only the `tests-cmdstanr` job runs them in full.

## Helpers

- `helper-references.R` gives references for the Bayes factor that share
  no code with the package. The prior is written as a mixture of g-priors
  and integrated with `integrate()`. `flat_r` uses a direction form, which
  is more precise when k is large.
- `helper-skips.R` decides when the Stan tests skip themselves.

## 1. `test-1-priors.R`: the priors themselves

| # | Test | independent |
|---|---|---|
| 1.1 | each prior has its shape, and joint_partial_r is the default | ✓ (it has no shape) |
| 1.2 | the tests' sample R-squared is that of the least-squares fit | |
| 1.3 | independent_partial_r is a prior; flat_partial_r, its old name, is not | ✓ |

## 2. `test-2-prior-draws.R`: the prior as a recipe (`draw_coefficients_from_prior()`)

| # | Test | independent |
|---|---|---|
| 2.1 | flat R: implied R-squared is Beta(1/2, 1) for k = 1, 3, 10, correlated or not | |
| 2.2 | with one predictor the correlation is flat from -1 to 1 | |
| 2.3 | joint flat partial r: implied R-squared is Beta(k/2, 1), correlated or not | |
| 2.4 | joint flat partial r: every coefficient's partial correlation is flat, for any k | |
| 2.5 | flat R-squared: implied R-squared is flat for every k, correlated or not | |
| 2.6 | the default is joint_partial_r, and an unknown prior is refused | |
| 2.7 | the prior draws carry what each set was drawn from | |
| 2.8 | each partial correlation is flat, for independent and joint partial r | ✓ |
| 2.9 | the partial correlations are unlinked for independent, linked for joint | ✓ |
| 2.10 | independent_partial_r carries what it was drawn from | ✓ |
| 2.11 | with one predictor, independent and joint partial r are the same prior | ✓ |
| 2.12 | Stan, sampling the prior only, gives the draws of the recipe | ✓ (both partial-r priors) |
| 2.13 | partial_r is each coefficient's partial correlation, by the textbook definition | ✓ (all four priors) |

Test 2.12 links the recipe to the Stan model. Test 2.13 links the recipe's
`partial_r` to the textbook partial correlation, computed from the precision
matrix of (y, x).

## 3. `test-3-model.R`: the Stan model and its data

| # | Test | independent |
|---|---|---|
| 3.1 | the Stan program and its data are what the functions give Stan | ✓ (the `independent_effects` flag) |
| 3.2 | BridgeStan's parameter names become cmdstanr's | |
| 3.3 | a model of the wrong kind is refused | |
| 3.4 | given the model, or building it themselves, the functions give the same result | *bridgestan* |
| 3.5 | the cache folder can be chosen, and the builds are found there | *bridgestan* |
| 3.6 | the Bayes factor from BridgeStan equals the one from model methods | *bridgestan* |

## 4. `test-4-fit.R`: fitting (`fit_regression()`)

| # | Test | independent |
|---|---|---|
| 4.1 | fit_regression() runs and works with print, loo and pp_check | ✓ (print, among other priors) |
| 4.2 | cmdstanr fits leave no files in R's temporary folder, and keep working | |

## 5. `test-5-bayes-factor.R`: Bayes factors (`bayes_factor_models()`)

| # | Test | independent |
|---|---|---|
| 5.1 | the inputs are checked before anything is fitted | |
| 5.2 | the intercept-only model's marginal likelihood matches its closed form | *bridgestan* |
| 5.3 | every prior's Bayes factor against the intercept-only model matches its reference | *bridgestan*, no reference (see below) |
| 5.4 | rescaling the predictors and the outcome leaves the Bayes factor unchanged | ✓ *bridgestan* |
| 5.5 | a seed, or set.seed() beforehand, makes the Bayes factor reproducible | *bridgestan* |
| 5.6 | everything computed for a Bayes factor comes out, per model | *bridgestan* |
| 5.7 | keep_fits adds the fits and bridge sampling's objects, and changes nothing else | *bridgestan* |
| 5.8 | without a seed, the seed that comes out reproduces the Bayes factor | *bridgestan* |
| 5.9 | with one predictor, independent_partial_r gives the Bayes factor of flat_r | ✓ *bridgestan* |

`independent_partial_r` has no closed-form R-squared for k > 1, so its Bayes
factor has no independent reference there (5.3). Its correctness rests on
four checks together:

- The Stan prior equals the recipe (2.12).
- The recipe equals the textbook partial correlation (2.13).
- Bridge sampling is the same code for every prior (5.3 for the others).
- With one predictor its Bayes factor equals that of `flat_r` (5.9).

## 6. `test-6-rstan.R`: the rstan backend

| # | Test | independent |
|---|---|---|
| 6.1 | fit_regression() and bayes_factor_models() run with rstan | ✓ (an independent fit) |

rstan runs in its own R process, because rstan and cmdstanr cannot share
one R session for Bayes factors.
