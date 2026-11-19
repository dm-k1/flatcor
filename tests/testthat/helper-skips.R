## When the Stan tests skip themselves. R CMD check (the check job) can have
## the cmdstanr and rstan packages installed without CmdStan, and allows at
## most 2 cores (_R_CHECK_LIMIT_CORES_); the test jobs for cmdstanr and rstan
## run these tests in full.

## CmdStan itself, not only the cmdstanr package
skip_without_cmdstan <- function() {
  skip_on_cran()
  skip_if_not(requireNamespace("cmdstanr", quietly = TRUE), "needs cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)), "needs CmdStan")
}

## Bayes factors with cmdstanr take the log density from BridgeStan
skip_without_bridgestan <- function() {
  skip_without_cmdstan()
  skip_if_not(requireNamespace("bridgestan", quietly = TRUE), "needs bridgestan")
}

## 4 chains run in parallel, more than R CMD check allows
skip_if_cores_limited <- function() {
  limit <- toupper(Sys.getenv("_R_CHECK_LIMIT_CORES_"))
  skip_if(nzchar(limit) && limit != "FALSE", "R CMD check allows 2 cores; 4 chains run at once")
}
