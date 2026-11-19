## =====================================================================
##  The Stan model, its log density, and what goes into them
## =====================================================================
## Two builds of the one Stan program (inst/stan/flatcor.stan), both kept as
## files in a cache folder, so they are built once and found by every R
## process that looks there:
##
##   flatcor_model()    the model, for sampling: cmdstanr's executable (or,
##                      with backend = "rstan", rstan's model, kept in memory)
##   flatcor_density()  the model's log density on Stan's unconstrained
##                      scale, for the Bayes factors: BridgeStan's library
##                      file. Bridge sampling evaluates it at many points.
##
## fit_regression() and bayes_factor_models() call these themselves, so daily
## use needs neither. Calling them first (with the same cache folder) gives
## the functions a model already built: a simulation builds once, before its
## workers start, and every worker finds the files. The cache folder is
## flatcor's user cache folder, or getOption("flatcor.cache"), or `cache`.
## flatcor_stan_file() and flatcor_stan_data() show what Stan is given.

compiled_models <- new.env()   # what this R session has loaded, by cache folder

## The folder for this flatcor version's builds: the builds of another version
## are never reused
cache_folder <- function(cache = NULL) {
  if (is.null(cache)) {
    cache <- getOption("flatcor.cache", tools::R_user_dir("flatcor", which = "cache"))
  }
  file.path(cache, as.character(utils::packageVersion("flatcor")))
}

## Evaluate `code` while holding a lock on `folder`, so that R processes
## starting at the same time build into it one after another. The lock is a
## folder (creating one is atomic); one older than 30 minutes is taken to be
## left by a process that stopped, and removed.
with_build_lock <- function(folder, code) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  lock <- file.path(folder, "build.lock")
  repeat {
    if (dir.create(lock, showWarnings = FALSE)) break
    age <- difftime(Sys.time(), file.mtime(lock), units = "mins")
    if (!is.na(age) && age > 30) unlink(lock, recursive = TRUE)
    Sys.sleep(2)
  }
  on.exit(unlink(lock, recursive = TRUE))
  code
}

#' The Stan program flatcor uses
#'
#' The path of \code{flatcor.stan}: one Stan program for the three priors (the
#' prior's shape is data), the posterior and the intercept-only model.
#' @return The path of the file.
#' @examples
#' cat(readLines(flatcor_stan_file()), sep = "\n")
#' @export
flatcor_stan_file <- function() {
  system.file("stan", "flatcor.stan", package = "flatcor", mustWork = TRUE)
}

## The data Stan gets for one model: predictors (a matrix, possibly with no
## columns) and outcome. per_observation: also draw log_lik and y_rep (for loo
## and pp_check; a Bayes factor does not need them).
stan_data_of <- function(predictors, outcome, prior, per_observation = FALSE) {
  list(
    N = length(outcome),
    K = ncol(predictors),
    X = predictors,
    y = outcome,
    prior_only = 0L,
    proper_stand_ins = 0L,
    per_observation = as.integer(per_observation),
    r_squared_shape = prior_shape(prior, ncol(predictors))
  )
}

#' The data Stan gets for a model
#'
#' What \code{fit_regression()} and \code{bayes_factor_models()} give the Stan
#' program (\code{flatcor_stan_file()}) for one formula: the rows with no
#' missing value, the predictor columns (without the intercept), the outcome
#' and the prior's shape.
#' @inheritParams fit_regression
#' @param per_observation Also ask for the log likelihood and predicted
#'   outcome of each observation (as \code{fit_regression()} does, for
#'   \code{loo()} and \code{pp_check()}).
#' @return A list, the data block of \code{flatcor.stan}.
#' @examples
#' str(flatcor_stan_data(mpg ~ wt + hp, mtcars, prior = "flat_r"))
#' @export
flatcor_stan_data <- function(formula, data, prior = "flat_partial_r", per_observation = FALSE) {
  prior <- check_prior(prior, defaulted = missing(prior))
  parts <- predictors_and_outcome(formula, complete_rows(data, list(formula)))
  check_full_rank(parts$predictors)
  stan_data_of(parts$predictors, parts$outcome, prior, per_observation)
}

#' flatcor's Stan model, built once
#'
#' The model \code{fit_regression()} and \code{bayes_factor_models()} sample.
#' With cmdstanr, its executable is built once in the cache folder and found
#' there by every later R session or process; with rstan, it is built once
#' per R session. The functions call this themselves; call it first to build
#' the model before parallel workers start, and pass the same \code{cache}
#' (or set \code{options(flatcor.cache = )}) so they find it.
#' @param backend \code{"cmdstanr"} (the default when CmdStan is installed)
#'   or \code{"rstan"}.
#' @param cache The cache folder: by default flatcor's user cache folder
#'   (\code{tools::R_user_dir("flatcor", "cache")}), or
#'   \code{getOption("flatcor.cache")} when set. Each flatcor version has its
#'   own folder inside it.
#' @return A \code{cmdstanr::CmdStanModel} or an rstan \code{stanmodel}.
#' @examples
#' # needs CmdStan (cmdstanr::install_cmdstan())
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
#'   model <- flatcor_model()
#'   fit_regression(mpg ~ wt + hp, mtcars, prior = "flat_r", model = model, seed = 1)
#' }
#' @export
flatcor_model <- function(backend = default_backend(), cache = NULL) {
  backend <- match.arg(backend, c("cmdstanr", "rstan"))
  folder <- cache_folder(cache)
  key <- paste(backend, folder)
  if (is.null(compiled_models[[key]])) {
    compiled_models[[key]] <- if (backend == "cmdstanr") {
      # cmdstanr reuses the executable when it is there and up to date
      with_build_lock(folder, cmdstanr::cmdstan_model(
        flatcor_stan_file(),
        dir = folder,
        quiet = TRUE
      ))
    } else {
      rstan::stan_model(flatcor_stan_file())
    }
  }
  compiled_models[[key]]
}

#' flatcor's log density, built once, for the Bayes factors
#'
#' Bridge sampling, which gives each model's marginal likelihood in
#' \code{bayes_factor_models()}, evaluates the model's log density at many
#' points. This builds that log density once with BridgeStan, into a library
#' file in the cache folder that every later R session or process loads. It is
#' used with the cmdstanr backend (with rstan, rstan's own fit provides the
#' log density). \code{bayes_factor_models()} calls this itself; call it first
#' to build it before parallel workers start, with the same \code{cache}.
#'
#' Needs the bridgestan R package
#' (\code{install.packages("bridgestan", repos =
#' "https://r-multiverse.r-universe.dev")}). The first build also downloads
#' BridgeStan's C++ sources to \code{~/.bridgestan}, unless the environment
#' variable \code{BRIDGESTAN} points to them.
#' @inheritParams flatcor_model
#' @return A \code{flatcor_density}: the path of the library file
#'   (\code{$library}).
#' @examples
#' # needs CmdStan and bridgestan
#' if (requireNamespace("bridgestan", quietly = TRUE) &&
#'     requireNamespace("cmdstanr", quietly = TRUE) &&
#'     !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
#'   density <- flatcor_density()
#'   bayes_factor_models(mpg ~ wt + hp, mpg ~ wt, mtcars, prior = "flat_r",
#'                       density = density, seed = 1)
#' }
#' @export
flatcor_density <- function(cache = NULL) {
  if (!requireNamespace("bridgestan", quietly = TRUE)) {
    stop("Bayes factors with cmdstanr need the bridgestan package: ",
         "install.packages(\"bridgestan\", repos = \"https://r-multiverse.r-universe.dev\"). ",
         "Or use backend = \"rstan\".", call. = FALSE)
  }
  bridgestan_version <- as.character(utils::packageVersion("bridgestan"))
  folder <- file.path(cache_folder(cache), paste0("bridgestan-", bridgestan_version))
  key <- paste("density", folder)
  if (is.null(compiled_models[[key]])) {
    stan_copy <- file.path(folder, "flatcor.stan")
    library_file <- file.path(folder, "flatcor_model.so")
    with_build_lock(folder, {
      # built again when the Stan program differs from the one it was built from
      unchanged <- file.exists(library_file) && file.exists(stan_copy) &&
        identical(unname(tools::md5sum(stan_copy)), unname(tools::md5sum(flatcor_stan_file())))
      if (!unchanged) {
        file.copy(flatcor_stan_file(), stan_copy, overwrite = TRUE)
        build_with_valid_tmpdir(bridgestan::compile_model(stan_copy))
      }
    })
    compiled_models[[key]] <- structure(list(library = library_file), class = "flatcor_density")
  }
  compiled_models[[key]]
}

## Evaluate `code` with TMPDIR naming a folder that exists: the compiler's
## tools write their temporary files there, and fail when it names a folder
## that does not exist (R itself then falls back quietly to another)
build_with_valid_tmpdir <- function(code) {
  temporary <- Sys.getenv("TMPDIR")
  if (nzchar(temporary) && dir.exists(temporary)) return(code)
  Sys.setenv(TMPDIR = tempdir())
  on.exit(if (nzchar(temporary)) Sys.setenv(TMPDIR = temporary) else Sys.unsetenv("TMPDIR"))
  code
}

#' @export
print.flatcor_density <- function(x, ...) {
  cat("flatcor's log density (BridgeStan library):", x$library, "\n")
  invisible(x)
}

## The model must be one the backend can sample
check_model <- function(model, backend) {
  expected <- if (backend == "cmdstanr") "CmdStanModel" else "stanmodel"
  if (!inherits(model, expected)) {
    stop("model = must be a ", backend, " model (", expected, "), such as flatcor_model(\"",
         backend, "\").", call. = FALSE)
  }
  invisible(TRUE)
}
