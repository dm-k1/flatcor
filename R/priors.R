## =====================================================================
##  The four priors: one with independent partial correlations,
##  and three R-squared ~ Beta(a, 1) with a uniform direction
## =====================================================================
## With one predictor, Jeffreys' flat correlation is the only choice. With k
## predictors, not everything can stay flat (flatcor-paper/docs/1_why_a_choice.md): keeping each
## partial correlation flat makes R-squared grow with k, keeping the total flat
## makes each partial correlation shrink with k. The options compared:
##
##   independent_partial_r
##                    each predictor's partial correlation with y gets its own
##                    flat prior on (-1, 1), not linked to the others (joint
##                    density 2^-k on the cube). Nothing keeps their total in
##                    check, so R-squared is whatever they imply: it has no
##                    closed form for k > 1, depends on how correlated the
##                    predictors are, and is shown by simulation
##                    (draw_coefficients_from_prior())
##   joint_partial_r  one prior over all valid correlation matrices of
##                    (y, x1 ... xk), every one equally likely. Each partial
##                    correlation is flat on (-1, 1) on its own, but they are
##                    linked: together they must make a valid matrix
##                    R-squared ~ Beta(k/2, 1)
##   flat_r           the multiple correlation R flat on (0, 1) for any k
##                    R-squared ~ Beta(1/2, 1)
##   flat_r_squared   R-squared flat on (0, 1) for any k
##                    R-squared ~ Beta(1, 1)
##
## The three R-squared ~ Beta(a, 1) priors divide the total over the predictors with no
## preference (a uniform direction), so they differ only in the shape a of
## R-squared ~ Beta(a, 1). On R that is the density 2a R^(2a - 1). With one
## predictor, independent_partial_r, joint_partial_r and flat_r are all flat
## on r.
## The default is joint_partial_r. A function that uses it because no prior
## was given says so, once per session, so the choice stays visible. The
## validation of the priors: https://github.com/dm-k1/flatcor-simulation-studies.

prior_names   <- c("independent_partial_r", "joint_partial_r", "flat_r", "flat_r_squared")
prior_label_text <- c("independent flat partial r", "joint flat partial r",
                      "flat R", "flat R-squared")

#' Readable labels of the four priors
#'
#' One short label per prior, in the order of the prior names, for tables
#' and plot legends.
#' @return A named character vector: the names are the values of the
#'   \code{prior} argument (\code{"independent_partial_r"},
#'   \code{"joint_partial_r"}, \code{"flat_r"},
#'   \code{"flat_r_squared"}), the values their labels.
#' @examples
#' prior_labels()
#' prior_labels()[["flat_r"]]
#' @export
prior_labels <- function() {
  stopifnot(length(prior_label_text) == length(prior_names))
  stats::setNames(prior_label_text, prior_names)
}
default_prior <- "joint_partial_r"
default_prior_announced <- new.env()

## Is this the prior with independent partial correlations (not R-squared ~
## Beta(a, 1))? The Stan model and the prior draws treat it apart.
prior_is_independent <- function(prior) {
  check_prior(prior) == "independent_partial_r"
}

## The prior to use: the default (with a one-time message) when none was
## given, otherwise prior itself, which must be one of the four names.
check_prior <- function(prior, defaulted = FALSE) {
  if (defaulted) {
    if (is.null(default_prior_announced$done)) {
      message("Using the default prior joint_partial_r (each partial correlation flat; ",
              "R-squared ~ Beta(k/2, 1)). Set prior = \"independent_partial_r\", ",
              "\"flat_r\" or \"flat_r_squared\" for the others.")
      default_prior_announced$done <- TRUE
    }
    return(default_prior)
  }
  if (is.null(prior))
    stop("Choose a prior: prior = \"", paste(prior_names, collapse = "\", \""), "\".",
         call. = FALSE)
  if (length(prior) != 1 || !prior %in% prior_names)
    stop("Unknown prior \"", paste(prior, collapse = ", "), "\". Choose one of \"",
         paste(prior_names, collapse = "\", \""), "\".", call. = FALSE)
  prior
}

## The shape a of R-squared ~ Beta(a, 1) for k predictor columns. Only the
## three R-based priors have one; independent_partial_r has no such shape.
prior_shape <- function(prior, number_of_predictors) {
  if (prior_is_independent(prior))
    stop("The prior independent_partial_r has no R-squared ~ Beta(a, 1) shape.", call. = FALSE)
  switch(
    check_prior(prior),
    joint_partial_r = number_of_predictors / 2,
    flat_r         = 1 / 2,
    flat_r_squared = 1
  )
}
