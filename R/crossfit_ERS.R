#' Estimate a cross-fitted causal exposure-response surface
#'
#' `crossfit_ERS()` is retained as a compatibility wrapper around
#' [estimate_ers()]. New code should call [estimate_ers()] directly.
#'
#' @param Y Numeric vector of observed outcomes.
#' @param X Numeric matrix or data frame of observed exposures/treatments.
#' @param C Numeric matrix or data frame of observed covariates.
#' @param x_eval Optional matrix or data frame of exposure values where the
#'   exposure-response surface should be evaluated.
#' @param estimator Estimator to use: `"DR"`, `"RA"`, or `"IPW"`.
#' @param L Number of folds for cross-fitting. `L = 1` disables cross-fitting.
#' @param outcome_fitter Function used to fit the outcome nuisance.
#' @param gps_fitter Function used to fit the generalized propensity score.
#' @param seed Optional random seed for fold construction.
#' @param h Optional bandwidth vector.
#' @param c_multiplier Multiplier used for the default bandwidth.
#' @param delta_n Positive floor applied to GPS estimates.
#' @param optimize_bw Deprecated; bandwidth optimization is not implemented in
#'   this wrapper.
#' @param args_outcome Additional arguments passed to `outcome_fitter`.
#' @param args_gps Additional arguments passed to `gps_fitter`.
#'
#' @return An object of class `"ers_fit"` with a legacy `metadata` element.
#'
#' @export
crossfit_ERS <- function(Y, X, C, x_eval = NULL,
                         estimator = c("DR", "RA", "IPW"),
                         L = 5,
                         outcome_fitter = SL_outcome_fitter,
                         gps_fitter = mvn_fitter,
                         seed = 42,
                         h = NULL,
                         c_multiplier = 1.25,
                         delta_n = 1e-16,
                         optimize_bw = FALSE,
                         args_outcome = list(),
                         args_gps = list()) {
  estimator <- match.arg(estimator)

  if (isTRUE(optimize_bw)) {
    warning(
      "'optimize_bw' is not implemented in the estimate_ers() wrapper and is ignored.",
      call. = FALSE
    )
  }

  fit <- estimate_ers(
    Y = Y,
    A = X,
    C = C,
    a_eval = x_eval,
    estimator = estimator,
    L = L,
    outcome_fitter = outcome_fitter,
    gps_fitter = gps_fitter,
    h = h,
    c_multiplier = c_multiplier,
    gps_floor = delta_n,
    seed = seed,
    verbose = TRUE,
    args_outcome = args_outcome,
    args_gps = args_gps
  )

  fit$metadata <- list(
    estimator = fit$estimator,
    n_total = length(Y),
    L_folds = fit$L,
    seed = seed,
    optimized_bw = FALSE
  )

  fit
}
