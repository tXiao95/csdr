#' Causal sufficient dimension reduction
#'
#' The package entry point is [csdr()], which constructs causal SDR targets and
#' fits MAVE for one or more variants. Advanced users can call [csdr_target()]
#' for target construction only, [estimate_ers()] for exposure-response surface
#' estimation, [estimate_pseudo_outcomes()] for pseudo-outcome construction,
#' and [estimate_residualized_pair()] for residualized-pair construction.
#'
#' Learners are specified with [csdr_learners()], [sl_regression()],
#' [mvn_gps()], [custom_regression()], and [custom_gps()].
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
