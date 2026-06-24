#' Compute New Response and Exposure for Causal SDR via Cross-Fitting
#'
#' @param Y Numeric vector of outcomes (length n).
#' @param X Numeric matrix or data frame of observed treatments (n x p).
#' @param C Numeric matrix or data frame of observed confounders (n x q).
#' @param method String or vector indicating methods to use (e.g., "DR" or c("RA", "DR", "PO", "RP")).
#' @param L Integer indicating the number of folds for cross-fitting. Defaults to 5.
#' @param folds Optional fold assignments.
#' @param outcome_fitter Function to train outcome models. Required for RA, DR, PO.
#' @param C_fitter Function to train nuisance models. Required for RP.
#' @param gps_fitter Function to train GPS models. Defaults to mvn_fitter. Required for DR, PO.
#' @param seed Random seed for cross-fitting folds.
#' @param args_outcome List of additional arguments to pass to outcome_fitter.
#' @param args_C List of additional arguments to pass to C_fitter.
#' @param args_gps List of additional arguments to pass to gps_fitter (e.g., list(method="linear")).
#' @param args_ers List of additional arguments for ERS target construction,
#'   such as `h`, `c_multiplier`, or `gps_floor`.
#' @param po_marginalization Marginalization mode for PO targets. Defaults to
#'   `"crossfit"` through [csdr_target()].
#' @param verbose Logical; if `TRUE`, print progress messages.
#' @return A list containing `new_Y`, `new_X`, and `metadata`.

compute_new_response_and_exposure <- function(Y, X, C, 
                                              method = c("RA"), # Default to DR
                                              L = 5,
                                              folds = NULL,
                                              outcome_fitter = SL_outcome_fitter,
                                              C_fitter = SL_nuisance_fitter,
                                              gps_fitter = mvn_fitter,
                                              seed = 42,
                                              args_outcome = list(),
                                              args_C = list(),
                                              args_gps = list(),
                                              args_ers = list(),
                                              po_marginalization = NULL,
                                              verbose = TRUE) {
  
  # Allow multiple methods to be computed in a single cross-fitting pass
  valid_methods <- c("RA", "DR", "PO", "RP")
  method <- match.arg(method, choices = valid_methods, several.ok = TRUE)

  target_args_ers <- args_ers
  if (is.null(po_marginalization)) {
    po_marginalization <- if (!is.null(target_args_ers$po_marginalization)) {
      target_args_ers$po_marginalization
    } else if (!is.null(target_args_ers$marginalization)) {
      target_args_ers$marginalization
    } else {
      "crossfit"
    }
  }
  target_args_ers$po_marginalization <- NULL
  target_args_ers$marginalization <- NULL

  if (L == 1L) {
    warning("L=1 disables cross-fitting.")
  }

  target <- csdr_target(
    Y = Y,
    A = X,
    C = C,
    methods = method,
    L = L,
    folds = folds,
    outcome_fitter = outcome_fitter,
    gps_fitter = gps_fitter,
    C_fitter = C_fitter,
    args_outcome = args_outcome,
    args_gps = args_gps,
    args_C = args_C,
    args_ers = target_args_ers,
    po_marginalization = po_marginalization,
    seed = seed,
    return_nuisance = FALSE,
    verbose = verbose
  )
  
  # ---------------------------------------------------------
  # 4. Final Output Formatting
  # ---------------------------------------------------------
  meta <- list(
    method = method,
    L_folds = target$diagnostics$L,
    n = target$diagnostics$n,
    p = target$diagnostics$p,
    seed = seed,
    target_diagnostics = target$diagnostics
  )
  
  # Backwards compatibility: If only 1 method was requested, return flat lists
  if (length(method) == 1) {
    return(list(
      new_Y = target$new_Y[[method]],
      new_X = as.matrix(target$new_A[[method]]),
      metadata = meta
    ))
  } else {
    return(list(
      new_Y = target$new_Y,
      new_X = lapply(target$new_A, as.matrix),
      metadata = meta
    ))
  }
}
