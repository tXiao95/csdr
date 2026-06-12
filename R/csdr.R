#' Causal Sufficient Dimension Reduction
#'
#' @param Y Numeric vector of outcomes (length n).
#' @param X Numeric matrix or data frame of observed treatments (n x p).
#' @param C Numeric matrix or data frame of observed confounders (n x q).
#' @param method String or vector indicating which causal transformation(s) to use ("RA", "DR", "PO", "RP").
#' @param args_compute_new_response List of top-level arguments for `compute_new_response_and_exposure`
#'        (e.g., list(L = 10, outcome_fitter = SL_outcome_fitter, seed = 123)).
#' @param args_outcome List of tuning parameters for the outcome fitter.
#' @param args_C List of tuning parameters for the C fitter.
#' @param args_gps List of tuning parameters for the GPS fitter.
#' @param args_ers List of tuning parameters for the ERS estimation.
#' @param args_MAVE List of arguments to pass to `MAVE::mave` (e.g., list(method = "meanOPG")).
#' @return A list containing the MAVE fit, dimension estimation object, selected dimension,
#'         the generated pseudo-data, and run metadata. If multiple methods are requested,
#'         returns a named list of such objects.
#'
#' @export

csdr <- function(Y, X, C,
                   method = c("DR"), # Defaulting to DR, but allowing multiple
                   args_compute_new_response = list(),
                   args_outcome = list(),
                   args_C = list(),
                   args_gps = list(),
                   args_ers = list(),
                   args_MAVE = list()) {

  valid_methods <- c("RA", "DR", "PO", "RP")
  method <- match.arg(method, choices = valid_methods, several.ok = TRUE)

  # ---------------------------------------------------------
  # 1. Compute the New Response and Exposure (Causal Transformation)
  # ---------------------------------------------------------

  # Build the base argument list mapping to the exact arguments required
  cre_base_args <- list(
    Y = Y,
    X = X,
    C = C,
    method = method,
    args_outcome = args_outcome,
    args_C = args_C,
    args_gps = args_gps,
    args_ers = args_ers
  )

  # Merge with any user-supplied top-level args (like L, seed, fitters).
  cre_final_args <- utils::modifyList(cre_base_args, args_compute_new_response)

  # Execute the cross-fitting pipeline
  message("Computing new response and exposure...")
  new_data_obj <- do.call(compute_new_response_and_exposure, cre_final_args)

  # ---------------------------------------------------------
  # 2. Internal Function to Fit MAVE per Method
  # ---------------------------------------------------------

  fit_mave_for_method <- function(m_name, y_vec, x_mat) {
    df <- data.frame(newY = y_vec, x_mat)

    mave_base_args <- list(formula = newY ~ ., data = df, method = "meanMAVE")
    mave_final_args <- utils::modifyList(mave_base_args, args_MAVE)

    fit_mave <- do.call(MAVE::mave, mave_final_args)
    dhat_obj <- MAVE::mave.dim(fit_mave, max.dim = ncol(x_mat))
    d_hat <- if (!is.null(dhat_obj$dim)) dhat_obj$dim.min else NA

    list(
      mave_fit     = fit_mave,
      mave_dim_obj = dhat_obj,
      d_hat        = d_hat,
      new_data     = list(new_Y = y_vec, new_X = x_mat),
      metadata     = list(
        causal_method   = m_name,
        mave_method     = mave_final_args$method,
        n_observations  = nrow(x_mat),
        p_exposures     = ncol(x_mat),
        cre_pipeline    = new_data_obj$metadata
      )
    )
  }

  # ---------------------------------------------------------
  # 3. Fit MAVE for all requested methods
  # ---------------------------------------------------------
  message("Running MAVE...")

  if (length(method) == 1) {
    # Single method: Return a flat, standard object for backwards compatibility
    res <- fit_mave_for_method(method, new_data_obj$new_Y, new_data_obj$new_X)
  } else {
    # Multiple methods: Loop over the lists returned by compute_new_response_and_exposure
    res <- lapply(method, function(m) {
      fit_mave_for_method(m, new_data_obj$new_Y[[m]], new_data_obj$new_X[[m]])
    })
    names(res) <- method
  }

  message("DONE!")
  return(res)
}
