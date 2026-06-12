#' Compute Doubly Robust Pseudo-Outcomes from Kennedy et al. (2017) JRSS-B
#'
#' @param Y Numeric vector of observed outcomes (length n).
#' @param X Numeric matrix or data frame of observed treatments (n x p).
#' @param C Numeric matrix or data frame of observed confounders (n x q).
#' @param out_model An S3 object of class "outcome_model".
#' @param gps_model An S3 object representing a global conditional density model.
#' @param delta_n Small positive threshold to floor propensity scores. Defaults to 1e-16.
#' @return A numeric vector of pseudo-outcomes (length n).

estimate_pseudo_outcomes <- function(Y, X, C, out_model, gps_model, delta_n = 1e-16) {

  # ---------------------------------------------------------
  # Input Validation & Safe Column Naming
  # ---------------------------------------------------------
  orig_X_names <- if(is.matrix(X) || is.data.frame(X)) colnames(X) else NULL
  orig_C_names <- if(is.matrix(C) || is.data.frame(C)) colnames(C) else NULL

  X_df <- as.data.frame(X)
  C_df <- as.data.frame(C)
  n <- length(Y)

  if (is.null(orig_X_names) && !is.null(out_model$X_names)) colnames(X_df) <- out_model$X_names
  if (is.null(orig_C_names) && !is.null(out_model$C_names)) colnames(C_df) <- out_model$C_names

  # ---------------------------------------------------------
  # Pre-computations (Batch predicting observed data)
  # ---------------------------------------------------------
  df_observed <- cbind(X_df, C_df)

  m_obs  <- predict(out_model, newdata = df_observed)
  pi_obs <- predict(gps_model, newdata = df_observed)
  pi_obs <- pmax(pi_obs, delta_n)

  # ---------------------------------------------------------
  # Zero-Copy Evaluation Grid Setup
  # ---------------------------------------------------------

  # 1. Base the grid entirely on C (which never changes for the marginalization)
  dt_grid <- as.data.table(C_df)

  # 2. Pre-allocate the X columns with dummy values (0.0)
  x_names <- colnames(X_df)
  for (x_col in x_names) {
    dt_grid[, (x_col) := 0.0]
  }

  # 3. Force the exact column order the models expect
  setcolorder(dt_grid, c(x_names, colnames(C_df)))

  # ---------------------------------------------------------
  # Main Loop over Individuals (Hyper-Optimized)
  # ---------------------------------------------------------
  pseudo_outcomes <- vapply(1:n, function(j) {
    message("Evaluating pseudo-outcome row ", j)

    # ZERO-COPY IN-PLACE UPDATE:
    # Instead of `rep()` and `cbind()`, we instantly overwrite the memory
    # of the X columns with the values for individual 'j'.
    for (x_col in x_names) {
      data.table::set(dt_grid, j = x_col, value = X_df[[x_col]][j])
    }

    # Predict
    m_grid  <- predict(out_model, newdata = dt_grid)
    pi_grid <- predict(gps_model, newdata = dt_grid)

    # Calculate the empirical expectations
    mean_pi <- mean(pi_grid)
    mean_m  <- mean(m_grid)

    # Assemble the pseudo-outcome
    xi_j <- ((Y[j] - m_obs[j]) / pi_obs[j]) * mean_pi + mean_m

    return(xi_j)

  }, numeric(1L))

  return(pseudo_outcomes)
}
