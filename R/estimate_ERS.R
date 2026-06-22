#' Estimate Causal Mean using RA, IPW, or DR (Self-Normalized) from Colangelo and Lee (2026)
#'
#' @param Y Numeric vector of observed outcomes (length n).
#' @param X Numeric matrix or data frame of observed treatments (n x p).
#' @param C Numeric matrix or data frame of observed confounders (n x q).
#' @param x_eval Matrix or data frame of target treatment values to evaluate (m x p).
#' @param estimator String indicating the estimator to use: "DR", "RA", or "IPW". Defaults to "DR".
#' @param out_model An S3 object of class "outcome_model" for the outcome regression. Required for RA and DR.
#' @param gps_model An S3 object of class "gps_model" for the GPS. Required for IPW and DR.
#' @param h Optional numeric vector of bandwidths (length p). Defaults to rule-of-thumb.
#' @param c_multiplier Numeric scalar for rule-of-thumb bandwidth calculation. Defaults to 1.25.
#' @param delta_n Small positive threshold to floor propensity scores. Defaults to 1e-4.
#' @param optimize_bw Logical. If TRUE and estimator is "DR", dynamically calculates the AMSE-optimal
#'        bandwidth using pilot estimates, undersmooths it, and re-estimates. Defaults to FALSE.
#' @param return_vector Logical. If TRUE, forces DR to return only a numeric vector of point estimates.
#'        RA and IPW always return a numeric vector regardless of this setting. Defaults to FALSE.
#'
#' @return A list containing a 'results' data frame (estimates and CIs) and 'metadata' (diagnostics),
#'         OR a numeric vector of estimates if return_vector is TRUE or estimator is RA/IPW.
#'
#' @export

estimate_ERS <- function(Y, X, C, x_eval = NULL,
                             estimator = c("DR", "RA", "IPW"),
                             out_model = NULL,
                             gps_model = NULL,
                             h = NULL,
                             c_multiplier = 1.25,
                             delta_n = 1e-16,
                             optimize_bw = FALSE,
                             return_vector = FALSE) {

  estimator <- match.arg(estimator)

  # Capture original names BEFORE coercion
  orig_X_names <- if(is.matrix(X) || is.data.frame(X)) colnames(X) else NULL
  orig_C_names <- if(is.matrix(C) || is.data.frame(C)) colnames(C) else NULL

  X_df <- as.data.frame(X)
  C_df <- as.data.frame(C)

  n <- length(Y)
  p <- ncol(X_df)

  # ---------------------------------------------------------
  # Input Validation & Safe Column Naming
  # ---------------------------------------------------------
  if (estimator %in% c("RA", "DR")) {
    if (is.null(out_model)) stop("An 'out_model' object is required for the RA and DR estimators.")

    # Safely apply names only if the ORIGINAL input lacked them
    if (is.null(orig_X_names) && !is.null(out_model$X_names)) colnames(X_df) <- out_model$X_names
    if (is.null(orig_C_names) && !is.null(out_model$C_names)) colnames(C_df) <- out_model$C_names
  }

  if (estimator %in% c("IPW", "DR")) {
    if (is.null(gps_model)) stop("A 'gps_model' object is required for the IPW and DR estimators.")

    if (is.null(orig_X_names) && !is.null(gps_model$X_names)) colnames(X_df) <- gps_model$X_names
    if (is.null(orig_C_names) && !is.null(gps_model$C_names)) colnames(C_df) <- gps_model$C_names
  }

  if (optimize_bw && estimator != "DR") {
    warning("Bandwidth optimization is only theoretically implemented for the DR estimator. optimize_bw set to FALSE.")
    optimize_bw <- FALSE
  }

  # ---------------------------------------------------------
  # x_eval Setup
  # ---------------------------------------------------------
  if(is.null(x_eval)){
    x_eval_df <- X_df
  } else {
    orig_x_eval_names <- if(is.matrix(x_eval) || is.data.frame(x_eval)) colnames(x_eval) else NULL
    x_eval_df <- as.data.frame(x_eval)
    if (ncol(x_eval_df) != p) stop("x_eval must have the same number of columns as the training X.")

    # Inherit names from X_df if x_eval lacks them
    if (is.null(orig_x_eval_names)) colnames(x_eval_df) <- colnames(X_df)
  }
  m <- nrow(x_eval_df)

  # ---------------------------------------------------------
  # Pre-computations (Outside the Evaluation Loop)
  # ---------------------------------------------------------
  if (estimator %in% c("IPW", "DR")) {
    # 1. Bandwidth Setup
    if (is.null(h)) {
      active_c <- ifelse(optimize_bw, 3.0, c_multiplier)
      h_pilot <- active_c * apply(X_df, 2, stats::sd) * (n^(-0.2))
    } else {
      h_pilot <- h
    }
    if (length(h_pilot) == 1 && p > 1) h_pilot <- rep(h_pilot, p)

    # 2. Predict Propensity Score
    df_observed <- cbind(X_df, C_df)
    pi_hat <- stats::predict(gps_model, newdata = df_observed)
    pi_hat <- pmax(pi_hat, delta_n)

    # 3. Cache inverse values for speed
    inv_pi <- 1 / pi_hat
    ipw_Y_weighted <- Y * inv_pi
  }

  h_names <- paste0("h_used_", 1:p)

  # ---------------------------------------------------------
  # Zero-Copy Evaluation Grid Setup (NEW)
  # ---------------------------------------------------------
  if (estimator %in% c("RA", "DR")) {

    # 1. Base the grid entirely on C
    dt_grid <- data.table::as.data.table(C_df)

    # 2. Pre-allocate the X columns with dummy values
    x_names <- colnames(X_df)
    for (x_col in x_names) {
      data.table::set(dt_grid, j = x_col, value = 0.0)
    }

    # 3. Force the exact column order the models expect
    data.table::setcolorder(dt_grid, c(x_names, colnames(C_df)))
  }

  # ---------------------------------------------------------
  # Main Loop over Target Treatment Values
  # ---------------------------------------------------------
  loop_out <- lapply(1:m, function(i) {
    message("Evaluating ERS row ", i)
    x_target <- as.numeric(x_eval_df[i, ])

    # Prepare dynamic bandwidth output vector
    h_out <- rep(NA_real_, p)
    names(h_out) <- h_names

    # --- Regression Adjustment Component (Hyper-Optimized) ---
    if (estimator %in% c("RA", "DR")) {

      # ZERO-COPY IN-PLACE UPDATE:
      # Omitting the 'i' argument in set() applies the update to ALL n rows instantly.
      for (j_col in seq_along(x_names)) {
        data.table::set(dt_grid, j = x_names[j_col], value = x_target[j_col])
      }

      m_hat <- stats::predict(out_model, newdata = dt_grid)
      ra_est <- mean(m_hat)
    }

    # ==========================================
    # ESTIMATOR: RA
    # ==========================================
    if (estimator == "RA") {
      return(list(
        row = c(estimate = ra_est, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_),
        h = h_out
      ))
    }

    # ==========================================
    # ESTIMATOR: IPW
    # ==========================================
    else if (estimator == "IPW") {
      K_weights <- rep(1, n)
      for (j in 1:p) {
        K_weights <- K_weights * stats::dnorm((X_df[, j] - x_target[j]) / h_pilot[j])
      }
      den <- sum(K_weights * inv_pi)

      h_out[] <- h_pilot
      if (den < 1e-12) {
        return(list(row = c(estimate = NA_real_, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_), h = h_out))
      }

      est <- sum(K_weights * ipw_Y_weighted) / den
      return(list(row = c(estimate = est, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_), h = h_out))
    }

    # ==========================================
    # ESTIMATOR: DR (Self-Normalized)
    # ==========================================
    else if (estimator == "DR") {

      # Abstracted DR math to easily run multiple bandwidth passes
      compute_dr <- function(h_vec) {
        K_weights <- rep(1, n)
        for (j in 1:p) K_weights <- K_weights * stats::dnorm((X_df[, j] - x_target[j]) / h_vec[j])
        den <- sum(K_weights * inv_pi)

        # Safety fallback if evaluated target is entirely outside observed support
        if (den < 1e-12) {
          psi_i <- m_hat - ra_est
          return(list(est = ra_est, se = stats::sd(psi_i)/sqrt(n), var_psi = stats::var(psi_i)))
        }

        w_i <- n * (K_weights * inv_pi) / den
        delta_est <- sum(w_i * (Y - m_hat)) / n
        est <- ra_est + delta_est

        # Centered IF as derived from Zhang & Chen (2025)
        psi_i <- w_i * (Y - m_hat) + (m_hat - est)
        return(list(est = est, se = stats::sd(psi_i)/sqrt(n), var_psi = stats::var(psi_i)))
      }

      if (optimize_bw) {
        # Step 1: Pilot Variance
        res_pilot <- compute_dr(h_pilot)
        V_t <- prod(h_pilot) * res_pilot$var_psi

        # Step 2: Leading Bias (using b = 2h and a = 0.5)
        b_vec <- 2 * h_pilot
        res_b <- compute_dr(b_vec)
        B_t <- (res_b$est - res_pilot$est) / ( (b_vec[1]^2) * 0.75 )

        # Step 3: Compute AMSE Optimal Bandwidth
        if (abs(B_t) < 1e-12) {
          h_opt_scalar <- h_pilot[1]
        } else {
          h_opt_scalar <- ( (p * V_t) / (4 * B_t^2) )^(1/(p+4)) * n^(-1/(p+4))
        }
        h_opt_vec <- h_opt_scalar * (h_pilot / h_pilot[1])

        # Step 4: Undersmoothing & Final Run
        h_final <- 0.8 * h_opt_vec
        h_out[] <- h_final

        final <- compute_dr(h_final)
        return(list(
          row = c(estimate = final$est, se = final$se, ci_lower = final$est - 1.96*final$se, ci_upper = final$est + 1.96*final$se),
          h = h_out
        ))

      } else {
        # Standard execution (No Optimization)
        h_out[] <- h_pilot
        final <- compute_dr(h_pilot)
        return(list(
          row = c(estimate = final$est, se = final$se, ci_lower = final$est - 1.96*final$se, ci_upper = final$est + 1.96*final$se),
          h = h_out
        ))
      }
    }
  })

  # ---------------------------------------------------------
  # Format and Return Final Output
  # ---------------------------------------------------------
  res_matrix <- do.call(rbind, lapply(loop_out, `[[`, "row"))

  # EXPLICIT BYPASS: Return simple vector for RA, IPW, or if return_vector is TRUE
  if (estimator %in% c("RA", "IPW") || return_vector) {
    return(as.numeric(res_matrix[, "estimate"]))
  }

  # Otherwise, assemble and return the full detailed DR object
  bw_matrix  <- do.call(rbind, lapply(loop_out, `[[`, "h"))

  results_df <- cbind(x_eval_df, as.data.frame(res_matrix), as.data.frame(bw_matrix))
  rownames(results_df) <- NULL

  return(list(
    results = results_df,
    metadata = list(
      estimator = estimator,
      n = n,
      optimized_bw = optimize_bw
    )
  ))
}
