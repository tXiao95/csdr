#' Estimate Cross-Fitted Causal Exposure Response Surface (Colangelo and Lee, 2026)
#'
#' @param Y Numeric vector of observed outcomes (length n).
#' @param X Numeric matrix or data frame of observed treatments (n x p).
#' @param C Numeric matrix or data frame of observed confounders (n x q).
#' @param x_eval Matrix or data frame of target treatment values to evaluate (m x p).
#' @param estimator String indicating the estimator to use: "DR", "RA", or "IPW". Defaults to "DR".
#' @param L Integer indicating the number of folds for cross-fitting. Defaults to 5.
#' @param outcome_fitter Function to train outcome models. Required for RA and DR.
#' @param gps_fitter Function to train GPS models. Required for IPW and DR.
#' @param seed Random seed for cross-fitting folds.
#' @param h Optional numeric vector of bandwidths (length p). Defaults to rule-of-thumb.
#' @param c_multiplier Numeric scalar for rule-of-thumb bandwidth calculation. Defaults to 1.25.
#' @param delta_n Small positive threshold to floor propensity scores. Defaults to 1e-16.
#' @param optimize_bw Logical. If TRUE and estimator is "DR", dynamically calculates the AMSE-optimal bandwidth.
#' @param args_outcome List of additional arguments to pass to outcome_fitter.
#' @param args_gps List of additional arguments to pass to gps_fitter.
#'
#' @return A list containing the globally aggregated 'results' data frame and 'metadata'.
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

  n <- length(Y)
  X_df <- as.data.frame(X)
  C_df <- as.data.frame(C)
  p <- ncol(X_df)

  if (is.null(x_eval)) x_eval <- X_df
  x_eval_df <- as.data.frame(x_eval)
  m <- nrow(x_eval_df)

  if (L == 1) warning("L=1 disables cross-fitting: nuisances and estimation use the same data.")

  set.seed(seed)
  folds <- sample(rep(1:L, length.out = n))

  # ---------------------------------------------------------
  # 1. Pre-allocate Global Prediction Structures
  # ---------------------------------------------------------
  pi_hat_full  <- rep(NA_real_, n)
  m_obs_full   <- rep(NA_real_, n)
  m_target_mat <- matrix(NA_real_, nrow = n, ncol = m)

  # ---------------------------------------------------------
  # 2. Cross-Fitting Loop (Prediction Only)
  # ---------------------------------------------------------
  for (k in 1:L) {
    if (L > 1) message(sprintf("Cross-fitting fold %d of %d...", k, L))

    train_idx <- if (L == 1L) seq_len(n) else which(folds != k)
    test_idx  <- if (L == 1L) seq_len(n) else which(folds == k)

    # --- Train GPS and Predict ---
    if (estimator %in% c("IPW", "DR")) {
      gps_req <- list(X = X_df[train_idx, , drop=FALSE], C = C_df[train_idx, , drop=FALSE], pi_fitter = gps_fitter)
      gps_mod <- do.call(gps_model, c(gps_req, args_gps))

      pi_hat_full[test_idx] <- stats::predict(gps_mod, newdata = cbind(X_df[test_idx, , drop=FALSE], C_df[test_idx, , drop=FALSE]))
    }

    # --- Train Outcome Model and Predict ---
    if (estimator %in% c("RA", "DR")) {
      out_req <- list(Y = Y[train_idx], X = X_df[train_idx, , drop=FALSE], C = C_df[train_idx, , drop=FALSE], mu_fitter = outcome_fitter)
      out_mod <- do.call(outcome_model, c(out_req, args_outcome))

      m_obs_full[test_idx] <- stats::predict(out_mod, newdata = cbind(X_df[test_idx, , drop=FALSE], C_df[test_idx, , drop=FALSE]))

      # Zero-Copy Target Evaluation for the test fold
      dt_grid <- data.table::as.data.table(C_df[test_idx, , drop=FALSE])
      x_names <- colnames(X_df)
      for (x_col in x_names) {
        data.table::set(dt_grid, j = x_col, value = 0.0)
      }
      data.table::setcolorder(dt_grid, c(x_names, colnames(C_df)))

      for (j in 1:m) {
        x_target <- as.numeric(x_eval_df[j, ])
        for (j_col in seq_along(x_names)) {
          data.table::set(dt_grid, j = x_names[j_col], value = x_target[j_col])
        }
          m_target_mat[test_idx, j] <- stats::predict(out_mod, newdata = dt_grid)
      }
    }
    # Clearing memory for nnet
    gc(verbose = FALSE)
  }

  # ---------------------------------------------------------
  # 3. Global Kernel Smoothing & Bandwidth Optimization
  # ---------------------------------------------------------
  message("Executing global causal integration...")

  if (estimator %in% c("IPW", "DR")) {
    pi_hat_full <- pmax(pi_hat_full, delta_n)
    inv_pi <- 1 / pi_hat_full
    ipw_Y_weighted <- Y * inv_pi

    if (is.null(h)) {
      active_c <- ifelse(optimize_bw, 3.0, c_multiplier)
      h_pilot <- active_c * apply(X_df, 2, stats::sd) * (n^(-0.2)) # Correctly scales by global N
    } else {
      h_pilot <- h
    }
    if (length(h_pilot) == 1 && p > 1) h_pilot <- rep(h_pilot, p)
  }

  h_names <- paste0("h_used_", 1:p)

  loop_out <- lapply(1:m, function(j) {
    x_target <- as.numeric(x_eval_df[j, ])
    h_out <- rep(NA_real_, p)
    names(h_out) <- h_names

    if (estimator %in% c("RA", "DR")) {
      ra_est <- mean(m_target_mat[, j])
    }

    if (estimator == "RA") {
      return(list(row = c(estimate = ra_est, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_), h = h_out))
    }

    else if (estimator == "IPW") {
      K_weights <- rep(1, n)
        for (dim_idx in 1:p) K_weights <- K_weights * stats::dnorm((X_df[, dim_idx] - x_target[dim_idx]) / h_pilot[dim_idx])
      den <- sum(K_weights * inv_pi)
      h_out[] <- h_pilot

      if (den < 1e-12) return(list(row = c(estimate = NA_real_, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_), h = h_out))

      est <- sum(K_weights * ipw_Y_weighted) / den
      return(list(row = c(estimate = est, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_), h = h_out))
    }

    else if (estimator == "DR") {
      m_hat_j <- m_target_mat[, j]

      compute_dr <- function(h_vec) {
          K_weights <- rep(1, n)
          for (dim_idx in 1:p) K_weights <- K_weights * stats::dnorm((X_df[, dim_idx] - x_target[dim_idx]) / h_vec[dim_idx])
          den <- sum(K_weights * inv_pi)

        if (den < 1e-12) {
          psi_i <- m_hat_j - ra_est
          return(list(est = ra_est, se = stats::sd(psi_i)/sqrt(n), var_psi = stats::var(psi_i)))
        }

        w_i <- n * (K_weights * inv_pi) / den
        delta_est <- sum(w_i * (Y - m_obs_full)) / n
        est <- ra_est + delta_est

        # The true DML influence function utilizing the full N cross-fitted residuals
        psi_i <- w_i * (Y - m_obs_full) + (m_hat_j - est)
        return(list(est = est, se = stats::sd(psi_i)/sqrt(n), var_psi = stats::var(psi_i)))
      }

      if (optimize_bw) {
        res_pilot <- compute_dr(h_pilot)
        V_t <- prod(h_pilot) * res_pilot$var_psi

        b_vec <- 2 * h_pilot
        res_b <- compute_dr(b_vec)
        B_t <- (res_b$est - res_pilot$est) / ( (b_vec[1]^2) * 0.75 )

        if (abs(B_t) < 1e-12) {
          h_opt_scalar <- h_pilot[1]
        } else {
          h_opt_scalar <- ( (p * V_t) / (4 * B_t^2) )^(1/(p+4)) * n^(-1/(p+4))
        }
        h_opt_vec <- h_opt_scalar * (h_pilot / h_pilot[1])

        h_final <- 0.8 * h_opt_vec
        h_out[] <- h_final
        final <- compute_dr(h_final)

      } else {
        h_out[] <- h_pilot
        final <- compute_dr(h_pilot)
      }

      return(list(
        row = c(estimate = final$est, se = final$se, ci_lower = final$est - 1.96*final$se, ci_upper = final$est + 1.96*final$se),
        h = h_out
      ))
    }
  })

  # ---------------------------------------------------------
  # 4. Final Output Formatting
  # ---------------------------------------------------------
  res_matrix <- do.call(rbind, lapply(loop_out, `[[`, "row"))
  bw_matrix  <- do.call(rbind, lapply(loop_out, `[[`, "h"))

  results_df <- cbind(x_eval_df, as.data.frame(res_matrix), as.data.frame(bw_matrix))
  rownames(results_df) <- NULL

  return(list(
    results = results_df,
    metadata = list(
      estimator = estimator,
      n_total = n,
      L_folds = L,
      seed = seed,
      optimized_bw = optimize_bw
    )
  ))
}
