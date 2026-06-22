#' Compute New Response and Exposure for Causal SDR via Cross-Fitting
#'
#' @param Y Numeric vector of outcomes (length n).
#' @param X Numeric matrix or data frame of observed treatments (n x p).
#' @param C Numeric matrix or data frame of observed confounders (n x q).
#' @param method String or vector indicating methods to use (e.g., "DR" or c("RA", "DR", "PO", "RP")).
#' @param L Integer indicating the number of folds for cross-fitting. Defaults to 5.
#' @param outcome_fitter Function to train outcome models. Required for RA, DR, PO.
#' @param C_fitter Function to train nuisance models. Required for RP.
#' @param gps_fitter Function to train GPS models. Defaults to mvn_fitter. Required for DR, PO.
#' @param seed Random seed for cross-fitting folds.
#' @param args_outcome List of additional arguments to pass to outcome_fitter.
#' @param args_C List of additional arguments to pass to C_fitter.
#' @param args_gps List of additional arguments to pass to gps_fitter (e.g., list(method="linear")).
#' @param args_ers List of additional arguments to pass to estimate_ERS (e.g., list(optimize_bw=TRUE)).
#' @return A list containing `new_Y`, `new_X`, and `metadata`.

compute_new_response_and_exposure <- function(Y, X, C, 
                                              method = c("RA"), # Default to DR
                                              L = 5,
                                              outcome_fitter = SL_outcome_fitter,
                                              C_fitter = SL_nuisance_fitter,
                                              gps_fitter = mvn_fitter,
                                              seed = 42,
                                              args_outcome = list(),
                                              args_C = list(),
                                              args_gps = list(),
                                              args_ers = list()) {
  
  # Allow multiple methods to be computed in a single cross-fitting pass
  valid_methods <- c("RA", "DR", "PO", "RP")
  method <- match.arg(method, choices = valid_methods, several.ok = TRUE)
  
  # ---------------------------------------------------------
  # 1. Input Validation and Setup
  # ---------------------------------------------------------
  n <- length(Y)
  X_mat <- as.matrix(X)
  C_mat <- as.matrix(C)
  p <- ncol(X_mat)
  
  if (is.null(colnames(X_mat))) colnames(X_mat) <- paste0("X", 1:p)
  
  if (any(c("RA", "DR", "PO") %in% method) && !is.function(outcome_fitter)) {
    stop("A valid 'outcome_fitter' must be provided for RA, DR, or PO.")
  }
  if ("RP" %in% method && !is.function(C_fitter)) {
    stop("A valid 'C_fitter' must be provided for RP.")
  }
  if (any(c("DR", "PO") %in% method) && !is.function(gps_fitter)) {
    stop("A valid 'gps_fitter' must be provided for DR or PO.")
  }
  if( L == 1) warning("L=1 disables cross-fitting.")
  po_gps_floor <- 1e-8
  if (!is.null(args_ers$gps_floor)) {
    po_gps_floor <- args_ers$gps_floor
  } else if (!is.null(args_ers$delta_n)) {
    po_gps_floor <- args_ers$delta_n
  }
  
  # ---------------------------------------------------------
  # 2. Pre-allocate Output Structures for ALL requested methods
  # ---------------------------------------------------------
  new_Y_list <- lapply(method, function(m) rep(NA_real_, n))
  names(new_Y_list) <- method

  new_X_list <- lapply(method, function(m) X_mat)
  names(new_X_list) <- method
  po_m_obs <- po_pi_obs <- rep(NA_real_, n)
  po_nuisance_folds <- if ("PO" %in% method) vector("list", L) else NULL
  
  set.seed(seed)
  folds <- sample(rep(1:L, length.out = n))
  
  # ---------------------------------------------------------
  # 3. Main Cross-Fitting Loop
  # ---------------------------------------------------------
  for (k in 1:L) {
    train_idx <- if (L == 1L) seq_len(n) else which(folds != k)
    test_idx  <- if (L == 1L) seq_len(n) else which(folds == k)
    
    Y_train <- Y[train_idx]; X_train <- X_mat[train_idx, , drop = FALSE]; C_train <- C_mat[train_idx, , drop = FALSE]
    Y_test  <- Y[test_idx];  X_test  <- X_mat[test_idx, , drop = FALSE];  C_test  <- C_mat[test_idx, , drop = FALSE]
    
    # --- Step A: Estimate Nuisances (Only once per requirement) ---
    if (any(c("RA", "DR", "PO") %in% method)) {
      out_req_args <- list(Y = Y_train, X = X_train, C = C_train, mu_fitter = outcome_fitter)
      out_mod <- do.call(outcome_model, c(out_req_args, args_outcome))
    }
    if (any(c("DR", "PO") %in% method)) {
      gps_req_args <- list(X = X_train, C = C_train, pi_fitter = gps_fitter)
      gps_mod <- do.call(gps_model, c(gps_req_args, args_gps))
    }
    if ("RP" %in% method) {
      C_mods <- fit_rp_nuisance(
        Y = Y_train,
        A = X_train,
        C = C_train,
        y_fitter = C_fitter,
        a_fitter = C_fitter,
        args_y = args_C,
        args_a = args_C,
        verbose = FALSE
      )
    }
    
    # --- Step B: Generate Pseudo-Data for each requested method ---
    if ("RA" %in% method) {
      new_Y_list[["RA"]][test_idx] <- do.call(estimate_ERS, c(list(Y=Y_test, X=X_test, C=C_test, estimator="RA", out_model=out_mod, return_vector=TRUE), args_ers))
    }
    if ("DR" %in% method) {
      new_Y_list[["DR"]][test_idx] <- do.call(estimate_ERS, c(list(Y=Y_test, X=X_test, C=C_test, estimator="DR", out_model=out_mod, gps_model=gps_mod, return_vector=TRUE), args_ers))
    }
    if ("PO" %in% method) {
      po_nuisance <- list(
        outcome_model = out_mod,
        gps_model = gps_mod,
        A_names = out_mod$X_names,
        C_names = out_mod$C_names
      )
      class(po_nuisance) <- "po_nuisance"
      po_nuisance_folds[[k]] <- po_nuisance
      po_pred <- predict_po_observed(
        nuisance = po_nuisance,
        A = X_test,
        C = C_test,
        gps_floor = po_gps_floor
      )
      po_m_obs[test_idx] <- po_pred$m_obs
      po_pi_obs[test_idx] <- po_pred$pi_obs
    }
    if ("RP" %in% method) {
      rp_pred <- predict_rp_nuisance(nuisance = C_mods, C = C_test)
      res_rp <- compute_residualized_pair(
        Y = Y_test,
        A = X_test,
        EY_C = rp_pred$EY_C,
        EA_C = rp_pred$EA_C
      )
      new_Y_list[["RP"]][test_idx] <- res_rp$Y_tilde
      new_X_list[["RP"]][test_idx, ] <- res_rp$A_tilde
    }
  }

  if ("PO" %in% method) {
    po_marginal <- predict_po_crossfit_marginals(
      nuisance_fits = list(folds = po_nuisance_folds),
      A_targets = X_mat,
      C = C_mat,
      folds = folds,
      gps_floor = po_gps_floor,
      verbose = FALSE
    )
    new_Y_list[["PO"]] <- compute_pseudo_outcomes(
      Y = Y,
      m_obs = po_m_obs,
      pi_obs = po_pi_obs,
      m_marginal = po_marginal$m_marginal,
      pi_marginal = po_marginal$pi_marginal,
      gps_floor = po_gps_floor
    )
  }
  
  # ---------------------------------------------------------
  # 4. Final Output Formatting
  # ---------------------------------------------------------
  meta <- list(method = method, L_folds = L, n = n, p = p, seed = seed)
  
  # Backwards compatibility: If only 1 method was requested, return flat lists
  if (length(method) == 1) {
    return(list(new_Y = new_Y_list[[1]], new_X = as.matrix(new_X_list[[1]]), metadata = meta))
  } else {
    return(list(new_Y = new_Y_list, new_X = lapply(new_X_list, as.matrix), metadata = meta))
  }
}
