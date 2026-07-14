#' Compute residualized response and exposure
#'
#' @param Y Numeric outcome vector.
#' @param A Numeric matrix or data frame of observed exposures.
#' @param EY_C Numeric vector of predictions for `E[Y | C]`.
#' @param EA_C Numeric matrix of predictions for `E[A | C]`.
#'
#' @return A list with `Y_tilde` and `A_tilde`.
compute_residualized_pair <- function(Y, A, EY_C, EA_C) {
  A_df <- as.data.frame(A)
  if (!is.numeric(Y) || !is.numeric(EY_C) ||
      !all(vapply(A_df, is.numeric, logical(1L))) || !is.numeric(EA_C)) {
    stop("'Y', 'A', 'EY_C', and 'EA_C' must be numeric.", call. = FALSE)
  }
  Y <- as.numeric(Y)
  A_mat <- as.matrix(A_df)
  EY_C <- as.numeric(EY_C)
  EA_C <- as.matrix(EA_C)

  if (length(Y) != length(EY_C)) {
    stop("'Y' and 'EY_C' must have the same length.", call. = FALSE)
  }
  if (nrow(A_mat) != length(Y)) {
    stop("'A' must have the same number of rows as 'Y'.", call. = FALSE)
  }
  if (!identical(dim(EA_C), dim(A_mat))) {
    stop("'EA_C' must have the same dimensions as 'A'.", call. = FALSE)
  }
  if (anyNA(Y) || anyNA(A_mat) || anyNA(EY_C) || anyNA(EA_C)) {
    stop("Residualized-pair inputs must not contain missing values.", call. = FALSE)
  }

  A_tilde <- A_mat - EA_C
  colnames(A_tilde) <- colnames(A_mat)
  list(
    Y_tilde = as.numeric(Y - EY_C),
    A_tilde = as.matrix(A_tilde)
  )
}

#' Fit nuisance models for residualized-pair construction
#'
#' @param Y Numeric outcome vector.
#' @param A Numeric matrix or data frame of observed exposures.
#' @param C Numeric matrix or data frame of covariates.
#' @param y_fitter Function used to fit `E[Y | C]`.
#' @param a_fitter Function used to fit each `E[A_j | C]`.
#' @param args_y Additional arguments passed to `y_fitter`.
#' @param args_a Additional arguments passed to `a_fitter`.
#' @param verbose Logical; if `TRUE`, print progress messages.
#'
#' @return A fitted residualized-pair nuisance object.
fit_rp_nuisance <- function(Y, A, C,
                            y_fitter = SL_nuisance_fitter,
                            a_fitter = SL_nuisance_fitter,
                            args_y = list(),
                            args_a = list(),
                            verbose = TRUE) {
  data <- normalize_ers_inputs(Y = Y, A = A, C = C, a_eval = A)

  if (!is.function(y_fitter)) {
    stop("'y_fitter' must be a function.", call. = FALSE)
  }
  if (!is.function(a_fitter)) {
    stop("'a_fitter' must be a function.", call. = FALSE)
  }

  if (isTRUE(verbose)) {
    message("Fitting residualized-pair nuisance models.")
  }

  Y_model <- do.call(
    nuisance_C_model,
    c(list(target = data$Y, C = data$C, fitter = y_fitter), args_y)
  )

  p <- ncol(data$A)
  A_models <- vector("list", p)
  for (j in seq_len(p)) {
    if (isTRUE(verbose)) {
      message(sprintf("Fitting E[A%d | C].", j))
    }
    A_models[[j]] <- do.call(
      nuisance_C_model,
      c(list(target = data$A[[j]], C = data$C, fitter = a_fitter), args_a)
    )
  }

  out <- list(
    Y_model = Y_model,
    A_models = A_models,
    A_names = colnames(data$A),
    C_names = colnames(data$C)
  )
  class(out) <- "rp_nuisance"
  out
}

#' Predict nuisance quantities for residualized-pair construction
#'
#' @param nuisance A fitted residualized-pair nuisance object.
#' @param C Numeric matrix or data frame of covariates.
#' @param verbose Logical; if `TRUE`, print progress messages.
#'
#' @return A list with `EY_C` and `EA_C`.
predict_rp_nuisance <- function(nuisance, C, verbose = FALSE) {
  if (is.null(nuisance$Y_model)) {
    stop("'nuisance' must contain 'Y_model'.", call. = FALSE)
  }
  A_models <- nuisance$A_models
  if (is.null(A_models) && !is.null(nuisance$X_models)) {
    A_models <- nuisance$X_models
  }
  if (!is.list(A_models) || length(A_models) < 1L) {
    stop("'nuisance' must contain a nonempty 'A_models' list.", call. = FALSE)
  }

  C_df <- normalize_target_table(C, n = NULL, arg = "C", prefix = "C")
  colnames(C_df) <- nuisance$Y_model$C_names

  EY_C <- as.numeric(stats::predict(nuisance$Y_model, newdata = C_df))
  EA_C <- matrix(NA_real_, nrow = nrow(C_df), ncol = length(A_models))
  for (j in seq_along(A_models)) {
    if (isTRUE(verbose)) {
      message(sprintf("Predicting E[A%d | C].", j))
    }
    EA_C[, j] <- as.numeric(stats::predict(A_models[[j]], newdata = C_df))
  }
  colnames(EA_C) <- if (!is.null(nuisance$A_names)) {
    nuisance$A_names
  } else {
    paste0("A", seq_along(A_models))
  }

  if (nrow(EA_C) != nrow(C_df)) {
    stop("'EA_C' must have one row per row of 'C'.", call. = FALSE)
  }

  list(EY_C = EY_C, EA_C = as.matrix(EA_C))
}

#' Estimate residualized response and exposure pairs for CSDR
#'
#' @param Y Numeric outcome vector.
#' @param A Numeric matrix or data frame of observed exposures.
#' @param C Numeric matrix or data frame of covariates.
#' @param L Number of folds. `L = 1` disables cross-fitting.
#' @param folds Optional fold assignments.
#' @param C_fitter Default fitter used for both residualization models.
#' @param y_fitter Function used to fit `E[Y | C]`.
#' @param a_fitter Function used to fit each `E[A_j | C]`.
#' @param nuisance Optional pre-fitted residualized-pair nuisance object.
#' @param seed Optional random seed for fold generation.
#' @param return_nuisance Logical; if `TRUE`, return nuisance fits.
#' @param verbose Logical; if `TRUE`, print progress messages.
#' @param args_y Additional arguments passed to `y_fitter`.
#' @param args_a Additional arguments passed to `a_fitter`.
#' @param X Deprecated compatibility alias for `A`.
#' @param C_models Deprecated compatibility alias for a pre-fitted RP nuisance.
#'
#' @return An S3 object of class `"rp_fit"` and `"csdr_target"`.
#'
#' @export
estimate_residualized_pair <- function(
  Y,
  A = NULL,
  C,
  L = 5,
  folds = NULL,
  C_fitter = SL_nuisance_fitter,
  y_fitter = C_fitter,
  a_fitter = C_fitter,
  nuisance = NULL,
  seed = NULL,
  return_nuisance = FALSE,
  verbose = TRUE,
  args_y = list(),
  args_a = list(),
  X = NULL,
  C_models = NULL
) {
  call <- match.call()
  legacy_return <- is.null(A) && !is.null(X) && !is.null(C_models)
  if (is.null(A)) {
    A <- X
  }
  if (is.null(nuisance) && !is.null(C_models)) {
    nuisance <- as_rp_nuisance(C_models)
  }

  data <- normalize_ers_inputs(Y = Y, A = A, C = C, a_eval = A)
  validate_ers_options(L = L, gps_floor = 1, normalize_ipw = TRUE)
  folds <- make_folds(n = length(data$Y), L = L, folds = folds, seed = seed)
  L_eff <- length(unique(folds))

  nuisance_fits <- prepare_rp_nuisance_fits(
    nuisance = nuisance,
    Y = data$Y,
    A = data$A,
    C = data$C,
    L = L_eff,
    folds = folds,
    y_fitter = y_fitter,
    a_fitter = a_fitter,
    args_y = args_y,
    args_a = args_a,
    verbose = verbose
  )

  n <- length(data$Y)
  p <- ncol(data$A)
  EY_C <- rep(NA_real_, n)
  EA_C <- matrix(NA_real_, nrow = n, ncol = p)
  colnames(EA_C) <- colnames(data$A)

  for (k in seq_len(L_eff)) {
    test_idx <- if (L_eff == 1L) seq_len(n) else which(folds == k)
    pred <- predict_rp_nuisance(
      nuisance = nuisance_fits$folds[[k]],
      C = data$C[test_idx, , drop = FALSE],
      verbose = verbose && L_eff > 1L
    )
    EY_C[test_idx] <- pred$EY_C
    EA_C[test_idx, ] <- pred$EA_C
  }

  residualized <- compute_residualized_pair(
    Y = data$Y,
    A = data$A,
    EY_C = EY_C,
    EA_C = EA_C
  )

  if (isTRUE(legacy_return)) {
    return(list(
      Ytilde = residualized$Y_tilde,
      Xtilde = residualized$A_tilde
    ))
  }

  out <- list(
    Y_tilde = residualized$Y_tilde,
    A_tilde = residualized$A_tilde,
    Ytilde = residualized$Y_tilde,
    Xtilde = residualized$A_tilde,
    folds = folds,
    nuisance = if (return_nuisance) nuisance_fits else NULL,
    call = call,
    diagnostics = list(
      method = "RP",
      L = L_eff,
      requested_L = L,
      n = n,
      p = p
    )
  )
  class(out) <- c("rp_fit", "csdr_target")
  out
}

#' Train all nuisance models for residualized-pair construction
#'
#' @param Y Numeric outcome vector.
#' @param X Deprecated exposure argument name.
#' @param C Numeric matrix or data frame of covariates.
#' @param fitter Fitter used for both outcome and exposure nuisance models.
#' @param ... Additional arguments passed to `fitter`.
#'
#' @return A residualized-pair nuisance object with legacy `X_models` alias.
train_nuisance_models <- function(Y, X, C, fitter, ...) {
  out <- fit_rp_nuisance(
    Y = Y,
    A = X,
    C = C,
    y_fitter = fitter,
    a_fitter = fitter,
    args_y = list(...),
    args_a = list(...),
    verbose = TRUE
  )
  out$X_models <- out$A_models
  out
}

prepare_rp_nuisance_fits <- function(nuisance, Y, A, C, L, folds,
                                     y_fitter, a_fitter,
                                     args_y, args_a, verbose) {
  if (!is.null(nuisance)) {
    nuisance <- as_rp_nuisance(nuisance)
    fold_fits <- if (!is.null(nuisance$folds)) {
      lapply(nuisance$folds, as_rp_nuisance)
    } else if (!is.null(nuisance$fold_nuisance)) {
      lapply(nuisance$fold_nuisance, as_rp_nuisance)
    } else {
      rep(list(nuisance), L)
    }
    if (length(fold_fits) != L) {
      stop("Supplied residualized-pair nuisance fits do not match the number of folds.", call. = FALSE)
    }
    return(list(folds = fold_fits, supplied = TRUE))
  }

  fold_fits <- vector("list", L)
  for (k in seq_len(L)) {
    if (isTRUE(verbose) && L > 1L) {
      message(sprintf("Fitting residualized-pair nuisance fold %d of %d.", k, L))
    }
    train_idx <- if (L == 1L) seq_along(Y) else which(folds != k)
    fold_fits[[k]] <- fit_rp_nuisance(
      Y = Y[train_idx],
      A = A[train_idx, , drop = FALSE],
      C = C[train_idx, , drop = FALSE],
      y_fitter = y_fitter,
      a_fitter = a_fitter,
      args_y = args_y,
      args_a = args_a,
      verbose = FALSE
    )
  }

  out <- list(folds = fold_fits, supplied = FALSE)
  class(out) <- "rp_crossfit_nuisance"
  out
}

as_rp_nuisance <- function(nuisance) {
  if (is.null(nuisance$Y_model)) {
    stop("'nuisance' must contain 'Y_model'.", call. = FALSE)
  }
  if (is.null(nuisance$A_models) && !is.null(nuisance$X_models)) {
    nuisance$A_models <- nuisance$X_models
  }
  if (is.null(nuisance$A_models)) {
    stop("'nuisance' must contain 'A_models'.", call. = FALSE)
  }
  if (is.null(nuisance$A_names)) {
    nuisance$A_names <- paste0("A", seq_along(nuisance$A_models))
  }
  class(nuisance) <- unique(c("rp_nuisance", class(nuisance)))
  nuisance
}
