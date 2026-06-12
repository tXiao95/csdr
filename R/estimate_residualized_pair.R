#' Compute Residualized Pair
#'
#' @param Y Numeric vector of observed outcomes (length n).
#' @param X Numeric matrix or data frame of observed treatments (n x p).
#' @param C Numeric matrix or data frame of observed confounders (n x q).
#' @param C_models A list containing pre-fitted nuisance models:
#'        \code{C_models$Y_model} (predicts E(Y | C)) and
#'        \code{C_models$X_models} (a list of p models predicting E(X_i | C)).
#' @return A list containing two objects: \code{Ytilde} (numeric vector) and \code{Xtilde} (numeric matrix).

estimate_residualized_pair <- function(Y, X, C, C_models) {

  # ---------------------------------------------------------
  # Input Validation & Safe Column Naming
  # ---------------------------------------------------------
  orig_X_names <- if(is.matrix(X) || is.data.frame(X)) colnames(X) else NULL

  X_df <- as.data.frame(X)
  C_df <- as.data.frame(C)
  p <- ncol(X_df)
  n <- length(Y)

  # Safely apply names to X
  if (is.null(orig_X_names)) {
    colnames(X_df) <- paste0("X", 1:p)
  } else {
    colnames(X_df) <- make.names(orig_X_names, unique = TRUE)
  }

  # Check that C_models contains the expected structure
  if (is.null(C_models$Y_model) || is.null(C_models$X_models) || length(C_models$X_models) != p) {
    stop("C_models must be a list with a '$Y_model' and an '$X_models' list of length p.")
  }

  # ---------------------------------------------------------
  # 1. Predict E(Y | C) and compute \tilde{Y}
  # ---------------------------------------------------------
  E_Y_C <- stats::predict(C_models$Y_model, newdata = C_df)
  Ytilde <- Y - E_Y_C

  # ---------------------------------------------------------
  # 2. Predict E(X_i | C) and compute \tilde{X} for i = 1,...,p
  # ---------------------------------------------------------
  Xtilde <- matrix(NA, nrow = n, ncol = p)
  colnames(Xtilde) <- colnames(X_df)

  for (i in 1:p) {
    E_X_i_C <- stats::predict(C_models$X_models[[i]], newdata = C_df)
    Xtilde[, i] <- X_df[, i] - E_X_i_C
  }

  # ---------------------------------------------------------
  # 3. Return the residualized pair
  # ---------------------------------------------------------
  return(list(
    Ytilde = as.numeric(Ytilde),
    Xtilde = as.matrix(Xtilde)
  ))
}

#' Train all Nuisance Models E(Y | C) and E(X_i | C)
#'
#' @param Y Numeric vector of outcomes.
#' @param X Numeric matrix of treatments.
#' @param C Numeric matrix of confounders.
#' @param fitter The wrapper function to train the models (e.g., SL_nuisance_fitter).
#' @param ... Additional arguments passed to the fitter.
#' @return A list containing $Y_model and $X_models.

train_nuisance_models <- function(Y, X, C, fitter, ...) {
  X_df <- as.data.frame(X)
  C_df <- as.data.frame(C)
  p <- ncol(X_df)

  # Train E(Y | C)
  Y_model <- nuisance_C_model(target = Y, C = C_df, fitter = fitter, ...)

  # Train E(X_i | C) for all p dimensions
  X_models <- vector("list", p)
  for (i in 1:p) {
    message("Training Model X", i, " | C")
    X_models[[i]] <- nuisance_C_model(target = X_df[, i], C = C_df, fitter = fitter, ...)
  }

  return(list(
    Y_model = Y_model,
    X_models = X_models
  ))
}
