#' Fit a Global Outcome Regression Object E(Y | X, C)
#'
#' @param Y Numeric vector of outcomes.
#' @param X Numeric matrix or data frame of observed treatments.
#' @param C Numeric matrix or data frame of observed confounders.
#' @param mu_fitter Function(Y, XC_df) that trains and returns a model.
#' @param ... Additional arguments passed to the inner fitter.
#' @return An S3 object of class "outcome_model".

outcome_model <- function(Y, X, C, mu_fitter, ...) {
  # 1. Capture original names BEFORE any coercion or processing
  orig_X_names <- if(is.matrix(X) || is.data.frame(X)) colnames(X) else NULL
  orig_C_names <- if(is.matrix(C) || is.data.frame(C)) colnames(C) else NULL

  # 2. Convert to data frames
  X_df <- as.data.frame(X)
  C_df <- as.data.frame(C)
  p <- ncol(X_df)
  q <- ncol(C_df)

  # 3. Apply the Contract: Use original names if they exist, otherwise generate ours
  # We use 'tx_' and 'conf_' prefixes to guarantee zero collision if we generate them
  if (is.null(orig_X_names)) {
    colnames(X_df) <- paste0("X", 1:p)
  } else {
    colnames(X_df) <- make.names(orig_X_names, unique = TRUE)
  }

  if (is.null(orig_C_names)) {
    colnames(C_df) <- paste0("C", 1:q)
  } else {
    colnames(C_df) <- make.names(orig_C_names, unique = TRUE)
  }

  # 4. Final safety check for overlaps
  overlapping <- intersect(colnames(X_df), colnames(C_df))
  if (length(overlapping) > 0) {
    stop("Overlapping column names detected: ", paste(overlapping, collapse = ", "))
  }

  # Fit using the consistent names
  df_train <- cbind(X_df, C_df)
  inner_fit <- mu_fitter(Y, df_train, ...)

  res <- list(
    inner_fit = inner_fit,
    X_names = colnames(X_df),
    C_names = colnames(C_df),
    p = p,
    q = q
  )
  class(res) <- "outcome_model"
  return(res)
}

#' Predict Method for Outcome Model
#'
#' @param object An object of class "outcome_model".
#' @param newdata A data frame containing the predictors.
#' @param ... Additional arguments passed to the inner predict method.
#'
#' @exportS3Method stats::predict
predict.outcome_model <- function(object, newdata, ...) {

  # Ensure newdata is a data frame
  newdata <- as.data.frame(newdata)

  # Combine required column names
  req_cols <- c(object$X_names, object$C_names)

  # 1. Check for missing columns
  missing_cols <- setdiff(req_cols, colnames(newdata))
  if (length(missing_cols) > 0) {
    stop(
      "The following required columns are missing from 'newdata': ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # 2. Subset and FORCE column order to match the training data exactly
  newdata <- newdata[, req_cols, drop = FALSE]

  # Predict using the inner model
  preds <- stats::predict(object$inner_fit, newdata = newdata, ...)

  # Extract numeric vector (SuperLearner returns a list with a $pred matrix)
  if (is.list(preds) && "pred" %in% names(preds)) {
    return(as.numeric(preds$pred))
  }

  return(as.numeric(preds))
}

# Low-level SuperLearner regression engine.
sl_regression_fitter <- function(Y,
                                 W,
                                 SL.library = csdr_default_sl_library(),
                                 family = stats::gaussian(),
                                 SL.lib = NULL,
                                 env = NULL,
                                 ...) {
  if (!is.null(SL.lib) && missing(SL.library)) {
    SL.library <- SL.lib
  }
  if (is.null(env)) {
    env <- asNamespace("SuperLearner")
  }
  SuperLearner::SuperLearner(
    Y = Y,
    X = as.data.frame(W),
    family = family,
    SL.library = resolve_available_sl_library(SL.library),
    env = env,
    ...
  )
}

# Backward-compatible wrapper.
SL_outcome_fitter <- function(Y,
                              XC_df,
                              SL.library = csdr_default_sl_library(),
                              family = stats::gaussian(),
                              SL.lib = NULL,
                              env = NULL,
                              ...) {
  sl_regression_fitter(
    Y = Y,
    W = XC_df,
    SL.library = SL.library,
    family = family,
    SL.lib = SL.lib,
    env = env,
    ...
  )
}

resolve_available_sl_library <- function(SL.library) {
  known_package_requirements <- c(
    SL.glmnet = "glmnet",
    SL.xgboost = "xgboost",
    SL.earth = "earth",
    SL.ranger = "ranger"
  )
  keep <- vapply(SL.library, function(wrapper) {
    pkg <- unname(known_package_requirements[wrapper])
    if (is.na(pkg)) {
      pkg <- NULL
    }
    if (!is.null(pkg) && !requireNamespace(pkg, quietly = TRUE)) {
      return(FALSE)
    }
    if (is.null(pkg)) {
      return(TRUE)
    }
    exists(wrapper, envir = asNamespace("SuperLearner"), inherits = FALSE)
  }, logical(1))
  out <- SL.library[keep]
  if (length(out) == 0L) {
    stop("No requested SuperLearner wrappers are available.", call. = FALSE)
  }
  out
}
