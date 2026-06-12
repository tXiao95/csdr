#' Fit a Marginal Nuisance Regression Object E(Target | C)
#'
#' @param target Numeric vector of the target variable (can be Y or a specific X dimension).
#' @param C Numeric matrix or data frame of observed confounders.
#' @param fitter Function(target, C_df) that trains and returns a model.
#' @param ... Additional arguments passed to the fitter.
#' @return An S3 object of class "nuisance_C_model".

nuisance_C_model <- function(target, C, fitter, ...) {
  # 1. Capture original names BEFORE any coercion or processing
  orig_C_names <- if(is.matrix(C) || is.data.frame(C)) colnames(C) else NULL

  # 2. Convert to data frame
  C_df <- as.data.frame(C)
  q <- ncol(C_df)

  # 3. Apply the Contract: Use original names if they exist, otherwise generate ours
  if (is.null(orig_C_names)) {
    colnames(C_df) <- paste0("C", 1:q)
  } else {
    colnames(C_df) <- make.names(orig_C_names, unique = TRUE)
  }

  # Fit using the consistent names
  inner_fit <- fitter(target, C_df, ...)

  res <- list(
    inner_fit = inner_fit,
    C_names = colnames(C_df),
    q = q
  )

  class(res) <- "nuisance_C_model"
  return(res)
}

#' Predict Method for Nuisance C Model
#'
#' @param object An object of class "nuisance_C_model".
#' @param newdata A data frame containing the confounders (C).
#' @param ... Additional arguments passed to the underlying predict method.
#' 
#' @return A numeric vector of predictions.

predict.nuisance_C_model <- function(object, newdata, ...) {

  # Ensure newdata is a data frame
  newdata <- as.data.frame(newdata)

  # Combine required column names
  req_cols <- object$C_names

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

# Define the wrapper for SuperLearner
SL_nuisance_fitter <- function(target, C_df, SL.lib = c("SL.glm", "SL.glmnet", "SL.xgboost", "SL.earth"), ...) {
  SuperLearner::SuperLearner(Y = target, X = C_df, family = stats::gaussian(), SL.lib = SL.lib, ...)
}
