#' Fit a Global Propensity Score Model f(X | C)
#'
#' @param X Numeric matrix or data frame of observed treatments.
#' @param C Numeric matrix or data frame of observed confounders.
#' @param pi_fitter Function(X_df, C_df, ...) that trains and returns a density model.
#' @param ... Additional arguments passed to the pi_fitter.
#' @return An S3 object of class "gps_model".

gps_model <- function(X, C, pi_fitter, ...) {
  # 1. Capture original names BEFORE any coercion or processing
  orig_X_names <- if(is.matrix(X) || is.data.frame(X)) colnames(X) else NULL
  orig_C_names <- if(is.matrix(C) || is.data.frame(C)) colnames(C) else NULL

  # 2. Convert to data frames
  X_df <- as.data.frame(X)
  C_df <- as.data.frame(C)
  p <- ncol(X_df)
  q <- ncol(C_df)

  # 3. Apply the Contract: Use original names if they exist, otherwise generate ours
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

  # 5. Fit using the consistent names
  # Note: Unlike outcome_model which binds cbind(X, C) to predict Y,
  # a GPS fitter needs to predict X from C, so we pass them separately.
  inner_fit <- pi_fitter(X = X_df, C = C_df, ...)

  res <- list(
    inner_fit = inner_fit,
    X_names = colnames(X_df),
    C_names = colnames(C_df),
    p = p,
    q = q
  )
  class(res) <- "gps_model"
  return(res)
}

#' Predict Method for GPS Model
#'
#' @param object An object of class "gps_model".
#' @param newdata A data frame containing predictors for the nuisance models.
#' @param ... Additional arguments passed to the inner predict method.
#'
#' @exportS3Method stats::predict
predict.gps_model <- function(object, newdata, ...) {
  newdata <- as.data.frame(newdata)
  req_cols <- c(object$X_names, object$C_names)

  # 1. Check for missing columns
  missing_cols <- setdiff(req_cols, colnames(newdata))
  if (length(missing_cols) > 0) {
    stop("The following required columns are missing from 'newdata': ",
         paste(missing_cols, collapse = ", "))
  }

  # 2. Subset and FORCE column order to match the training data exactly
  newdata <- newdata[, req_cols, drop = FALSE]

  # Predict using the inner model
  preds <- stats::predict(object$inner_fit, newdata = newdata, ...)

  return(as.numeric(preds))
}

#' Inner Fitter for MVN GPS
#'
#' @param A Numeric matrix or data frame of observed treatments.
#' @param C Numeric matrix or data frame of observed confounders.
#' @param X Compatibility alias for `A`.
#' @param method_gps String in `c("linear","SuperLearner")`.
#' @param ... Additional arguments passed to the density fitter.
mvn_gps_fitter <- function(A = NULL, C, X = NULL,
                           method_gps = c("linear", "SuperLearner"),
                           delta_n = 1e-16, ...) {
  if (is.null(A)) {
    A <- X
  }
  if (is.null(A)) {
    stop("'A' must be supplied.", call. = FALSE)
  }
  if (!is.numeric(delta_n) || length(delta_n) != 1L ||
      is.na(delta_n) || !is.finite(delta_n) || delta_n <= 0) {
    stop("'delta_n' must be a positive finite number.", call. = FALSE)
  }
  method_gps <- match.arg(method_gps)

  # Data is already clean from the wrapper
  p <- ncol(A)
  X_mat <- as.matrix(A)
  C_df <- C # SuperLearner and lm both prefer data frames for predictors

  if (method_gps == "linear") {
    inner_fit <- stats::lm(X_mat ~ ., data = C_df)
    resids <- stats::residuals(inner_fit)

  } else if (method_gps == "SuperLearner") {
    inner_fit <- list()
    resids <- matrix(NA, nrow = nrow(X_mat), ncol = p)

    for (j in 1:p) {
      sl_fit <- sl_regression_fitter(
        Y = X_mat[, j],
        W = C_df,
        ...
      )
      inner_fit[[paste0("X", j)]] <- sl_fit
      resids[, j] <- X_mat[, j] - sl_fit$SL.predict
    }
  }

  sigma_hat <- as.matrix(stats::cov(as.matrix(resids)))

  res <- list(
    inner_fit = inner_fit,
    sigma_hat = sigma_hat,
    method = method_gps,
    p = p,
    X_names = colnames(A),
    C_names = colnames(C),
    density_floor = delta_n
  )
  class(res) <- "mvn_inner"
  return(res)
}

# Backward-compatible MVN GPS fitter name.
mvn_fitter <- function(X, C, method_gps = c("linear", "SuperLearner"), ...) {
  mvn_gps_fitter(A = X, C = C, method_gps = method_gps, ...)
}

#' Predict Method for Inner MVN
#'
#' @param object An object of class "mvn_inner".
#' @param newdata A data frame containing predictors.
#' @param delta_n Floor for density predictions.
#' @param ... Additional arguments passed to the inner predict method.
#'
#' @exportS3Method stats::predict

predict.mvn_inner <- function(object, newdata, delta_n = object$density_floor %||% 1e-16, ...) {
  # newdata is already clean and correctly ordered by predict.gps_model
  X_new <- as.matrix(newdata[, object$X_names, drop = FALSE])
  C_new <- newdata[, object$C_names, drop = FALSE]

  if (object$method == "linear") {
    mu_hat <- as.matrix(stats::predict(object$inner_fit, newdata = C_new))

  } else if (object$method == "SuperLearner") {
    mu_hat <- matrix(NA, nrow = nrow(C_new), ncol = object$p)
    for (j in 1:object$p) {
      mu_hat[, j] <- stats::predict(object$inner_fit[[j]], newdata = C_new)$pred
    }
  }

  centered_x <- X_new - mu_hat

  raw_f_hat <- mvtnorm::dmvnorm(x = centered_x,
                                mean = rep(0, object$p),
                                sigma = object$sigma_hat)

  return(pmax(raw_f_hat, delta_n))
}
