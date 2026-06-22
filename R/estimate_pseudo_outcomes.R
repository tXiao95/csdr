#' Compute pseudo-outcomes from numeric nuisance predictions
#'
#' @param Y Numeric outcome vector.
#' @param m_obs Numeric vector of observed outcome regression predictions.
#' @param pi_obs Numeric vector of observed GPS predictions.
#' @param m_marginal Numeric vector of marginal outcome regression predictions.
#' @param pi_marginal Numeric vector of marginal GPS predictions.
#' @param gps_floor Positive floor applied to observed GPS predictions.
#'
#' @return A numeric vector of pseudo-outcomes.
compute_pseudo_outcomes <- function(Y, m_obs, pi_obs, m_marginal, pi_marginal,
                                    gps_floor = 1e-8) {
  inputs <- list(
    Y = Y,
    m_obs = m_obs,
    pi_obs = pi_obs,
    m_marginal = m_marginal,
    pi_marginal = pi_marginal
  )
  numeric_vectors <- vapply(
    inputs,
    function(x) is.numeric(x) && is.null(dim(x)),
    logical(1L)
  )
  if (!all(numeric_vectors)) {
    stop(
      "'Y', 'm_obs', 'pi_obs', 'm_marginal', and 'pi_marginal' must be numeric vectors.",
      call. = FALSE
    )
  }
  lengths <- vapply(inputs, length, integer(1L))

  if (length(unique(lengths)) != 1L) {
    stop(
      "'Y', 'm_obs', 'pi_obs', 'm_marginal', and 'pi_marginal' must have the same length.",
      call. = FALSE
    )
  }
  if (!is.numeric(gps_floor) || length(gps_floor) != 1L ||
      is.na(gps_floor) || !is.finite(gps_floor) || gps_floor <= 0) {
    stop("'gps_floor' must be a positive finite number.", call. = FALSE)
  }

  if (any(vapply(inputs, function(x) anyNA(x), logical(1L)))) {
    stop("Pseudo-outcome inputs must not contain missing values.", call. = FALSE)
  }

  with(inputs, {
    pi_marginal * (Y - m_obs) / pmax(pi_obs, gps_floor) + m_marginal
  })
}

#' Fit nuisance models for pseudo-outcome construction
#'
#' @param Y Numeric outcome vector.
#' @param A Numeric matrix or data frame of observed exposures.
#' @param C Numeric matrix or data frame of covariates.
#' @param outcome_fitter Function used by [outcome_model()].
#' @param gps_fitter Function used by [gps_model()].
#' @param args_outcome Additional arguments passed to `outcome_fitter`.
#' @param args_gps Additional arguments passed to `gps_fitter`.
#' @param verbose Logical; if `TRUE`, print progress messages.
#'
#' @return A fitted pseudo-outcome nuisance object.
fit_po_nuisance <- function(Y, A, C,
                            outcome_fitter = SL_outcome_fitter,
                            gps_fitter = mvn_fitter,
                            args_outcome = list(),
                            args_gps = list(),
                            verbose = TRUE) {
  data <- normalize_ers_inputs(Y = Y, A = A, C = C, a_eval = A)

  if (!is.function(outcome_fitter)) {
    stop("'outcome_fitter' must be a function.", call. = FALSE)
  }
  if (!is.function(gps_fitter)) {
    stop("'gps_fitter' must be a function.", call. = FALSE)
  }

  if (isTRUE(verbose)) {
    message("Fitting pseudo-outcome nuisance models.")
  }

  outcome_args <- list(
    Y = data$Y,
    X = data$A,
    C = data$C,
    mu_fitter = outcome_fitter
  )
  gps_args <- list(
    X = data$A,
    C = data$C,
    pi_fitter = gps_fitter
  )

  out <- list(
    outcome_model = do.call(outcome_model, c(outcome_args, args_outcome)),
    gps_model = do.call(gps_model, c(gps_args, args_gps)),
    A_names = colnames(data$A),
    C_names = colnames(data$C)
  )
  class(out) <- "po_nuisance"
  out
}

#' Predict nuisance quantities for pseudo-outcome construction
#'
#' @param nuisance A fitted pseudo-outcome nuisance object.
#' @param A Numeric matrix or data frame of observed exposures to predict.
#' @param C Numeric matrix or data frame of covariates paired with `A`.
#' @param C_marginal Covariate rows used for empirical marginalization.
#' @param gps_floor Positive floor applied to GPS predictions.
#' @param chunk_size Optional chunk size for future batched prediction.
#' @param verbose Logical; if `TRUE`, print progress messages.
#'
#' @return A list with `m_obs`, `pi_obs`, `m_marginal`, and `pi_marginal`.
predict_po_nuisance <- function(nuisance, A, C, C_marginal = C,
                                gps_floor = 1e-8, chunk_size = NULL,
                                verbose = FALSE) {
  if (is.null(nuisance$outcome_model) || is.null(nuisance$gps_model)) {
    stop("'nuisance' must contain 'outcome_model' and 'gps_model'.", call. = FALSE)
  }
  if (!is.numeric(gps_floor) || length(gps_floor) != 1L ||
      is.na(gps_floor) || !is.finite(gps_floor) || gps_floor <= 0) {
    stop("'gps_floor' must be a positive finite number.", call. = FALSE)
  }
  if (!is.null(chunk_size) &&
      (!is.numeric(chunk_size) || length(chunk_size) != 1L ||
       is.na(chunk_size) || chunk_size < 1L)) {
    stop("'chunk_size' must be NULL or a positive integer.", call. = FALSE)
  }

  C_df <- normalize_target_table(C, n = NULL, arg = "C", prefix = "C")
  A_df <- normalize_target_table(A, n = nrow(C_df), arg = "A", prefix = "A")
  C_marginal_df <- normalize_target_table(
    C_marginal,
    n = NULL,
    arg = "C_marginal",
    prefix = "C"
  )
  if (nrow(C_marginal_df) < 1L) {
    stop("'C_marginal' must contain at least one row.", call. = FALSE)
  }

  obs_newdata <- make_observed_newdata(
    model = nuisance$outcome_model,
    A = A_df,
    C = C_df
  )
  m_obs <- as.numeric(stats::predict(nuisance$outcome_model, newdata = obs_newdata))

  gps_newdata <- make_observed_newdata(
    model = nuisance$gps_model,
    A = A_df,
    C = C_df
  )
  pi_obs <- pmax(
    as.numeric(stats::predict(nuisance$gps_model, newdata = gps_newdata)),
    gps_floor
  )

  A_for_grid <- A_df
  C_for_grid <- C_marginal_df
  colnames(A_for_grid) <- nuisance$outcome_model$X_names
  colnames(C_for_grid) <- nuisance$outcome_model$C_names

  x_names <- nuisance$outcome_model$X_names
  dt_grid <- data.table::as.data.table(C_for_grid)
  for (x_col in x_names) {
    data.table::set(dt_grid, j = x_col, value = 0.0)
  }
  data.table::setcolorder(dt_grid, c(x_names, nuisance$outcome_model$C_names))

  n <- nrow(A_for_grid)
  m_marginal <- pi_marginal <- rep(NA_real_, n)
  for (j in seq_len(n)) {
    if (isTRUE(verbose) && (j == 1L || j == n || j %% 50L == 0L)) {
      message(sprintf("Predicting pseudo-outcome marginal row %d of %d.", j, n))
    }
    for (k in seq_along(x_names)) {
      data.table::set(dt_grid, j = x_names[k], value = A_for_grid[[k]][j])
    }
    m_marginal[j] <- mean(
      as.numeric(stats::predict(nuisance$outcome_model, newdata = dt_grid))
    )
    pi_marginal[j] <- mean(pmax(
      as.numeric(stats::predict(nuisance$gps_model, newdata = dt_grid)),
      gps_floor
    ))
  }

  list(
    m_obs = m_obs,
    pi_obs = pi_obs,
    m_marginal = m_marginal,
    pi_marginal = pi_marginal
  )
}

#' Estimate pseudo-outcomes for CSDR target construction
#'
#' @param Y Numeric outcome vector.
#' @param A Numeric matrix or data frame of observed exposures.
#' @param C Numeric matrix or data frame of covariates.
#' @param L Number of folds. `L = 1` disables cross-fitting.
#' @param folds Optional fold assignments.
#' @param outcome_fitter Function used by [outcome_model()].
#' @param gps_fitter Function used by [gps_model()].
#' @param nuisance Optional pre-fitted pseudo-outcome nuisance object.
#' @param C_marginal Covariate rows used for empirical marginalization.
#' @param gps_floor Positive floor applied to GPS predictions.
#' @param seed Optional random seed for fold generation.
#' @param return_nuisance Logical; if `TRUE`, return nuisance fits.
#' @param verbose Logical; if `TRUE`, print progress messages.
#' @param args_outcome Additional arguments passed to `outcome_fitter`.
#' @param args_gps Additional arguments passed to `gps_fitter`.
#' @param X Deprecated compatibility alias for `A`.
#' @param out_model Deprecated pre-fitted outcome model argument.
#' @param gps_model Deprecated pre-fitted GPS model argument.
#' @param delta_n Deprecated compatibility alias for `gps_floor`.
#' @param chunk_size Optional chunk size for future batched prediction.
#'
#' @return An S3 object of class `"po_fit"` and `"csdr_target"`.
#'
#' @export
estimate_pseudo_outcomes <- function(
  Y,
  A = NULL,
  C,
  L = 5,
  folds = NULL,
  outcome_fitter = SL_outcome_fitter,
  gps_fitter = mvn_fitter,
  nuisance = NULL,
  C_marginal = C,
  gps_floor = 1e-8,
  seed = NULL,
  return_nuisance = FALSE,
  verbose = TRUE,
  args_outcome = list(),
  args_gps = list(),
  X = NULL,
  out_model = NULL,
  gps_model = NULL,
  delta_n = NULL,
  chunk_size = NULL
) {
  call <- match.call()
  legacy_return <- is.null(A) && !is.null(X) &&
    !is.null(out_model) && !is.null(gps_model)
  if (is.null(A)) {
    A <- X
  }
  if (!is.null(delta_n)) {
    gps_floor <- delta_n
  }
  if (is.null(nuisance) && (!is.null(out_model) || !is.null(gps_model))) {
    nuisance <- list(
      outcome_model = out_model,
      gps_model = gps_model,
      A_names = out_model$X_names,
      C_names = out_model$C_names
    )
    class(nuisance) <- "po_nuisance"
  }

  data <- normalize_ers_inputs(Y = Y, A = A, C = C, a_eval = A)
  C_marginal_df <- normalize_target_table(
    C_marginal,
    n = NULL,
    arg = "C_marginal",
    prefix = "C"
  )
  validate_ers_options(L = L, gps_floor = gps_floor, normalize_ipw = TRUE)
  folds <- make_folds(n = length(data$Y), L = L, folds = folds, seed = seed)
  L_eff <- length(unique(folds))

  nuisance_fits <- prepare_po_nuisance_fits(
    nuisance = nuisance,
    Y = data$Y,
    A = data$A,
    C = data$C,
    L = L_eff,
    folds = folds,
    outcome_fitter = outcome_fitter,
    gps_fitter = gps_fitter,
    args_outcome = args_outcome,
    args_gps = args_gps,
    verbose = verbose
  )

  n <- length(data$Y)
  m_obs <- pi_obs <- m_marginal <- pi_marginal <- rep(NA_real_, n)
  for (k in seq_len(L_eff)) {
    test_idx <- if (L_eff == 1L) seq_len(n) else which(folds == k)
    pred <- predict_po_nuisance(
      nuisance = nuisance_fits$folds[[k]],
      A = data$A[test_idx, , drop = FALSE],
      C = data$C[test_idx, , drop = FALSE],
      C_marginal = C_marginal_df,
      gps_floor = gps_floor,
      chunk_size = chunk_size,
      verbose = verbose && L_eff > 1L
    )
    m_obs[test_idx] <- pred$m_obs
    pi_obs[test_idx] <- pred$pi_obs
    m_marginal[test_idx] <- pred$m_marginal
    pi_marginal[test_idx] <- pred$pi_marginal
  }

  xi_hat <- compute_pseudo_outcomes(
    Y = data$Y,
    m_obs = m_obs,
    pi_obs = pi_obs,
    m_marginal = m_marginal,
    pi_marginal = pi_marginal,
    gps_floor = gps_floor
  )

  if (isTRUE(legacy_return)) {
    return(xi_hat)
  }

  out <- list(
    pseudo_outcomes = xi_hat,
    Y_tilde = xi_hat,
    A_tilde = as.matrix(data$A),
    folds = folds,
    nuisance = if (return_nuisance) nuisance_fits else NULL,
    call = call,
    diagnostics = list(
      method = "PO",
      L = L_eff,
      requested_L = L,
      n = n,
      p = ncol(data$A),
      gps_summary = summary(pi_obs)
    )
  )
  class(out) <- c("po_fit", "csdr_target")
  out
}

prepare_po_nuisance_fits <- function(nuisance, Y, A, C, L, folds,
                                     outcome_fitter, gps_fitter,
                                     args_outcome, args_gps, verbose) {
  if (!is.null(nuisance)) {
    fold_fits <- if (!is.null(nuisance$folds)) {
      nuisance$folds
    } else if (!is.null(nuisance$fold_nuisance)) {
      nuisance$fold_nuisance
    } else {
      rep(list(nuisance), L)
    }
    if (length(fold_fits) != L) {
      stop("Supplied pseudo-outcome nuisance fits do not match the number of folds.", call. = FALSE)
    }
    return(list(folds = fold_fits, supplied = TRUE))
  }

  fold_fits <- vector("list", L)
  for (k in seq_len(L)) {
    if (isTRUE(verbose) && L > 1L) {
      message(sprintf("Fitting pseudo-outcome nuisance fold %d of %d.", k, L))
    }
    train_idx <- if (L == 1L) seq_along(Y) else which(folds != k)
    fold_fits[[k]] <- fit_po_nuisance(
      Y = Y[train_idx],
      A = A[train_idx, , drop = FALSE],
      C = C[train_idx, , drop = FALSE],
      outcome_fitter = outcome_fitter,
      gps_fitter = gps_fitter,
      args_outcome = args_outcome,
      args_gps = args_gps,
      verbose = FALSE
    )
  }

  out <- list(folds = fold_fits, supplied = FALSE)
  class(out) <- "po_crossfit_nuisance"
  out
}

normalize_target_table <- function(x, n = NULL, arg, prefix) {
  if (is.null(x)) {
    stop(sprintf("'%s' cannot be NULL.", arg), call. = FALSE)
  }
  if (is.null(dim(x))) {
    x <- matrix(x, ncol = 1L)
  }
  df <- as.data.frame(x)
  if (!is.null(n) && nrow(df) != n) {
    stop(sprintf("'%s' must have %d rows.", arg, n), call. = FALSE)
  }
  if (ncol(df) < 1L) {
    stop(sprintf("'%s' must have at least one column.", arg), call. = FALSE)
  }
  if (is.null(colnames(df))) {
    colnames(df) <- paste0(prefix, seq_len(ncol(df)))
  } else {
    colnames(df) <- make.names(colnames(df), unique = TRUE)
  }
  df[] <- lapply(df, as.numeric)
  if (anyNA(df)) {
    stop(sprintf("'%s' must contain only non-missing numeric values.", arg), call. = FALSE)
  }
  df
}
