#' Estimate a causal exposure-response surface
#'
#' `estimate_ers()` estimates the causal exposure-response surface
#' \eqn{\mu(a) = E[Y^a]} for a continuous, possibly multivariate exposure.
#' It supports regression adjustment (`"RA"`), inverse probability weighting
#' (`"IPW"`), and doubly robust (`"DR"`) estimators. Nuisance models are fit
#' either on the full data (`L = 1`) or by cross-fitting (`L > 1`).
#'
#' @param Y Numeric outcome vector of length `n`.
#' @param A Numeric matrix or data frame of observed exposures/treatments
#'   with `n` rows.
#' @param C Numeric matrix or data frame of covariates with `n` rows.
#' @param a_eval Optional matrix or data frame of exposure values where the
#'   exposure-response surface should be evaluated. Defaults to the observed
#'   rows of `A`.
#' @param estimator Estimator to use: `"DR"`, `"RA"`, or `"IPW"`.
#' @param L Number of folds for cross-fitting. `L = 1` disables cross-fitting
#'   and fits nuisances on the full data.
#' @param folds Optional fold assignments of length `n`. If supplied, these
#'   are used instead of randomly generated folds.
#' @param outcome_fitter Function used by [outcome_model()] to fit
#'   \eqn{E[Y \mid A, C]}. Required for `"RA"` and `"DR"`.
#' @param gps_fitter Function used by [gps_model()] to fit the conditional
#'   exposure density \eqn{f(A \mid C)}. Required for `"IPW"` and `"DR"`.
#' @param h Optional numeric bandwidth vector. If `NULL`, defaults to
#'   `c_multiplier * apply(A, 2, sd) * n^(-0.2)`.
#' @param c_multiplier Multiplier used for the default bandwidth.
#' @param gps_floor Positive floor applied to GPS estimates.
#' @param normalize_ipw Logical; if `TRUE`, use self-normalized IPW weights.
#' @param seed Optional random seed used when `folds` is not supplied.
#' @param return_nuisance Logical; if `TRUE`, return fitted nuisance objects
#'   and out-of-fold predictions.
#' @param verbose Logical; if `TRUE`, print progress messages.
#' @param args_outcome Additional arguments passed to `outcome_fitter`.
#' @param args_gps Additional arguments passed to `gps_fitter`.
#'
#' @return An S3 object of class `"ers_fit"`.
#'
#' @rdname estimate-ers
#' @export
estimate_ers <- function(
  Y,
  A,
  C,
  a_eval = NULL,
  estimator = c("DR", "RA", "IPW"),
  L = 5,
  folds = NULL,
  outcome_fitter = SL_outcome_fitter,
  gps_fitter = mvn_fitter,
  h = NULL,
  c_multiplier = 1.25,
  gps_floor = 1e-8,
  normalize_ipw = TRUE,
  seed = NULL,
  return_nuisance = FALSE,
  verbose = TRUE,
  args_outcome = list(),
  args_gps = list()
) {
  call <- match.call()
  estimator <- match.arg(estimator)

  data <- normalize_ers_inputs(Y = Y, A = A, C = C, a_eval = a_eval)
  validate_ers_options(
    L = L,
    gps_floor = gps_floor,
    normalize_ipw = normalize_ipw
  )
  folds <- make_folds(n = length(data$Y), L = L, folds = folds, seed = seed)
  L_eff <- length(unique(folds))

  if (verbose && L_eff == 1L) {
    message("Fitting ERS nuisances on the full data (L = 1).")
  } else if (verbose) {
    message(sprintf("Fitting ERS nuisances with %d-fold cross-fitting.", L_eff))
  }

  nuisance_fit <- fit_ers_nuisance(
    Y = data$Y,
    A = data$A,
    C = data$C,
    estimator = estimator,
    L = L_eff,
    folds = folds,
    outcome_fitter = outcome_fitter,
    gps_fitter = gps_fitter,
    args_outcome = args_outcome,
    args_gps = args_gps,
    verbose = verbose
  )

  nuisance_pred <- predict_ers_nuisance(
    nuisance = nuisance_fit,
    Y = data$Y,
    A = data$A,
    C = data$C,
    a_eval = data$a_eval,
    estimator = estimator,
    gps_floor = gps_floor,
    verbose = verbose
  )

  ers <- compute_ers(
    Y = data$Y,
    A = data$A,
    a_eval = data$a_eval,
    estimator = estimator,
    m_eval = nuisance_pred$m_eval,
    pi_obs = nuisance_pred$pi_obs,
    h = h,
    c_multiplier = c_multiplier,
    gps_floor = gps_floor,
    normalize_ipw = normalize_ipw
  )

  out <- list(
    results = ers$results,
    nuisance = if (return_nuisance) nuisance_pred else NULL,
    estimator = estimator,
    L = L_eff,
    folds = folds,
    h = ers$h,
    call = call,
    diagnostics = utils::modifyList(
      list(
        n = length(data$Y),
        n_eval = nrow(data$a_eval),
        p = ncol(data$A),
        requested_L = L,
        nuisance_strategy = if (L_eff == 1L) "fullfit" else "crossfit",
        seed = seed,
        return_nuisance = return_nuisance
      ),
      ers$diagnostics
    )
  )
  class(out) <- "ers_fit"
  out
}

#' @keywords internal
make_folds <- function(n, L, folds = NULL, seed = NULL) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1L) {
    stop("'n' must be a positive integer.", call. = FALSE)
  }
  n <- as.integer(n)

  if (!is.numeric(L) || length(L) != 1L || is.na(L) || L < 1L) {
    stop("'L' must be a positive integer.", call. = FALSE)
  }
  L <- as.integer(L)
  if (L > n) {
    stop("'L' cannot exceed the number of observations.", call. = FALSE)
  }

  if (!is.null(folds)) {
    if (length(folds) != n) {
      stop("'folds' must have length n.", call. = FALSE)
    }
    if (anyNA(folds)) {
      stop("'folds' cannot contain missing values.", call. = FALSE)
    }
    fold_levels <- unique(folds)
    if (L > 1L && length(fold_levels) < 2L) {
      stop("'folds' must contain at least two folds when L > 1.", call. = FALSE)
    }
    return(as.integer(match(folds, fold_levels)))
  }

  if (L == 1L) {
    return(rep.int(1L, n))
  }

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  sample(rep(seq_len(L), length.out = n))
}

#' @keywords internal
fit_ers_nuisance <- function(Y, A, C, estimator, L, folds,
                             outcome_fitter, gps_fitter,
                             args_outcome, args_gps, verbose) {
  estimator <- match.arg(estimator, choices = c("DR", "RA", "IPW"))

  if (estimator %in% c("RA", "DR") && !is.function(outcome_fitter)) {
    stop("'outcome_fitter' must be a function for RA and DR.", call. = FALSE)
  }
  if (estimator %in% c("IPW", "DR") && !is.function(gps_fitter)) {
    stop("'gps_fitter' must be a function for IPW and DR.", call. = FALSE)
  }

  outcome_models <- vector("list", L)
  gps_models <- vector("list", L)

  for (k in seq_len(L)) {
    if (verbose && L > 1L) {
      message(sprintf("Fitting nuisance fold %d of %d.", k, L))
    }

    train_idx <- if (L == 1L) seq_along(Y) else which(folds != k)

    if (estimator %in% c("RA", "DR")) {
      outcome_args <- list(
        Y = Y[train_idx],
        X = A[train_idx, , drop = FALSE],
        C = C[train_idx, , drop = FALSE],
        mu_fitter = outcome_fitter
      )
      outcome_models[[k]] <- do.call(outcome_model, c(outcome_args, args_outcome))
    }

    if (estimator %in% c("IPW", "DR")) {
      gps_args <- list(
        X = A[train_idx, , drop = FALSE],
        C = C[train_idx, , drop = FALSE],
        pi_fitter = gps_fitter
      )
      gps_models[[k]] <- do.call(gps_model, c(gps_args, args_gps))
    }
  }

  out <- list(
    estimator = estimator,
    L = L,
    folds = folds,
    outcome_models = outcome_models,
    gps_models = gps_models
  )
  class(out) <- "ers_nuisance_fit"
  out
}

#' @keywords internal
predict_ers_nuisance <- function(nuisance, Y, A, C, a_eval, estimator,
                                 gps_floor, verbose) {
  estimator <- match.arg(estimator, choices = c("DR", "RA", "IPW"))
  n <- length(Y)
  m <- nrow(a_eval)
  L <- nuisance$L
  folds <- nuisance$folds

  m_eval <- if (estimator %in% c("RA", "DR")) {
    matrix(NA_real_, nrow = n, ncol = m)
  } else {
    NULL
  }
  pi_obs <- if (estimator %in% c("IPW", "DR")) rep(NA_real_, n) else NULL

  for (k in seq_len(L)) {
    if (verbose && L > 1L) {
      message(sprintf("Predicting nuisance fold %d of %d.", k, L))
    }

    test_idx <- if (L == 1L) seq_len(n) else which(folds == k)

    if (estimator %in% c("RA", "DR")) {
      out_model <- nuisance$outcome_models[[k]]
      m_eval[test_idx, ] <- predict_outcome_eval_matrix(
        out_model = out_model,
        C = C[test_idx, , drop = FALSE],
        a_eval = a_eval
      )
    }

    if (estimator %in% c("IPW", "DR")) {
      gps_mod <- nuisance$gps_models[[k]]
      newdata <- make_observed_newdata(
        model = gps_mod,
        A = A[test_idx, , drop = FALSE],
        C = C[test_idx, , drop = FALSE]
      )
      pi_obs[test_idx] <- pmax(
        stats::predict(gps_mod, newdata = newdata),
        gps_floor
      )
    }
  }

  out <- list(
    fit = nuisance,
    m_eval = m_eval,
    pi_obs = pi_obs,
    diagnostics = list(
      min_gps = if (is.null(pi_obs)) NA_real_ else min(pi_obs, na.rm = TRUE)
    )
  )
  class(out) <- "ers_nuisance_predictions"
  out
}

#' @keywords internal
compute_ers <- function(Y, A, a_eval, estimator,
                        m_eval = NULL, pi_obs = NULL,
                        h = NULL, c_multiplier = 1.25,
                        gps_floor = 1e-8, normalize_ipw = TRUE,
                        level = 0.95) {
  estimator <- match.arg(estimator, choices = c("DR", "RA", "IPW"))
  Y <- as.numeric(Y)
  A_df <- as.data.frame(A)
  a_eval_df <- as.data.frame(a_eval)
  A_mat <- as.matrix(A_df)
  a_eval_mat <- as.matrix(a_eval_df)
  n <- length(Y)
  p <- ncol(A_mat)
  m <- nrow(a_eval_mat)

  if (estimator %in% c("RA", "DR")) {
    if (is.null(m_eval) || !identical(dim(m_eval), c(n, m))) {
      stop("'m_eval' must be an n by nrow(a_eval) matrix for RA and DR.", call. = FALSE)
    }
  }
  if (estimator %in% c("IPW", "DR")) {
    if (is.null(pi_obs) || length(pi_obs) != n) {
      stop("'pi_obs' must have length n for IPW and DR.", call. = FALSE)
    }
    pi_obs <- pmax(as.numeric(pi_obs), gps_floor)
  }

  h_used <- resolve_ers_bandwidth(
    A = A_mat,
    h = h,
    c_multiplier = c_multiplier,
    require_positive = estimator %in% c("IPW", "DR")
  )
  h_names <- paste0("h_used_", seq_len(p))
  z_value <- stats::qnorm(1 - (1 - level) / 2)

  estimate <- se <- ci_lower <- ci_upper <- effective_n <- rep(NA_real_, m)
  min_gps <- if (estimator %in% c("IPW", "DR")) min(pi_obs, na.rm = TRUE) else NA_real_

  for (j in seq_len(m)) {
    if (estimator == "RA") {
      m_j <- m_eval[, j]
      estimate[j] <- mean(m_j)
      psi <- m_j - estimate[j]
      se[j] <- stats::sd(psi) / sqrt(n)
      effective_n[j] <- n
    } else {
      kernel_weights <- product_gaussian_kernel(
        A = A_mat,
        a = as.numeric(a_eval_mat[j, ]),
        h = h_used
      )
      weighted_kernel <- kernel_weights / pi_obs
      denom <- sum(weighted_kernel)
      effective_n[j] <- if (sum(weighted_kernel^2) > 0) {
        denom^2 / sum(weighted_kernel^2)
      } else {
        0
      }

      if (!is.finite(denom) || denom <= .Machine$double.eps) {
        estimate[j] <- NA_real_
        se[j] <- NA_real_
      } else if (estimator == "IPW") {
        if (isTRUE(normalize_ipw)) {
          w <- n * weighted_kernel / denom
          estimate[j] <- mean(w * Y)
          psi <- w * (Y - estimate[j])
        } else {
          estimate[j] <- mean(weighted_kernel * Y)
          psi <- weighted_kernel * Y - estimate[j]
        }
        se[j] <- stats::sd(psi) / sqrt(n)
      } else if (estimator == "DR") {
        m_j <- m_eval[, j]
        ra_est <- mean(m_j)
        if (isTRUE(normalize_ipw)) {
          w <- n * weighted_kernel / denom
          delta <- mean(w * (Y - m_j))
          estimate[j] <- ra_est + delta
          psi <- w * (Y - m_j) + (m_j - estimate[j])
        } else {
          delta <- mean(weighted_kernel * (Y - m_j))
          estimate[j] <- ra_est + delta
          psi <- (m_j - ra_est) + weighted_kernel * (Y - m_j) - delta
        }
        se[j] <- stats::sd(psi) / sqrt(n)
      }
    }

    ci_lower[j] <- estimate[j] - z_value * se[j]
    ci_upper[j] <- estimate[j] + z_value * se[j]
  }

  h_df <- as.data.frame(matrix(rep(h_used, each = m), nrow = m))
  names(h_df) <- h_names

  results <- cbind(
    a_eval_df,
    data.frame(
      estimate = estimate,
      se = se,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      effective_n = effective_n,
      min_gps = rep(min_gps, m)
    ),
    h_df
  )
  rownames(results) <- NULL

  list(
    results = results,
    h = h_used,
    diagnostics = list(
      level = level,
      gps_floor = gps_floor,
      normalize_ipw = normalize_ipw,
      min_gps = min_gps,
      min_effective_n = min(effective_n, na.rm = TRUE),
      max_effective_n = max(effective_n, na.rm = TRUE)
    )
  )
}

#' @export
print.ers_fit <- function(x, ...) {
  cat("Causal exposure-response surface fit\n")
  cat("Estimator:", x$estimator, "\n")
  cat("Folds:", x$L, "\n")
  cat("Evaluation points:", nrow(x$results), "\n\n")
  print(utils::head(x$results), row.names = FALSE)
  invisible(x)
}

#' @export
summary.ers_fit <- function(object, ...) {
  out <- list(
    estimator = object$estimator,
    L = object$L,
    n_eval = nrow(object$results),
    h = object$h,
    diagnostics = object$diagnostics,
    results = object$results
  )
  class(out) <- "summary.ers_fit"
  out
}

#' @export
print.summary.ers_fit <- function(x, ...) {
  cat("Summary of causal exposure-response surface fit\n")
  cat("Estimator:", x$estimator, "\n")
  cat("Folds:", x$L, "\n")
  cat("Evaluation points:", x$n_eval, "\n")
  cat("Bandwidth:", paste(signif(x$h, 4), collapse = ", "), "\n\n")
  print(utils::head(x$results), row.names = FALSE)
  invisible(x)
}

normalize_ers_inputs <- function(Y, A, C, a_eval) {
  Y <- as.numeric(Y)
  n <- length(Y)
  if (n < 1L || anyNA(Y)) {
    stop("'Y' must be a non-missing numeric vector.", call. = FALSE)
  }

  A_df <- normalize_ers_table(A, n = n, arg = "A", prefix = "A")
  C_df <- normalize_ers_table(C, n = n, arg = "C", prefix = "C")
  a_eval_df <- normalize_ers_eval(a_eval, A_df)

  list(Y = Y, A = A_df, C = C_df, a_eval = a_eval_df)
}

validate_ers_options <- function(L, gps_floor, normalize_ipw) {
  if (!is.numeric(L) || length(L) != 1L || is.na(L) || L < 1L) {
    stop("'L' must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(gps_floor) || length(gps_floor) != 1L ||
      is.na(gps_floor) || !is.finite(gps_floor) || gps_floor <= 0) {
    stop("'gps_floor' must be a positive finite number.", call. = FALSE)
  }
  if (!is.logical(normalize_ipw) || length(normalize_ipw) != 1L ||
      is.na(normalize_ipw)) {
    stop("'normalize_ipw' must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(TRUE)
}

normalize_ers_table <- function(x, n, arg, prefix) {
  if (is.null(x)) {
    stop(sprintf("'%s' cannot be NULL.", arg), call. = FALSE)
  }
  orig_names <- if (is.matrix(x) || is.data.frame(x)) colnames(x) else NULL
  if (is.null(dim(x))) {
    x <- matrix(x, ncol = 1L)
  }
  df <- as.data.frame(x)
  if (nrow(df) != n) {
    stop(sprintf("'%s' must have the same number of rows as 'Y'.", arg), call. = FALSE)
  }
  if (ncol(df) < 1L) {
    stop(sprintf("'%s' must have at least one column.", arg), call. = FALSE)
  }
  if (is.null(orig_names)) {
    colnames(df) <- paste0(prefix, seq_len(ncol(df)))
  } else {
    colnames(df) <- make.names(orig_names, unique = TRUE)
  }
  df[] <- lapply(df, as.numeric)
  if (anyNA(df)) {
    stop(sprintf("'%s' must contain only non-missing numeric values.", arg), call. = FALSE)
  }
  df
}

normalize_ers_eval <- function(a_eval, A_df) {
  p <- ncol(A_df)
  if (is.null(a_eval)) {
    out <- A_df
  } else {
    if (is.null(dim(a_eval))) {
      if (p == 1L) {
        a_eval <- matrix(a_eval, ncol = 1L)
      } else if (length(a_eval) == p) {
        a_eval <- matrix(a_eval, nrow = 1L)
      } else {
        stop("'a_eval' must have one column per exposure.", call. = FALSE)
      }
    }
    out <- as.data.frame(a_eval)
    if (ncol(out) != p) {
      stop("'a_eval' must have the same number of columns as 'A'.", call. = FALSE)
    }
    colnames(out) <- colnames(A_df)
    out[] <- lapply(out, as.numeric)
    if (anyNA(out)) {
      stop("'a_eval' must contain only non-missing numeric values.", call. = FALSE)
    }
  }
  out
}

make_observed_newdata <- function(model, A, C) {
  A_df <- as.data.frame(A)
  C_df <- as.data.frame(C)
  if (ncol(A_df) != length(model$X_names) || ncol(C_df) != length(model$C_names)) {
    stop("New data dimensions do not match the fitted nuisance model.", call. = FALSE)
  }
  colnames(A_df) <- model$X_names
  colnames(C_df) <- model$C_names
  cbind(A_df, C_df)
}

predict_outcome_eval_matrix <- function(out_model, C, a_eval) {
  C_df <- as.data.frame(C)
  a_eval_df <- as.data.frame(a_eval)
  colnames(C_df) <- out_model$C_names
  colnames(a_eval_df) <- out_model$X_names

  n <- nrow(C_df)
  m <- nrow(a_eval_df)
  x_names <- out_model$X_names
  out <- matrix(NA_real_, nrow = n, ncol = m)

  dt_grid <- data.table::as.data.table(C_df)
  for (x_col in x_names) {
    data.table::set(dt_grid, j = x_col, value = 0.0)
  }
  data.table::setcolorder(dt_grid, c(x_names, out_model$C_names))

  for (j in seq_len(m)) {
    a_target <- as.numeric(a_eval_df[j, ])
    for (k in seq_along(x_names)) {
      data.table::set(dt_grid, j = x_names[k], value = a_target[k])
    }
    out[, j] <- stats::predict(out_model, newdata = dt_grid)
  }

  out
}

resolve_ers_bandwidth <- function(A, h, c_multiplier, require_positive) {
  p <- ncol(A)
  if (is.null(h)) {
    h <- c_multiplier * apply(A, 2, stats::sd) * (nrow(A)^(-0.2))
    if (any(!is.finite(h) | h <= 0)) {
      stop(
        "Default bandwidth is nonpositive or nonfinite; check for zero-variance exposure columns.",
        call. = FALSE
      )
    }
  }
  h <- as.numeric(h)
  if (length(h) == 1L && p > 1L) {
    h <- rep(h, p)
  }
  if (length(h) != p) {
    stop("'h' must have length 1 or ncol(A).", call. = FALSE)
  }
  if (any(!is.finite(h) | h <= 0)) {
    stop("'h' must contain positive finite bandwidths.", call. = FALSE)
  }
  h
}

product_gaussian_kernel <- function(A, a, h) {
  z <- sweep(A, 2, a, FUN = "-")
  z <- sweep(z, 2, h, FUN = "/")
  apply(stats::dnorm(z), 1, prod)
}
