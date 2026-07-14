#' Construct CSDR target response and exposure pairs
#'
#' `csdr_target()` constructs generated response/exposure pairs for causal SDR
#' methods while sharing nuisance fits across compatible targets.
#'
#' @param Y Numeric outcome vector.
#' @param A Numeric matrix or data frame of observed exposures.
#' @param C Numeric matrix or data frame of covariates.
#' @param methods Methods to construct: `"RA"`, `"DR"`, `"PO"`, and/or `"RP"`.
#' @param L Number of folds. `L = 1` disables cross-fitting.
#' @param folds Optional fold assignments.
#' @param outcome_fitter Function used by [outcome_model()] for RA, DR, and PO.
#' @param gps_fitter Function used by [gps_model()] for DR and PO.
#' @param rp_y_fitter Function used by [fit_rp_nuisance()] for `E[Y | C]`
#'   in RP.
#' @param rp_a_fitter Function used by [fit_rp_nuisance()] for `E[A_j | C]`
#'   in RP.
#' @param C_fitter Compatibility alias for using one RP fitter for both
#'   `rp_y_fitter` and `rp_a_fitter`.
#' @param args_outcome Additional arguments passed to `outcome_fitter`.
#' @param args_gps Additional arguments passed to `gps_fitter`.
#' @param args_rp_y Additional arguments passed to `rp_y_fitter`.
#' @param args_rp_a Additional arguments passed to `rp_a_fitter`.
#' @param args_C Compatibility alias for using one RP argument list for both
#'   `args_rp_y` and `args_rp_a`.
#' @param args_ers Additional arguments passed to ERS target construction for
#'   RA and DR. Supported entries include `h`, `c_multiplier`, `gps_floor`,
#'   `delta_n` as a compatibility alias for `gps_floor`, `normalize_ipw`, and
#'   `level`.
#' @param po_marginalization Marginalization mode for PO. See
#'   [estimate_pseudo_outcomes()].
#' @param seed Optional random seed for fold generation.
#' @param return_nuisance Logical; if `TRUE`, return fitted nuisance objects.
#' @param verbose Logical; if `TRUE`, print progress messages.
#'
#' @return An S3 object of class `"csdr_target"`.
#'
#' @export
csdr_target <- function(
  Y,
  A,
  C,
  methods = c("DR", "RA", "PO", "RP"),
  L = 5,
  folds = NULL,
  outcome_fitter = SL_outcome_fitter,
  gps_fitter = mvn_fitter,
  rp_y_fitter = SL_nuisance_fitter,
  rp_a_fitter = SL_nuisance_fitter,
  C_fitter = NULL,
  args_outcome = list(),
  args_gps = list(),
  args_rp_y = list(),
  args_rp_a = list(),
  args_C = NULL,
  args_ers = list(),
  po_marginalization = c("crossfit", "fold", "all"),
  seed = NULL,
  return_nuisance = FALSE,
  verbose = TRUE
) {
  call <- match.call()
  rp_y_fitter_missing <- missing(rp_y_fitter)
  rp_a_fitter_missing <- missing(rp_a_fitter)
  args_rp_y_missing <- missing(args_rp_y)
  args_rp_a_missing <- missing(args_rp_a)
  valid_methods <- c("RA", "DR", "PO", "RP")
  methods <- match.arg(methods, choices = valid_methods, several.ok = TRUE)
  po_marginalization <- match.arg(po_marginalization)

  if (!is.null(C_fitter)) {
    if (rp_y_fitter_missing) {
      rp_y_fitter <- C_fitter
    }
    if (rp_a_fitter_missing) {
      rp_a_fitter <- C_fitter
    }
  }
  if (!is.null(args_C)) {
    if (args_rp_y_missing || identical(args_rp_y, list())) {
      args_rp_y <- args_C
    }
    if (args_rp_a_missing || identical(args_rp_a, list())) {
      args_rp_a <- args_C
    }
  }

  data <- normalize_ers_inputs(Y = Y, A = A, C = C, a_eval = A)
  folds <- make_folds(n = length(data$Y), L = L, folds = folds, seed = seed)
  L_eff <- length(unique(folds))
  n <- length(data$Y)
  p <- ncol(data$A)
  A_mat <- as.matrix(data$A)
  C_mat <- as.matrix(data$C)

  needs_outcome <- any(c("RA", "DR", "PO") %in% methods)
  needs_gps <- any(c("DR", "PO") %in% methods)
  needs_rp <- "RP" %in% methods

  if (needs_outcome && !is.function(outcome_fitter)) {
    stop("'outcome_fitter' must be a function for RA, DR, and PO.", call. = FALSE)
  }
  if (needs_gps && !is.function(gps_fitter)) {
    stop("'gps_fitter' must be a function for DR and PO.", call. = FALSE)
  }
  if (needs_rp && !is.function(rp_y_fitter)) {
    stop("'rp_y_fitter' must be a function for RP.", call. = FALSE)
  }
  if (needs_rp && !is.function(rp_a_fitter)) {
    stop("'rp_a_fitter' must be a function for RP.", call. = FALSE)
  }

  po_gps_floor <- resolve_target_gps_floor(args_ers)

  new_Y <- lapply(methods, function(x) rep(NA_real_, n))
  names(new_Y) <- methods
  new_A <- lapply(methods, function(x) A_mat)
  names(new_A) <- methods

  nuisance_store <- list(
    outcome_models = if (needs_outcome) vector("list", L_eff) else NULL,
    gps_models = if (needs_gps) vector("list", L_eff) else NULL,
    po_models = if ("PO" %in% methods) vector("list", L_eff) else NULL,
    rp_models = if (needs_rp) vector("list", L_eff) else NULL
  )

  po_m_obs <- po_pi_obs <- if ("PO" %in% methods) rep(NA_real_, n) else NULL
  modular_ers_args <- prepare_modular_ers_args(args_ers)

  for (k in seq_len(L_eff)) {
    if (isTRUE(verbose) && L_eff > 1L) {
      message(sprintf("Constructing CSDR targets for fold %d of %d.", k, L_eff))
    }

    train_idx <- if (L_eff == 1L) seq_len(n) else which(folds != k)
    test_idx <- if (L_eff == 1L) seq_len(n) else which(folds == k)

    Y_train <- data$Y[train_idx]
    A_train <- data$A[train_idx, , drop = FALSE]
    C_train <- data$C[train_idx, , drop = FALSE]
    Y_test <- data$Y[test_idx]
    A_test <- data$A[test_idx, , drop = FALSE]
    C_test <- data$C[test_idx, , drop = FALSE]

    out_mod <- gps_mod <- NULL
    if (needs_outcome) {
      out_mod <- fit_csdr_target_outcome_nuisance(
        Y = Y_train,
        A = A_train,
        C = C_train,
        outcome_fitter = outcome_fitter,
        args_outcome = args_outcome
      )
      nuisance_store$outcome_models[[k]] <- out_mod
    }
    if (needs_gps) {
      gps_mod <- fit_csdr_target_gps_nuisance(
        A = A_train,
        C = C_train,
        gps_fitter = gps_fitter,
        args_gps = args_gps
      )
      nuisance_store$gps_models[[k]] <- gps_mod
    }

    if ("RA" %in% methods) {
      new_Y[["RA"]][test_idx] <- compute_modular_ers_target(
        Y = Y_test,
        A = A_test,
        C = C_test,
        estimator = "RA",
        out_model = out_mod,
        gps_model = NULL,
        args_ers = modular_ers_args
      )
    }
    if ("DR" %in% methods) {
      new_Y[["DR"]][test_idx] <- compute_modular_ers_target(
        Y = Y_test,
        A = A_test,
        C = C_test,
        estimator = "DR",
        out_model = out_mod,
        gps_model = gps_mod,
        args_ers = modular_ers_args
      )
    }
    if ("PO" %in% methods) {
      po_nuisance <- list(
        outcome_model = out_mod,
        gps_model = gps_mod,
        A_names = out_mod$X_names,
        C_names = out_mod$C_names
      )
      class(po_nuisance) <- "po_nuisance"
      nuisance_store$po_models[[k]] <- po_nuisance

      if (po_marginalization == "crossfit") {
        po_obs <- predict_po_observed(
          nuisance = po_nuisance,
          A = A_test,
          C = C_test,
          gps_floor = po_gps_floor
        )
        po_m_obs[test_idx] <- po_obs$m_obs
        po_pi_obs[test_idx] <- po_obs$pi_obs
      } else {
        po_pred <- predict_po_nuisance(
          nuisance = po_nuisance,
          A = A_test,
          C = C_test,
          C_marginal = if (po_marginalization == "fold") C_test else data$C,
          gps_floor = po_gps_floor,
          verbose = FALSE
        )
        po_m_obs[test_idx] <- po_pred$m_obs
        po_pi_obs[test_idx] <- po_pred$pi_obs
        new_Y[["PO"]][test_idx] <- compute_pseudo_outcomes(
          Y = Y_test,
          m_obs = po_pred$m_obs,
          pi_obs = po_pred$pi_obs,
          m_marginal = po_pred$m_marginal,
          pi_marginal = po_pred$pi_marginal,
          gps_floor = po_gps_floor
        )
      }
    }
    if (needs_rp) {
      rp_mod <- fit_rp_nuisance(
        Y = Y_train,
        A = A_train,
        C = C_train,
        y_fitter = rp_y_fitter,
        a_fitter = rp_a_fitter,
        args_y = args_rp_y,
        args_a = args_rp_a,
        verbose = FALSE
      )
      nuisance_store$rp_models[[k]] <- rp_mod
      rp_pred <- predict_rp_nuisance(nuisance = rp_mod, C = C_test)
      rp_res <- compute_residualized_pair(
        Y = Y_test,
        A = A_test,
        EY_C = rp_pred$EY_C,
        EA_C = rp_pred$EA_C
      )
      new_Y[["RP"]][test_idx] <- rp_res$Y_tilde
      new_A[["RP"]][test_idx, ] <- rp_res$A_tilde
    }
  }

  if ("PO" %in% methods && po_marginalization == "crossfit") {
    po_marginal <- predict_po_crossfit_marginals(
      nuisance_fits = list(folds = nuisance_store$po_models),
      A_targets = data$A,
      C = data$C,
      folds = folds,
      gps_floor = po_gps_floor,
      verbose = FALSE
    )
    new_Y[["PO"]] <- compute_pseudo_outcomes(
      Y = data$Y,
      m_obs = po_m_obs,
      pi_obs = po_pi_obs,
      m_marginal = po_marginal$m_marginal,
      pi_marginal = po_marginal$pi_marginal,
      gps_floor = po_gps_floor
    )
  }

  out <- list(
    new_Y = new_Y,
    new_A = lapply(new_A, as.matrix),
    folds = folds,
    methods = methods,
    nuisance = if (return_nuisance) nuisance_store else NULL,
    diagnostics = list(
      n = n,
      p = p,
      L = L_eff,
      requested_L = L,
      methods = methods,
      po_marginalization = po_marginalization
    ),
    call = call
  )
  class(out) <- "csdr_target"
  out
}

fit_csdr_target_outcome_nuisance <- function(Y, A, C, outcome_fitter, args_outcome) {
  outcome_args <- list(Y = Y, X = A, C = C, mu_fitter = outcome_fitter)
  do.call(outcome_model, c(outcome_args, args_outcome))
}

fit_csdr_target_gps_nuisance <- function(A, C, gps_fitter, args_gps) {
  gps_args <- list(X = A, C = C, pi_fitter = gps_fitter)
  do.call(gps_model, c(gps_args, args_gps))
}

compute_modular_ers_target <- function(Y, A, C, estimator, out_model, gps_model,
                                       args_ers) {
  estimator <- match.arg(estimator, choices = c("RA", "DR"))
  n <- length(Y)
  nuisance <- list(
    estimator = estimator,
    L = 1L,
    folds = rep.int(1L, n),
    outcome_models = list(out_model),
    gps_models = if (estimator == "DR") list(gps_model) else list(NULL)
  )
  class(nuisance) <- "ers_nuisance_fit"

  nuisance_pred <- predict_ers_nuisance(
    nuisance = nuisance,
    Y = Y,
    A = A,
    C = C,
    a_eval = A,
    estimator = estimator,
    gps_floor = args_ers$gps_floor,
    verbose = FALSE
  )

  compute_args <- c(
    list(
      Y = Y,
      A = A,
      a_eval = A,
      estimator = estimator,
      m_eval = nuisance_pred$m_eval,
      pi_obs = nuisance_pred$pi_obs
    ),
    args_ers
  )
  ers <- do.call(compute_ers, compute_args)
  as.numeric(ers$results$estimate)
}

prepare_modular_ers_args <- function(args_ers) {
  out <- args_ers
  if (!is.null(out$gps_floor) && is.null(out$delta_n)) {
    out$delta_n <- out$gps_floor
  }
  if (!is.null(out$delta_n) && is.null(out$gps_floor)) {
    out$gps_floor <- out$delta_n
  }
  out$delta_n <- NULL
  out$po_marginalization <- NULL
  out$marginalization <- NULL
  if (isTRUE(out$optimize_bw)) {
    warning(
      "'optimize_bw' is not implemented in modular ERS target construction and is ignored.",
      call. = FALSE
    )
  }
  out$optimize_bw <- NULL
  out$return_vector <- NULL
  supported <- c("h", "c_multiplier", "gps_floor", "normalize_ipw", "level")
  unsupported <- setdiff(names(out), supported)
  if (length(unsupported) > 0L) {
    stop(
      "Unsupported entries in 'args_ers': ",
      paste(unsupported, collapse = ", "),
      call. = FALSE
    )
  }
  if (is.null(out$gps_floor)) {
    out$gps_floor <- 1e-8
  }
  if (is.null(out$normalize_ipw)) {
    out$normalize_ipw <- TRUE
  }
  out
}

resolve_target_gps_floor <- function(args_ers) {
  if (!is.null(args_ers$gps_floor)) {
    gps_floor <- args_ers$gps_floor
  } else if (!is.null(args_ers$delta_n)) {
    gps_floor <- args_ers$delta_n
  } else {
    gps_floor <- 1e-8
  }
  if (!is.numeric(gps_floor) || length(gps_floor) != 1L ||
      is.na(gps_floor) || !is.finite(gps_floor) || gps_floor <= 0) {
    stop("'gps_floor' must be a positive finite number.", call. = FALSE)
  }
  gps_floor
}
