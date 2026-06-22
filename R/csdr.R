#' Causal Sufficient Dimension Reduction
#'
#' `csdr()` is the primary user-facing function for estimating a causal
#' sufficient dimension reduction basis. It first constructs generated
#' response/exposure targets with [csdr_target()], then fits MAVE separately for
#' each requested target variant.
#'
#' @param Y Numeric outcome vector of length `n`.
#' @param A Numeric matrix or data frame of observed continuous exposures.
#' @param C Numeric matrix or data frame of observed covariates.
#' @param variants Target variants to fit. One or more of `"DR"`, `"RA"`,
#'   `"PO"`, and `"RP"`.
#' @param d Optional structural dimension. If `NULL`, [MAVE::mave.dim()] is
#'   used to select the dimension for each variant.
#' @param max_dim Optional maximum dimension passed to [MAVE::mave.dim()].
#'   Defaults to the number of exposure columns for each target.
#' @param L Number of folds for target construction. `L = 1` disables
#'   cross-fitting.
#' @param seed Optional random seed used when folds are generated.
#' @param fitters Optional named list of nuisance fitters. Supported names are
#'   `outcome`, `gps`, and `C`.
#' @param target_control Optional named list passed to [csdr_target()]. Supported
#'   entries include `folds`, `args_outcome`, `args_gps`, `args_C`, `args_ers`,
#'   and `po_marginalization`.
#' @param mave_control Optional named list of arguments passed to
#'   [MAVE::mave()].
#' @param keep_targets Logical; if `TRUE`, keep generated target responses and
#'   exposures in each variant fit.
#' @param keep_mave Logical; if `TRUE`, keep raw MAVE fit and dimension-selection
#'   objects.
#' @param keep_nuisance Logical; if `TRUE`, keep fitted nuisance objects inside
#'   the returned target object.
#' @param verbose Logical; if `TRUE`, print progress messages.
#'
#' @return An S3 object of class `"csdr_fit"`.
#'
#' @examples
#' \dontrun{
#' fit <- csdr(Y, A, C)
#' fit_all <- csdr(Y, A, C, variants = c("RA", "DR", "PO", "RP"))
#' coef(fit, variant = "DR")
#' }
#'
#' @export
csdr <- function(
  Y,
  A,
  C,
  variants = "DR",
  d = NULL,
  max_dim = NULL,
  L = 5,
  seed = NULL,
  fitters = list(),
  target_control = list(),
  mave_control = list(),
  keep_targets = TRUE,
  keep_mave = TRUE,
  keep_nuisance = FALSE,
  verbose = TRUE
) {
  call <- match.call()
  valid_variants <- c("RA", "DR", "PO", "RP")
  variants <- match.arg(variants, choices = valid_variants, several.ok = TRUE)

  validate_csdr_options(
    Y = Y,
    d = d,
    max_dim = max_dim,
    L = L,
    fitters = fitters,
    target_control = target_control,
    mave_control = mave_control,
    keep_targets = keep_targets,
    keep_mave = keep_mave,
    keep_nuisance = keep_nuisance,
    verbose = verbose
  )

  data <- normalize_ers_inputs(Y = Y, A = A, C = C, a_eval = A)
  A_mat <- as.matrix(data$A)
  C_mat <- as.matrix(data$C)

  if (!is.null(d) && d > ncol(A_mat)) {
    stop("'d' cannot exceed the number of exposure columns.", call. = FALSE)
  }
  if (!is.null(max_dim) && max_dim > ncol(A_mat)) {
    stop("'max_dim' cannot exceed the number of exposure columns.", call. = FALSE)
  }

  target_obj <- csdr_target(
    Y = data$Y,
    A = data$A,
    C = data$C,
    methods = variants,
    L = L,
    folds = target_control$folds,
    outcome_fitter = fitters$outcome %||% SL_outcome_fitter,
    gps_fitter = fitters$gps %||% mvn_fitter,
    C_fitter = fitters$C %||% SL_nuisance_fitter,
    args_outcome = target_control$args_outcome %||% list(),
    args_gps = target_control$args_gps %||% list(),
    args_C = target_control$args_C %||% list(),
    args_ers = target_control$args_ers %||% list(),
    po_marginalization = target_control$po_marginalization %||% "crossfit",
    seed = seed,
    return_nuisance = keep_nuisance,
    verbose = verbose
  )

  fits <- lapply(variants, function(variant) {
    fit_csdr_variant(
      variant = variant,
      target_Y = target_obj$new_Y[[variant]],
      target_A = target_obj$new_A[[variant]],
      d = d,
      max_dim = max_dim,
      mave_control = mave_control,
      keep_targets = keep_targets,
      keep_mave = keep_mave,
      target_diagnostics = target_obj$diagnostics,
      verbose = verbose
    )
  })
  names(fits) <- variants

  summary_table <- do.call(
    rbind,
    lapply(fits, function(x) {
      data.frame(
        variant = x$variant,
        d_hat = x$d_hat,
        n = x$diagnostics$n,
        p = x$diagnostics$p,
        mave_method = x$diagnostics$mave_method,
        target_exposure = x$diagnostics$target_exposure,
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(summary_table) <- NULL

  out <- list(
    call = call,
    variants = variants,
    input = list(
      n = length(data$Y),
      p = ncol(A_mat),
      q = ncol(C_mat),
      A_names = colnames(A_mat),
      C_names = colnames(C_mat)
    ),
    target = target_obj,
    fits = fits,
    summary = summary_table,
    control = list(
      L = L,
      seed = seed,
      d = d,
      max_dim = max_dim,
      keep_targets = keep_targets,
      keep_mave = keep_mave,
      keep_nuisance = keep_nuisance,
      target_control = target_control,
      mave_control = mave_control
    )
  )
  class(out) <- "csdr_fit"
  out
}

fit_csdr_variant <- function(variant, target_Y, target_A, d, max_dim,
                             mave_control, keep_targets, keep_mave,
                             target_diagnostics, verbose) {
  target_A <- as.matrix(target_A)
  target_Y <- as.numeric(target_Y)
  if (length(target_Y) != nrow(target_A)) {
    stop("Target response and exposure have incompatible dimensions.", call. = FALSE)
  }
  if (anyNA(target_Y) || anyNA(target_A)) {
    stop("Target response and exposure cannot contain missing values.", call. = FALSE)
  }

  max_dim_variant <- max_dim %||% ncol(target_A)
  if (max_dim_variant > ncol(target_A)) {
    stop("'max_dim' cannot exceed the number of target exposure columns.", call. = FALSE)
  }
  if (!is.null(d) && d > ncol(target_A)) {
    stop("'d' cannot exceed the number of target exposure columns.", call. = FALSE)
  }

  if (isTRUE(verbose)) {
    message(sprintf("Running MAVE for %s target.", variant))
  }

  df <- data.frame(target_Y = target_Y, as.data.frame(target_A))
  mave_args <- utils::modifyList(
    list(formula = target_Y ~ ., data = df, method = "meanMAVE"),
    mave_control
  )
  mave_fit <- do.call(MAVE::mave, mave_args)

  dim_obj <- NULL
  if (is.null(d)) {
    dim_obj <- MAVE::mave.dim(mave_fit, max.dim = max_dim_variant)
    d_hat <- extract_mave_dim(dim_obj)
  } else {
    d_hat <- as.integer(d)
  }

  beta <- extract_mave_beta(mave_fit, d_hat = d_hat)
  score <- target_A %*% beta
  colnames(score) <- paste0("score", seq_len(ncol(score)))

  list(
    variant = variant,
    beta = beta,
    d_hat = d_hat,
    score = score,
    target_Y = if (keep_targets) target_Y else NULL,
    target_A = if (keep_targets) target_A else NULL,
    mave_fit = if (keep_mave) mave_fit else NULL,
    mave_dim_obj = if (keep_mave) dim_obj else NULL,
    diagnostics = list(
      n = nrow(target_A),
      p = ncol(target_A),
      mave_method = mave_args$method,
      target_exposure = if (variant == "RP") "residualized_A" else "A",
      target_diagnostics = target_diagnostics
    )
  )
}

extract_mave_dim <- function(dim_obj) {
  d_hat <- dim_obj$dim.min %||% dim_obj$dim %||% NULL
  if (is.null(d_hat) || length(d_hat) != 1L || is.na(d_hat)) {
    stop("Could not extract selected dimension from the MAVE dimension object.", call. = FALSE)
  }
  d_hat <- as.integer(d_hat)
  if (d_hat < 1L) {
    stop("MAVE selected a nonpositive dimension.", call. = FALSE)
  }
  d_hat
}

extract_mave_beta <- function(mave_fit, d_hat) {
  if (!is.numeric(d_hat) || length(d_hat) != 1L || is.na(d_hat) || d_hat < 1L) {
    stop("'d_hat' must be a positive integer.", call. = FALSE)
  }
  d_hat <- as.integer(d_hat)
  p <- if (!is.null(mave_fit$x)) ncol(as.matrix(mave_fit$x)) else NULL

  candidates <- list()
  if (!is.null(mave_fit$dir)) {
    if (is.list(mave_fit$dir)) {
      if (length(mave_fit$dir) >= d_hat) {
        candidates[[length(candidates) + 1L]] <- mave_fit$dir[[d_hat]]
      }
      candidates <- c(candidates, mave_fit$dir)
    } else {
      candidates[[length(candidates) + 1L]] <- mave_fit$dir
    }
  }
  for (field in c("beta", "basis", "B", "directions")) {
    if (!is.null(mave_fit[[field]])) {
      candidates[[length(candidates) + 1L]] <- mave_fit[[field]]
    }
  }

  for (candidate in candidates) {
    beta <- coerce_mave_beta(candidate, d_hat = d_hat, p = p)
    if (!is.null(beta)) {
      colnames(beta) <- paste0("dir", seq_len(ncol(beta)))
      return(beta)
    }
  }

  stop(
    "Could not extract a beta matrix from the MAVE object. ",
    "Expected a direction matrix with one row per exposure and at least 'd_hat' columns.",
    call. = FALSE
  )
}

coerce_mave_beta <- function(candidate, d_hat, p) {
  if (is.null(candidate)) {
    return(NULL)
  }
  mat <- tryCatch(as.matrix(candidate), error = function(e) NULL)
  if (is.null(mat) || !is.numeric(mat) || length(dim(mat)) != 2L) {
    return(NULL)
  }

  if (!is.null(p) && nrow(mat) == p && ncol(mat) >= d_hat) {
    return(mat[, seq_len(d_hat), drop = FALSE])
  }
  if (!is.null(p) && ncol(mat) == p && nrow(mat) >= d_hat) {
    return(t(mat[seq_len(d_hat), , drop = FALSE]))
  }
  if (is.null(p) && ncol(mat) >= d_hat) {
    return(mat[, seq_len(d_hat), drop = FALSE])
  }
  if (is.null(p) && nrow(mat) >= d_hat) {
    return(t(mat[seq_len(d_hat), , drop = FALSE]))
  }
  NULL
}

validate_csdr_options <- function(Y, d, max_dim, L, fitters, target_control,
                                  mave_control, keep_targets, keep_mave,
                                  keep_nuisance, verbose) {
  if (!is.numeric(Y) || length(Y) < 1L || anyNA(Y)) {
    stop("'Y' must be a non-missing numeric vector.", call. = FALSE)
  }
  if (!is.numeric(L) || length(L) != 1L || is.na(L) || L < 1L) {
    stop("'L' must be a positive integer.", call. = FALSE)
  }
  if (!is.null(d) &&
      (!is.numeric(d) || length(d) != 1L || is.na(d) || d < 1L || d != as.integer(d))) {
    stop("'d' must be a positive integer when supplied.", call. = FALSE)
  }
  if (!is.null(max_dim) &&
      (!is.numeric(max_dim) || length(max_dim) != 1L || is.na(max_dim) ||
       max_dim < 1L || max_dim != as.integer(max_dim))) {
    stop("'max_dim' must be a positive integer when supplied.", call. = FALSE)
  }
  if (!is.list(fitters)) {
    stop("'fitters' must be a list.", call. = FALSE)
  }
  if (!is.list(target_control)) {
    stop("'target_control' must be a list.", call. = FALSE)
  }
  if (!is.list(mave_control)) {
    stop("'mave_control' must be a list.", call. = FALSE)
  }
  for (arg in c("keep_targets", "keep_mave", "keep_nuisance", "verbose")) {
    value <- get(arg)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(sprintf("'%s' must be TRUE or FALSE.", arg), call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' @export
print.csdr_fit <- function(x, ...) {
  cat("Causal sufficient dimension reduction fit\n")
  cat("Variants:", paste(x$variants, collapse = ", "), "\n")
  cat("Observations:", x$input$n, "\n")
  cat("Exposures:", x$input$p, "\n")
  cat("Folds:", x$control$L, "\n\n")
  print(x$summary, row.names = FALSE)
  invisible(x)
}

#' @export
summary.csdr_fit <- function(object, ...) {
  out <- list(
    summary = object$summary,
    input = object$input,
    variants = object$variants,
    target_diagnostics = object$target$diagnostics
  )
  class(out) <- "summary.csdr_fit"
  out
}

#' @export
print.summary.csdr_fit <- function(x, ...) {
  cat("Summary of causal sufficient dimension reduction fit\n\n")
  print(x$summary, row.names = FALSE)
  invisible(x)
}

#' @export
coef.csdr_fit <- function(object, variant = NULL, ...) {
  selected <- select_csdr_variants(object, variant)
  out <- lapply(selected, function(v) object$fits[[v]]$beta)
  if (length(out) == 1L) out[[1L]] else out
}

#' Extract CSDR scores
#'
#' @param object A fitted object.
#' @param ... Additional arguments passed to methods.
#'
#' @export
scores <- function(object, ...) {
  UseMethod("scores")
}

#' @param variant Optional variant name. If `NULL`, return all variants.
#' @rdname scores
#' @export
scores.csdr_fit <- function(object, variant = NULL, ...) {
  selected <- select_csdr_variants(object, variant)
  out <- lapply(selected, function(v) object$fits[[v]]$score)
  if (length(out) == 1L) out[[1L]] else out
}

#' Extract generated CSDR targets
#'
#' @param object A fitted object.
#' @param ... Additional arguments passed to methods.
#'
#' @export
targets <- function(object, ...) {
  UseMethod("targets")
}

#' @param variant Optional variant name. If `NULL`, return all variants.
#' @rdname targets
#' @export
targets.csdr_fit <- function(object, variant = NULL, ...) {
  selected <- select_csdr_variants(object, variant)
  out <- lapply(selected, function(v) {
    fit <- object$fits[[v]]
    if (is.null(fit$target_Y) || is.null(fit$target_A)) {
      stop("Targets were not retained. Refit with 'keep_targets = TRUE'.", call. = FALSE)
    }
    list(target_Y = fit$target_Y, target_A = fit$target_A)
  })
  if (length(out) == 1L) out[[1L]] else out
}

select_csdr_variants <- function(object, variant) {
  if (is.null(variant)) {
    return(object$variants)
  }
  variant <- match.arg(variant, choices = object$variants, several.ok = TRUE)
  variant
}
