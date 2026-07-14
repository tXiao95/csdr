#' Causal Sufficient Dimension Reduction
#'
#' `csdr()` is the primary user-facing function for estimating a causal
#' sufficient dimension reduction basis. It first constructs generated
#' response/exposure targets with [csdr_target()], then fits MAVE separately for
#' each requested target variant.
#'
#' @details
#' The available target variants are:
#'
#' - `"RA"`: regression-adjustment target using an outcome regression.
#' - `"DR"`: doubly robust target using an outcome regression and generalized
#'   propensity score (GPS).
#' - `"PO"`: pseudo-outcome target using both nuisance components and an
#'   additional marginalization step.
#' - `"RP"`: residualized-pair target using regressions of the outcome and each
#'   exposure on the covariates.
#'
#' `"DR"` is the default. `"PO"` is generally the most computationally
#' intensive because of its marginalization step. See [csdr_target()] for the
#' lower-level target-construction interface.
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
#' @param learners Learner specifications created by [csdr_learners()].
#'   `learners$rp_y` controls `E[Y | C]` for RP targets and `learners$rp_a`
#'   controls `E[A_j | C]` for RP targets.
#' @param target_control Optional named list passed to [csdr_target()]. Supported
#'   entries include `folds`, `args_outcome`, `args_gps`, `args_rp_y`,
#'   `args_rp_a`, `args_C`, `args_ers`, and `po_marginalization`.
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
#'
#' fit_glm <- csdr(
#'   Y, A, C,
#'   learners = csdr_learners(sl_library = c("SL.glm", "SL.ranger", "SL.xgboost"))
#' )
#'
#' fit_outcome <- csdr(
#'   Y, A, C,
#'   learners = csdr_learners(
#'     outcome = sl_regression(SL.library = c("SL.glm", "SL.ranger"))
#'   )
#' )
#'
#' fit_gps <- csdr(
#'   Y, A, C,
#'   learners = csdr_learners(
#'     gps = mvn_gps(method = "SuperLearner", SL.library = c("SL.glm", "SL.ranger"))
#'   )
#' )
#'
#' my_fitter <- function(Y, W, ...) {
#'   fit <- stats::lm(Y ~ ., data = data.frame(Y = Y, W))
#'   structure(list(fit = fit), class = "my_regression")
#' }
#' predict.my_regression <- function(object, newdata, ...) {
#'   as.numeric(stats::predict(object$fit, newdata = as.data.frame(newdata)))
#' }
#' fit_custom <- csdr(
#'   Y, A, C,
#'   learners = csdr_learners(
#'     outcome = custom_regression(my_fitter, label = "custom lm")
#'   )
#' )
#'
#' # Custom GPS learners should return fitted density objects whose predict()
#' # method returns numeric conditional density estimates.
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
  learners = csdr_learners(),
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
    A = A,
    C = C,
    d = d,
    max_dim = max_dim,
    L = L,
    learners = learners,
    target_control = target_control,
    mave_control = mave_control,
    keep_targets = keep_targets,
    keep_mave = keep_mave,
    keep_nuisance = keep_nuisance,
    verbose = verbose
  )
  validate_csdr_learners(learners)

  data <- normalize_ers_inputs(Y = Y, A = A, C = C, a_eval = A)
  A_mat <- as.matrix(data$A)
  C_mat <- as.matrix(data$C)

  if (!is.null(d) && d > ncol(A_mat)) {
    stop("'d' cannot exceed the number of exposure columns.", call. = FALSE)
  }
  if (!is.null(max_dim) && max_dim > ncol(A_mat)) {
    stop("'max_dim' cannot exceed the number of exposure columns.", call. = FALSE)
  }

  learner_summary <- summarize_csdr_learners(learners = learners, variants = variants)
  if (isTRUE(verbose)) {
    announce_csdr_learners(learner_summary)
  }
  learner_args <- as_csdr_target_args(learners = learners, target_control = target_control)

  target_obj <- csdr_target(
    Y = data$Y,
    A = data$A,
    C = data$C,
    methods = variants,
    L = L,
    folds = target_control$folds,
    outcome_fitter = learner_args$outcome_fitter,
    gps_fitter = learner_args$gps_fitter,
    rp_y_fitter = learner_args$rp_y_fitter,
    rp_a_fitter = learner_args$rp_a_fitter,
    args_outcome = learner_args$args_outcome,
    args_gps = learner_args$args_gps,
    args_rp_y = learner_args$args_rp_y,
    args_rp_a = learner_args$args_rp_a,
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
    learners = learners,
    learner_summary = learner_summary,
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

as_csdr_target_args <- function(learners, target_control) {
  list(
    outcome_fitter = learners$outcome$fitter,
    gps_fitter = learners$gps$fitter,
    rp_y_fitter = learners$rp_y$fitter,
    rp_a_fitter = learners$rp_a$fitter,
    args_outcome = utils::modifyList(
      learners$outcome$args,
      target_control$args_outcome %||% list()
    ),
    args_gps = utils::modifyList(
      learners$gps$args,
      target_control$args_gps %||% list()
    ),
    args_rp_y = utils::modifyList(
      learners$rp_y$args,
      target_control$args_rp_y %||% target_control$args_C %||% list()
    ),
    args_rp_a = utils::modifyList(
      learners$rp_a$args,
      target_control$args_rp_a %||% target_control$args_C %||% list()
    )
  )
}

summarize_csdr_learners <- function(learners, variants) {
  used <- list(
    outcome = any(c("RA", "DR", "PO") %in% variants),
    gps = any(c("DR", "PO") %in% variants),
    rp_y = "RP" %in% variants,
    rp_a = "RP" %in% variants
  )

  out <- do.call(
    rbind,
    lapply(names(used), function(role) {
      learner <- learners[[role]]
      details <- describe_csdr_learner(learner)
      data.frame(
        role = role,
        used = used[[role]],
        engine = learner$engine,
        label = learner$label,
        details = details,
        is_default = isTRUE(learner$is_default),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  out
}

describe_csdr_learner <- function(learner) {
  if (identical(learner$engine, "SuperLearner")) {
    return(paste(
      "SL.library:",
      paste(learner$summary$SL.library, collapse = ", ")
    ))
  }
  if (identical(learner$engine, "MVN GPS")) {
    return(sprintf(
      "method_gps: %s; density_floor: %s",
      learner$summary$method_gps,
      format(learner$summary$density_floor, scientific = TRUE)
    ))
  }
  arg_names <- learner$summary$args
  if (length(arg_names) > 0L) {
    return(paste("args:", paste(arg_names, collapse = ", ")))
  }
  "no stored args"
}

announce_csdr_learners <- function(learner_summary) {
  message("CSDR learner choices:")
  for (i in seq_len(nrow(learner_summary))) {
    role <- learner_summary$role[[i]]
    label <- switch(
      role,
      outcome = "Outcome regression E[Y | A, C]",
      gps = "GPS f(A | C)",
      rp_y = "RP regression E[Y | C]",
      rp_a = "RP regression E[A_j | C]",
      role
    )
    if (!learner_summary$used[[i]]) {
      message(sprintf("- %s: not used", label))
    } else {
      message(sprintf("- %s: %s", label, learner_summary$label[[i]]))
      message(sprintf("  %s", learner_summary$details[[i]]))
    }
  }
  invisible(learner_summary)
}

validate_csdr_options <- function(Y, A, C, d, max_dim, L, learners, target_control,
                                  mave_control, keep_targets, keep_mave,
                                  keep_nuisance, verbose) {
  validate_csdr_data(Y = Y, A = A, C = C)
  n <- length(Y)
  if (!is.numeric(L) || length(L) != 1L || is.na(L) || !is.finite(L) ||
      L < 1L || L != as.integer(L)) {
    stop("'L' must be a positive integer.", call. = FALSE)
  }
  if (L > n) {
    stop(sprintf("'L' (%d) cannot exceed the number of observations (%d).", L, n),
         call. = FALSE)
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
  if (!inherits(learners, "csdr_learners")) {
    stop("'learners' must be created by csdr_learners().", call. = FALSE)
  }
  if (!is.list(target_control)) {
    stop("'target_control' must be a list.", call. = FALSE)
  }
  if (!is.list(mave_control)) {
    stop("'mave_control' must be a list.", call. = FALSE)
  }
  validate_control_names(
    target_control,
    valid = c(
      "folds", "args_outcome", "args_gps", "args_rp_y", "args_rp_a",
      "args_C", "args_ers", "po_marginalization"
    ),
    arg = "target_control"
  )
  reserved_mave <- intersect(names(mave_control), c("formula", "data"))
  if (length(reserved_mave) > 0L) {
    stop(
      sprintf(
        "'mave_control' cannot override internally managed argument%s: %s.",
        if (length(reserved_mave) == 1L) "" else "s",
        paste(reserved_mave, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  mave_names <- setdiff(names(formals(MAVE::mave)), c("formula", "data"))
  validate_control_names(mave_control, valid = mave_names, arg = "mave_control")
  for (arg in c("keep_targets", "keep_mave", "keep_nuisance", "verbose")) {
    value <- get(arg)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(sprintf("'%s' must be TRUE or FALSE.", arg), call. = FALSE)
    }
  }
  invisible(TRUE)
}

validate_csdr_data <- function(Y, A, C) {
  if (!is.numeric(Y) || !is.null(dim(Y)) || length(Y) < 1L ||
      any(!is.finite(Y))) {
    stop("'Y' must be a non-missing finite numeric vector.", call. = FALSE)
  }
  validate_csdr_table(A, n = length(Y), arg = "A")
  validate_csdr_table(C, n = length(Y), arg = "C")

  overlap <- intersect(colnames(A), colnames(C))
  if (length(overlap) > 0L) {
    stop(
      sprintf(
        "'A' and 'C' must have distinct column names; duplicated across inputs: %s.",
        paste(overlap, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

validate_csdr_table <- function(x, n, arg) {
  if (!(is.matrix(x) || is.data.frame(x))) {
    stop(sprintf("'%s' must be a numeric matrix or data frame.", arg), call. = FALSE)
  }
  if (nrow(x) != n) {
    stop(
      sprintf("'%s' must have %d rows to match 'Y'; received %d.", arg, n, nrow(x)),
      call. = FALSE
    )
  }
  if (ncol(x) < 1L) {
    stop(sprintf("'%s' must have at least one column.", arg), call. = FALSE)
  }

  numeric_columns <- vapply(as.data.frame(x), is.numeric, logical(1))
  if (!all(numeric_columns)) {
    stop(
      sprintf(
        "'%s' contains nonnumeric column%s: %s.",
        arg,
        if (sum(!numeric_columns) == 1L) "" else "s",
        paste(names(numeric_columns)[!numeric_columns], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  column_names <- colnames(x)
  if (is.null(column_names) || anyNA(column_names) || any(!nzchar(column_names))) {
    stop(sprintf("'%s' must have nonempty column names.", arg), call. = FALSE)
  }
  duplicated_names <- unique(column_names[duplicated(column_names)])
  if (length(duplicated_names) > 0L) {
    stop(
      sprintf("'%s' has duplicated column names: %s.",
              arg, paste(duplicated_names, collapse = ", ")),
      call. = FALSE
    )
  }

  values <- as.matrix(x)
  if (any(!is.finite(values))) {
    bad_columns <- column_names[colSums(!is.finite(values)) > 0L]
    stop(
      sprintf(
        "'%s' contains missing or non-finite values in column%s: %s.",
        arg,
        if (length(bad_columns) == 1L) "" else "s",
        paste(bad_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  constant <- vapply(as.data.frame(x), function(column) {
    length(unique(column)) < 2L
  }, logical(1))
  if (any(constant)) {
    stop(
      sprintf(
        "'%s' contains constant column%s: %s.",
        arg,
        if (sum(constant) == 1L) "" else "s",
        paste(column_names[constant], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

validate_control_names <- function(control, valid, arg) {
  supplied <- names(control)
  if (length(control) > 0L &&
      (is.null(supplied) || anyNA(supplied) || any(!nzchar(supplied)))) {
    stop(sprintf("'%s' must be a fully named list.", arg), call. = FALSE)
  }
  if (anyDuplicated(supplied)) {
    duplicated_names <- unique(supplied[duplicated(supplied)])
    stop(sprintf("'%s' has duplicated entries: %s.",
                 arg, paste(duplicated_names, collapse = ", ")),
         call. = FALSE)
  }
  unknown <- setdiff(supplied, valid)
  if (length(unknown) > 0L) {
    suggestions <- vapply(unknown, function(entry) {
      distance <- utils::adist(entry, valid)
      closest <- valid[which.min(distance)]
      if (min(distance) <= max(2L, floor(nchar(entry) / 3L))) {
        sprintf(" Did you mean '%s'?", closest)
      } else {
        ""
      }
    }, character(1))
    stop(
      paste0(
        "Unknown entry in '", arg, "': '", unknown[[1L]], "'.",
        suggestions[[1L]]
      ),
      call. = FALSE
    )
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
    target_diagnostics = object$target$diagnostics,
    learner_summary = object$learner_summary
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
