#' Default SuperLearner regression library for CSDR
#'
#' @return A character vector of SuperLearner wrapper names.
#'
#' @export
csdr_default_sl_library <- function() {
  c("SL.glm", "SL.glmnet", "SL.xgboost", "SL.earth")
}

#' Construct CSDR learner specifications
#'
#' @param SL.library Character vector of SuperLearner wrapper names.
#' @param family Regression family passed to [SuperLearner::SuperLearner()].
#' @param ... Additional arguments stored for the low-level fitter.
#'
#' @return A CSDR learner specification object.
#'
#' @export
sl_regression <- function(SL.library = csdr_default_sl_library(),
                          family = stats::gaussian(), ...) {
  out <- list(
    type = "regression",
    engine = "SuperLearner",
    label = "SuperLearner regression",
    fitter = sl_regression_fitter,
    args = c(list(SL.library = SL.library, family = family), list(...)),
    summary = list(
      SL.library = SL.library,
      family = family$family %||% "gaussian"
    ),
    is_default = FALSE
  )
  class(out) <- c("csdr_regression_learner", "csdr_learner")
  out
}

#' @param method Conditional mean model for the multivariate normal GPS.
#' @param density_floor Positive lower bound for predicted densities.
#'
#' @rdname sl_regression
#' @export
mvn_gps <- function(method = c("linear", "SuperLearner"),
                    density_floor = 1e-16, ...) {
  method <- match.arg(method)
  if (!is.numeric(density_floor) || length(density_floor) != 1L ||
      is.na(density_floor) || !is.finite(density_floor) || density_floor <= 0) {
    stop("'density_floor' must be a positive finite number.", call. = FALSE)
  }

  out <- list(
    type = "gps",
    engine = "MVN GPS",
    label = paste("MVN GPS with", method, "conditional mean"),
    fitter = mvn_gps_fitter,
    args = c(list(method_gps = method, delta_n = density_floor), list(...)),
    summary = list(
      method_gps = method,
      density_floor = density_floor
    ),
    is_default = FALSE
  )
  class(out) <- c("csdr_gps_learner", "csdr_learner")
  out
}

#' @param fitter User-supplied fitting function.
#' @param args Additional arguments stored for the fitter.
#' @param label Optional human-readable learner label.
#'
#' @details
#' `custom_regression()` expects `fitter(Y, W, ...)` to return a fitted object
#' with a working `predict(object, newdata = ...)` method. `custom_gps()`
#' expects `fitter(A, C, ...)` to return a fitted density object whose
#' `predict(object, newdata = ...)` method returns numeric conditional density
#' estimates.
#'
#' @rdname sl_regression
#' @export
custom_regression <- function(fitter, args = list(), label = NULL) {
  if (!is.function(fitter)) {
    stop("'fitter' must be a function.", call. = FALSE)
  }
  if (!is.list(args)) {
    stop("'args' must be a list.", call. = FALSE)
  }
  out <- list(
    type = "regression",
    engine = "custom",
    label = label %||% "custom regression",
    fitter = fitter,
    args = args,
    summary = list(args = names(args)),
    is_default = FALSE
  )
  class(out) <- c("csdr_regression_learner", "csdr_learner")
  out
}

#' @rdname sl_regression
#' @export
custom_gps <- function(fitter, args = list(), label = NULL) {
  if (!is.function(fitter)) {
    stop("'fitter' must be a function.", call. = FALSE)
  }
  if (!is.list(args)) {
    stop("'args' must be a list.", call. = FALSE)
  }
  gps_fitter <- function(A = NULL, C, X = NULL, ...) {
    if (is.null(A)) {
      A <- X
    }
    if (is.null(A)) {
      stop("'A' must be supplied.", call. = FALSE)
    }
    fitter(A = A, C = C, ...)
  }
  out <- list(
    type = "gps",
    engine = "custom",
    label = label %||% "custom GPS",
    fitter = gps_fitter,
    args = args,
    summary = list(args = names(args)),
    is_default = FALSE
  )
  class(out) <- c("csdr_gps_learner", "csdr_learner")
  out
}

#' Specify learners for CSDR nuisance fitting
#'
#' @param outcome Learner for `E[Y | A, C]`.
#' @param gps Learner for `f(A | C)`.
#' @param rp_y Learner for `E[Y | C]` in residualized-pair targets.
#' @param rp_a Learner for `E[A_j | C]` in residualized-pair targets.
#' @param sl_library Optional SuperLearner library used by unspecified
#'   SuperLearner regression learners.
#'
#' @return A list of CSDR learner specifications with class `"csdr_learners"`.
#'
#' @export
csdr_learners <- function(outcome = NULL,
                          gps = NULL,
                          rp_y = NULL,
                          rp_a = NULL,
                          sl_library = NULL) {
  default_sl <- sl_library %||% csdr_default_sl_library()

  if (is.null(outcome)) {
    outcome <- mark_default_learner(sl_regression(SL.library = default_sl))
  }
  if (is.null(gps)) {
    gps <- mark_default_learner(mvn_gps(method = "linear"))
  }
  if (is.null(rp_y)) {
    rp_y <- mark_default_learner(sl_regression(SL.library = default_sl))
  }
  if (is.null(rp_a)) {
    rp_a <- mark_default_learner(sl_regression(SL.library = default_sl))
  }

  out <- list(outcome = outcome, gps = gps, rp_y = rp_y, rp_a = rp_a)
  class(out) <- "csdr_learners"
  validate_csdr_learners(out)
  out
}

is_csdr_learner <- function(x) {
  inherits(x, "csdr_learner")
}

validate_csdr_learners <- function(learners) {
  if (!inherits(learners, "csdr_learners")) {
    stop("'learners' must be created by csdr_learners().", call. = FALSE)
  }
  require_regression_learner(learners$outcome, "outcome")
  require_gps_learner(learners$gps, "gps")
  require_regression_learner(learners$rp_y, "rp_y")
  require_regression_learner(learners$rp_a, "rp_a")
  invisible(TRUE)
}

#' @export
print.csdr_learners <- function(x, ...) {
  cat("CSDR learners\n")
  for (role in c("outcome", "gps", "rp_y", "rp_a")) {
    learner <- x[[role]]
    cat(sprintf("- %s: %s", role, learner$label))
    if (isTRUE(learner$is_default)) {
      cat(" (default)")
    }
    cat("\n")
  }
  invisible(x)
}

mark_default_learner <- function(x) {
  x$is_default <- TRUE
  x
}

require_regression_learner <- function(x, role) {
  if (!inherits(x, "csdr_regression_learner") || !is.function(x$fitter) ||
      !is.list(x$args)) {
    stop(sprintf("'%s' must be a CSDR regression learner.", role), call. = FALSE)
  }
  invisible(TRUE)
}

require_gps_learner <- function(x, role) {
  if (!inherits(x, "csdr_gps_learner") || !is.function(x$fitter) ||
      !is.list(x$args)) {
    stop(sprintf("'%s' must be a CSDR GPS learner.", role), call. = FALSE)
  }
  invisible(TRUE)
}
