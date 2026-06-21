ers_sample_data <- function(n = 50, p = 2, q = 1, seed = 2026) {
  set.seed(seed)
  C <- matrix(rnorm(n * q), n, q, dimnames = list(NULL, paste0("C", seq_len(q))))
  A <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("A", seq_len(p))))
  A[, 1] <- 0.5 * C[, 1] + A[, 1]

  Y <- 1 + 0.8 * A[, 1] - 0.3 * A[, 2] + 0.5 * C[, 1] + rnorm(n, sd = 0.1)

  list(Y = Y, A = A, C = C)
}

ers_fit_linear <- function(Y, XC_df, ...) {
  lm(Y ~ ., data = as.data.frame(XC_df))
}

test_that("estimate_ers L = 1 runs without cross-fitting", {
  d <- ers_sample_data(seed = 101)
  fit <- estimate_ers(
    Y = d$Y,
    A = d$A,
    C = d$C,
    a_eval = d$A[1:5, , drop = FALSE],
    estimator = "RA",
    L = 1,
    outcome_fitter = ers_fit_linear,
    verbose = FALSE
  )

  expect_s3_class(fit, "ers_fit")
  expect_equal(fit$L, 1L)
  expect_true(all(fit$folds == 1L))
  expect_equal(fit$diagnostics$nuisance_strategy, "fullfit")
  expect_equal(nrow(fit$results), 5L)
  expect_true(all(is.finite(fit$results$estimate)))
})

test_that("estimate_ers L > 1 runs with cross-fitting", {
  d <- ers_sample_data(seed = 102)
  fit <- estimate_ers(
    Y = d$Y,
    A = d$A,
    C = d$C,
    a_eval = d$A[1:6, , drop = FALSE],
    estimator = "RA",
    L = 2,
    outcome_fitter = ers_fit_linear,
    seed = 102,
    verbose = FALSE
  )

  expect_s3_class(fit, "ers_fit")
  expect_equal(fit$L, 2L)
  expect_equal(sort(unique(fit$folds)), 1:2)
  expect_equal(fit$diagnostics$nuisance_strategy, "crossfit")
  expect_equal(nrow(fit$results), 6L)
  expect_true(all(is.finite(fit$results$estimate)))
})

test_that("estimate_ers a_eval NULL evaluates observed exposure rows", {
  d <- ers_sample_data(n = 30, seed = 103)
  fit <- estimate_ers(
    Y = d$Y,
    A = d$A,
    C = d$C,
    estimator = "RA",
    L = 1,
    outcome_fitter = ers_fit_linear,
    verbose = FALSE
  )

  expect_equal(nrow(fit$results), nrow(d$A))
  expect_equal(names(fit$results)[1:2], colnames(d$A))
})

test_that("estimate_ers RA does not require GPS fitting", {
  d <- ers_sample_data(seed = 104)
  gps_stop <- function(...) stop("GPS should not be fit")

  fit <- estimate_ers(
    Y = d$Y,
    A = d$A,
    C = d$C,
    a_eval = d$A[1:4, , drop = FALSE],
    estimator = "RA",
    L = 2,
    outcome_fitter = ers_fit_linear,
    gps_fitter = gps_stop,
    seed = 104,
    verbose = FALSE
  )

  expect_s3_class(fit, "ers_fit")
  expect_true(all(is.finite(fit$results$estimate)))
  expect_true(all(is.na(fit$results$min_gps)))
})

test_that("estimate_ers IPW does not require outcome fitting", {
  d <- ers_sample_data(seed = 105)
  outcome_stop <- function(...) stop("Outcome should not be fit")

  fit <- estimate_ers(
    Y = d$Y,
    A = d$A,
    C = d$C,
    a_eval = d$A[1:4, , drop = FALSE],
    estimator = "IPW",
    L = 2,
    outcome_fitter = outcome_stop,
    gps_fitter = mvn_fitter,
    seed = 105,
    verbose = FALSE,
    args_gps = list(method_gps = "linear")
  )

  expect_s3_class(fit, "ers_fit")
  expect_equal(nrow(fit$results), 4L)
  expect_true(all(is.finite(fit$results$estimate)))
  expect_true(all(is.finite(fit$results$min_gps)))
})

test_that("estimate_ers DR uses both nuisances", {
  d <- ers_sample_data(seed = 106)
  fit <- estimate_ers(
    Y = d$Y,
    A = d$A,
    C = d$C,
    a_eval = d$A[1:4, , drop = FALSE],
    estimator = "DR",
    L = 2,
    outcome_fitter = ers_fit_linear,
    gps_fitter = mvn_fitter,
    seed = 106,
    verbose = FALSE,
    args_gps = list(method_gps = "linear")
  )

  expect_s3_class(fit, "ers_fit")
  expect_equal(nrow(fit$results), 4L)
  expect_true(all(is.finite(fit$results$estimate)))
  expect_true(all(is.finite(fit$results$se)))
})

test_that("estimate_ers returns classed object with results", {
  d <- ers_sample_data(seed = 107)
  fit <- estimate_ers(
    Y = d$Y,
    A = d$A,
    C = d$C,
    a_eval = d$A[1:3, , drop = FALSE],
    estimator = "RA",
    L = 1,
    outcome_fitter = ers_fit_linear,
    verbose = FALSE
  )

  expect_s3_class(fit, "ers_fit")
  expect_type(fit, "list")
  expect_s3_class(fit$results, "data.frame")
  expect_named(fit$results, c(
    "A1", "A2", "estimate", "se", "ci_lower", "ci_upper",
    "effective_n", "min_gps", "h_used_1", "h_used_2"
  ))
})

test_that("crossfit_ERS remains a compatibility wrapper", {
  d <- ers_sample_data(seed = 108)
  fit <- crossfit_ERS(
    Y = d$Y,
    X = d$A,
    C = d$C,
    x_eval = d$A[1:4, , drop = FALSE],
    estimator = "RA",
    L = 2,
    outcome_fitter = ers_fit_linear,
    seed = 108
  )

  expect_s3_class(fit, "ers_fit")
  expect_equal(nrow(fit$results), 4L)
  expect_equal(fit$metadata$L_folds, 2L)
  expect_equal(fit$metadata$seed, 108L)
  expect_true(all(is.finite(fit$results$estimate)))
})
