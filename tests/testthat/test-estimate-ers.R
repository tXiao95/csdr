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

ers_linear_truth_data <- function(n = 400, seed = 2027, noise_sd = 0) {
  set.seed(seed)
  C <- matrix(rnorm(n), n, 1, dimnames = list(NULL, "C1"))
  A <- matrix(rnorm(n * 2), n, 2, dimnames = list(NULL, c("A1", "A2")))
  A[, 1] <- 0.7 * C[, 1] + A[, 1]
  A[, 2] <- -0.4 * C[, 1] + A[, 2]
  Y <- 1 + 2 * A[, 1] - A[, 2] + 3 * C[, 1] + rnorm(n, sd = noise_sd)

  list(Y = Y, A = A, C = C)
}

ers_linear_truth <- function(a_eval, C) {
  a_eval <- as.matrix(a_eval)
  as.numeric(1 + 2 * a_eval[, 1] - a_eval[, 2] + 3 * mean(C[, 1]))
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

test_that("estimate_ers RA matches known linear g-formula", {
  d <- ers_linear_truth_data(seed = 201, noise_sd = 0.03)
  a_eval <- matrix(
    c(-1, 0.5,
       0, 0,
       1.5, -0.5),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("A1", "A2"))
  )
  truth <- ers_linear_truth(a_eval, d$C)

  fit <- estimate_ers(
    Y = d$Y,
    A = d$A,
    C = d$C,
    a_eval = a_eval,
    estimator = "RA",
    L = 1,
    outcome_fitter = ers_fit_linear,
    verbose = FALSE
  )

  expect_equal(fit$results$estimate, truth, tolerance = 0.05)
})

test_that("estimate_ers DR agrees with RA when outcome model is correct", {
  d <- ers_linear_truth_data(seed = 202, noise_sd = 0)
  a_eval <- matrix(
    c(-0.5, 0.25,
       0.75, -1),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("A1", "A2"))
  )
  truth <- ers_linear_truth(a_eval, d$C)

  fit_ra <- estimate_ers(
    Y = d$Y,
    A = d$A,
    C = d$C,
    a_eval = a_eval,
    estimator = "RA",
    L = 1,
    outcome_fitter = ers_fit_linear,
    verbose = FALSE
  )
  fit_dr <- estimate_ers(
    Y = d$Y,
    A = d$A,
    C = d$C,
    a_eval = a_eval,
    estimator = "DR",
    L = 1,
    outcome_fitter = ers_fit_linear,
    gps_fitter = mvn_fitter,
    h = c(0.15, 0.15),
    verbose = FALSE,
    args_gps = list(method_gps = "linear")
  )

  expect_equal(fit_dr$results$estimate, fit_ra$results$estimate, tolerance = 0.02)
  expect_equal(fit_ra$results$estimate, truth, tolerance = 1e-8)
  expect_equal(fit_dr$results$estimate, truth, tolerance = 0.02)
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
  expect_equal(as.matrix(fit$results[, 1:2]), d$A, tolerance = 0)
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

test_that("estimate_ers supports one-dimensional exposure input forms", {
  set.seed(301)
  n <- 80
  C <- data.frame(C1 = rnorm(n))
  A_vec <- 0.4 * C$C1 + rnorm(n)
  Y <- 1 + 1.5 * A_vec + 2 * C$C1

  a_vec <- c(-1, 0, 1)
  a_mat <- matrix(a_vec, ncol = 1, dimnames = list(NULL, "A1"))
  a_df <- data.frame(A1 = a_vec)

  fit_vec <- estimate_ers(
    Y = Y,
    A = A_vec,
    C = C,
    a_eval = a_vec,
    estimator = "RA",
    L = 1,
    outcome_fitter = ers_fit_linear,
    verbose = FALSE
  )
  fit_mat <- estimate_ers(
    Y = Y,
    A = matrix(A_vec, ncol = 1, dimnames = list(NULL, "A1")),
    C = C,
    a_eval = a_mat,
    estimator = "RA",
    L = 1,
    outcome_fitter = ers_fit_linear,
    verbose = FALSE
  )
  fit_df <- estimate_ers(
    Y = Y,
    A = data.frame(A1 = A_vec),
    C = C,
    a_eval = a_df,
    estimator = "RA",
    L = 1,
    outcome_fitter = ers_fit_linear,
    verbose = FALSE
  )

  expect_equal(nrow(fit_vec$results), length(a_vec))
  expect_equal(nrow(fit_mat$results), length(a_vec))
  expect_equal(nrow(fit_df$results), length(a_vec))
  expect_equal(fit_vec$results$estimate, fit_mat$results$estimate, tolerance = 1e-8)
  expect_equal(fit_mat$results$estimate, fit_df$results$estimate, tolerance = 1e-8)
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

test_that("estimate_ers validates inputs explicitly", {
  d <- ers_sample_data(n = 20, seed = 401)

  expect_error(
    estimate_ers(d$Y[-1], d$A, d$C, estimator = "RA", L = 1,
                 outcome_fitter = ers_fit_linear, verbose = FALSE),
    "'A' must have the same number of rows as 'Y'",
    fixed = TRUE
  )
  expect_error(
    estimate_ers(d$Y, d$A, d$C[-1, , drop = FALSE], estimator = "RA", L = 1,
                 outcome_fitter = ers_fit_linear, verbose = FALSE),
    "'C' must have the same number of rows as 'Y'",
    fixed = TRUE
  )
  expect_error(
    estimate_ers(d$Y, d$A, d$C, a_eval = matrix(0, nrow = 2, ncol = 1),
                 estimator = "RA", L = 1, outcome_fitter = ers_fit_linear,
                 verbose = FALSE),
    "'a_eval' must have the same number of columns as 'A'",
    fixed = TRUE
  )
  expect_error(
    estimate_ers(d$Y, d$A, d$C, estimator = "RA", L = 1,
                 outcome_fitter = ers_fit_linear, h = c(1, 1, 1),
                 verbose = FALSE),
    "'h' must have length 1 or ncol(A)",
    fixed = TRUE
  )
  expect_error(
    estimate_ers(d$Y, d$A, d$C, estimator = "RA", L = 1,
                 outcome_fitter = ers_fit_linear, h = c(1, 0),
                 verbose = FALSE),
    "'h' must contain positive finite bandwidths",
    fixed = TRUE
  )
  expect_error(
    estimate_ers(d$Y, d$A, d$C, estimator = "RA", L = 0,
                 outcome_fitter = ers_fit_linear, verbose = FALSE),
    "'L' must be a positive integer",
    fixed = TRUE
  )
  expect_error(
    estimate_ers(d$Y, d$A, d$C, estimator = "RA", L = nrow(d$A) + 1,
                 outcome_fitter = ers_fit_linear, verbose = FALSE),
    "'L' cannot exceed the number of observations",
    fixed = TRUE
  )
  expect_error(
    estimate_ers(d$Y, d$A, d$C, estimator = "RA", L = 2,
                 folds = rep(1:2, length.out = nrow(d$A) - 1),
                 outcome_fitter = ers_fit_linear, verbose = FALSE),
    "'folds' must have length n",
    fixed = TRUE
  )
  expect_error(
    estimate_ers(d$Y, d$A, d$C, estimator = "RA", L = 2,
                 folds = rep(1, nrow(d$A)),
                 outcome_fitter = ers_fit_linear, verbose = FALSE),
    "'folds' must contain at least two folds when L > 1",
    fixed = TRUE
  )
  expect_error(
    estimate_ers(d$Y, d$A, d$C, estimator = "IPW", L = 1,
                 gps_fitter = mvn_fitter, gps_floor = 0, verbose = FALSE,
                 args_gps = list(method_gps = "linear")),
    "'gps_floor' must be a positive finite number",
    fixed = TRUE
  )

  A_zero_var <- d$A
  A_zero_var[, 2] <- 1
  expect_error(
    estimate_ers(d$Y, A_zero_var, d$C, estimator = "RA", L = 1,
                 outcome_fitter = ers_fit_linear, verbose = FALSE),
    "Default bandwidth is nonpositive or nonfinite",
    fixed = TRUE
  )
})

test_that("estimate_ers diagnostics are populated by estimator", {
  d <- ers_sample_data(seed = 501)

  fit_ra <- estimate_ers(
    Y = d$Y,
    A = d$A,
    C = d$C,
    a_eval = d$A[1:4, , drop = FALSE],
    estimator = "RA",
    L = 1,
    outcome_fitter = ers_fit_linear,
    verbose = FALSE
  )
  fit_dr <- estimate_ers(
    Y = d$Y,
    A = d$A,
    C = d$C,
    a_eval = d$A[1:4, , drop = FALSE],
    estimator = "DR",
    L = 1,
    outcome_fitter = ers_fit_linear,
    gps_fitter = mvn_fitter,
    verbose = FALSE,
    args_gps = list(method_gps = "linear")
  )

  required_cols <- c(
    "estimate", "se", "ci_lower", "ci_upper",
    "effective_n", "min_gps", "h_used_1", "h_used_2"
  )
  expect_true(all(required_cols %in% names(fit_ra$results)))
  expect_true(all(required_cols %in% names(fit_dr$results)))
  expect_true(all(is.na(fit_ra$results$min_gps)))
  expect_true(all(is.finite(fit_dr$results$min_gps)))
  expect_true(all(is.finite(fit_ra$results$effective_n)))
  expect_true(all(is.finite(fit_dr$results$effective_n)))
  expect_true(all(fit_ra$results$effective_n >= 0))
  expect_true(all(fit_dr$results$effective_n >= 0))
})
