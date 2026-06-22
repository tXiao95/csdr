target_sample_data <- function(n = 50, p = 2, seed = 6101) {
  set.seed(seed)
  C <- matrix(rnorm(n), n, 1, dimnames = list(NULL, "C1"))
  A <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("A", seq_len(p))))
  A[, 1] <- 0.5 * C[, 1] + A[, 1]
  if (p > 1L) {
    A[, 2] <- -0.25 * C[, 1] + A[, 2]
  }
  Y <- 1 + 0.75 * A[, 1] + 0.4 * C[, 1] + rnorm(n, sd = 0.05)
  list(Y = Y, A = A, C = C)
}

target_fit_linear <- function(Y, XC_df, ...) {
  lm(Y ~ ., data = as.data.frame(XC_df))
}

target_fit_c_linear <- function(target, C_df, ...) {
  lm(target ~ ., data = as.data.frame(C_df))
}

test_that("compute_pseudo_outcomes matches manual arithmetic", {
  Y <- c(2, 3, 5)
  m_obs <- c(1, 2, 4)
  pi_obs <- c(0.5, 0.25, 0)
  m_marginal <- c(1.5, 2.5, 4.5)
  pi_marginal <- c(0.4, 0.2, 0.1)
  gps_floor <- 0.01

  out <- compute_pseudo_outcomes(
    Y = Y,
    m_obs = m_obs,
    pi_obs = pi_obs,
    m_marginal = m_marginal,
    pi_marginal = pi_marginal,
    gps_floor = gps_floor
  )
  manual <- pi_marginal * (Y - m_obs) / pmax(pi_obs, gps_floor) + m_marginal

  expect_equal(out, manual)
})

test_that("predict_po_nuisance returns numeric vectors of length n", {
  d <- target_sample_data(seed = 6102)
  nuis <- fit_po_nuisance(
    Y = d$Y,
    A = d$A,
    C = d$C,
    outcome_fitter = target_fit_linear,
    gps_fitter = mvn_fitter,
    args_gps = list(method_gps = "linear"),
    verbose = FALSE
  )

  pred <- predict_po_nuisance(
    nuisance = nuis,
    A = d$A,
    C = d$C,
    gps_floor = 1e-8,
    verbose = FALSE
  )

  expect_equal(names(pred), c("m_obs", "pi_obs", "m_marginal", "pi_marginal"))
  expect_equal(length(pred$m_obs), nrow(d$A))
  expect_equal(length(pred$pi_obs), nrow(d$A))
  expect_equal(length(pred$m_marginal), nrow(d$A))
  expect_equal(length(pred$pi_marginal), nrow(d$A))
  expect_true(all(is.finite(pred$m_obs)))
  expect_true(all(is.finite(pred$pi_obs)))
})

test_that("estimate_pseudo_outcomes supports full-fit and cross-fit modes", {
  d <- target_sample_data(seed = 6103)

  fit_l1 <- estimate_pseudo_outcomes(
    Y = d$Y,
    A = d$A,
    C = d$C,
    L = 1,
    outcome_fitter = target_fit_linear,
    gps_fitter = mvn_fitter,
    args_gps = list(method_gps = "linear"),
    return_nuisance = TRUE,
    verbose = FALSE
  )
  fit_l2 <- estimate_pseudo_outcomes(
    Y = d$Y,
    A = d$A,
    C = d$C,
    L = 2,
    outcome_fitter = target_fit_linear,
    gps_fitter = mvn_fitter,
    args_gps = list(method_gps = "linear"),
    seed = 6103,
    verbose = FALSE
  )

  expect_s3_class(fit_l1, "po_fit")
  expect_s3_class(fit_l1, "csdr_target")
  expect_equal(length(fit_l1$pseudo_outcomes), nrow(d$A))
  expect_equal(length(fit_l2$pseudo_outcomes), nrow(d$A))
  expect_equal(length(unique(fit_l2$folds)), 2L)
  expect_true(all(is.finite(fit_l1$pseudo_outcomes)))
  expect_type(fit_l1$nuisance, "list")
})

test_that("estimate_pseudo_outcomes supports pseudo-outcome marginalization modes", {
  d <- target_sample_data(n = 36, seed = 61031)
  folds <- rep(1:2, length.out = nrow(d$A))

  fit_default <- estimate_pseudo_outcomes(
    Y = d$Y,
    A = d$A,
    C = d$C,
    L = 2,
    folds = folds,
    outcome_fitter = target_fit_linear,
    gps_fitter = mvn_fitter,
    args_gps = list(method_gps = "linear"),
    verbose = FALSE
  )
  fit_crossfit <- estimate_pseudo_outcomes(
    Y = d$Y,
    A = d$A,
    C = d$C,
    L = 2,
    folds = folds,
    outcome_fitter = target_fit_linear,
    gps_fitter = mvn_fitter,
    marginalization = "crossfit",
    args_gps = list(method_gps = "linear"),
    verbose = FALSE
  )
  fit_fold <- estimate_pseudo_outcomes(
    Y = d$Y,
    A = d$A,
    C = d$C,
    L = 2,
    folds = folds,
    outcome_fitter = target_fit_linear,
    gps_fitter = mvn_fitter,
    marginalization = "fold",
    args_gps = list(method_gps = "linear"),
    verbose = FALSE
  )
  fit_all <- estimate_pseudo_outcomes(
    Y = d$Y,
    A = d$A,
    C = d$C,
    L = 2,
    folds = folds,
    outcome_fitter = target_fit_linear,
    gps_fitter = mvn_fitter,
    marginalization = "all",
    args_gps = list(method_gps = "linear"),
    verbose = FALSE
  )

  expect_equal(length(fit_crossfit$pseudo_outcomes), nrow(d$A))
  expect_equal(length(fit_fold$pseudo_outcomes), nrow(d$A))
  expect_equal(length(fit_all$pseudo_outcomes), nrow(d$A))
  expect_equal(fit_default$pseudo_outcomes, fit_crossfit$pseudo_outcomes)
  expect_equal(fit_crossfit$diagnostics$marginalization, "crossfit")
  expect_equal(fit_fold$diagnostics$marginalization, "fold")
  expect_equal(fit_all$diagnostics$marginalization, "all")
})

test_that("estimate_pseudo_outcomes validates inputs", {
  d <- target_sample_data(n = 20, seed = 6104)

  expect_error(
    compute_pseudo_outcomes(d$Y, d$Y[-1], d$Y, d$Y, d$Y),
    "must have the same length",
    fixed = TRUE
  )
  expect_error(
    estimate_pseudo_outcomes(
      Y = d$Y[-1],
      A = d$A,
      C = d$C,
      L = 1,
      outcome_fitter = target_fit_linear,
      gps_fitter = mvn_fitter,
      verbose = FALSE
    ),
    "'A' must have the same number of rows as 'Y'",
    fixed = TRUE
  )
  expect_error(
    estimate_pseudo_outcomes(
      Y = d$Y,
      A = d$A,
      C = d$C,
      L = 0,
      outcome_fitter = target_fit_linear,
      gps_fitter = mvn_fitter,
      verbose = FALSE
    ),
    "'L' must be a positive integer",
    fixed = TRUE
  )
  expect_error(
    estimate_pseudo_outcomes(
      Y = d$Y,
      A = d$A,
      C = d$C,
      L = 1,
      outcome_fitter = target_fit_linear,
      gps_fitter = mvn_fitter,
      gps_floor = 0,
      verbose = FALSE
    ),
    "'gps_floor' must be a positive finite number",
    fixed = TRUE
  )
  expect_error(
    estimate_pseudo_outcomes(
      Y = d$Y,
      A = d$A,
      C = d$C,
      L = 1,
      outcome_fitter = target_fit_linear,
      gps_fitter = mvn_fitter,
      marginalization = "bad",
      verbose = FALSE
    ),
    "'arg' should be one of",
    fixed = TRUE
  )
})

test_that("legacy pseudo-outcome call path still returns a vector", {
  d <- target_sample_data(seed = 6105)
  out_mod <- outcome_model(Y = d$Y, X = d$A, C = d$C, mu_fitter = target_fit_linear)
  gps_mod <- gps_model(
    X = d$A,
    C = d$C,
    pi_fitter = mvn_fitter,
    method_gps = "linear"
  )

  xi <- estimate_pseudo_outcomes(
    Y = d$Y,
    X = d$A,
    C = d$C,
    out_model = out_mod,
    gps_model = gps_mod,
    delta_n = 1e-8,
    verbose = FALSE
  )

  expect_type(xi, "double")
  expect_equal(length(xi), nrow(d$A))
})

test_that("compute_residualized_pair matches manual subtraction", {
  Y <- c(1, 3, 6)
  A <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3)
  EY_C <- c(0.5, 1, 2)
  EA_C <- matrix(1, nrow = 3, ncol = 2)

  out <- compute_residualized_pair(Y = Y, A = A, EY_C = EY_C, EA_C = EA_C)

  expect_equal(out$Y_tilde, Y - EY_C)
  expect_equal(unname(out$A_tilde), unname(A - EA_C))
})

test_that("fit and predict residualized-pair nuisances return expected shapes", {
  d <- target_sample_data(seed = 6201)
  nuis <- fit_rp_nuisance(
    Y = d$Y,
    A = d$A,
    C = d$C,
    y_fitter = target_fit_c_linear,
    a_fitter = target_fit_c_linear,
    verbose = FALSE
  )
  pred <- predict_rp_nuisance(nuisance = nuis, C = d$C)

  expect_s3_class(nuis, "rp_nuisance")
  expect_equal(length(nuis$A_models), ncol(d$A))
  expect_equal(length(pred$EY_C), nrow(d$A))
  expect_equal(dim(pred$EA_C), dim(d$A))
  expect_true(all(is.finite(pred$EY_C)))
  expect_true(all(is.finite(pred$EA_C)))
})

test_that("estimate_residualized_pair supports full-fit and cross-fit modes", {
  d <- target_sample_data(seed = 6202)

  fit_l1 <- estimate_residualized_pair(
    Y = d$Y,
    A = d$A,
    C = d$C,
    L = 1,
    C_fitter = target_fit_c_linear,
    return_nuisance = TRUE,
    verbose = FALSE
  )
  fit_l2 <- estimate_residualized_pair(
    Y = d$Y,
    A = d$A,
    C = d$C,
    L = 2,
    C_fitter = target_fit_c_linear,
    seed = 6202,
    verbose = FALSE
  )

  expect_s3_class(fit_l1, "rp_fit")
  expect_s3_class(fit_l1, "csdr_target")
  expect_equal(length(fit_l1$Y_tilde), nrow(d$A))
  expect_equal(dim(fit_l1$A_tilde), dim(d$A))
  expect_equal(length(fit_l2$Y_tilde), nrow(d$A))
  expect_equal(dim(fit_l2$A_tilde), dim(d$A))
  expect_equal(length(unique(fit_l2$folds)), 2L)
  expect_type(fit_l1$nuisance, "list")
})

test_that("estimate_residualized_pair validates inputs", {
  d <- target_sample_data(n = 20, seed = 6203)

  expect_error(
    compute_residualized_pair(
      Y = d$Y,
      A = d$A,
      EY_C = d$Y[-1],
      EA_C = d$A
    ),
    "'Y' and 'EY_C' must have the same length",
    fixed = TRUE
  )
  expect_error(
    compute_residualized_pair(
      Y = d$Y,
      A = d$A,
      EY_C = d$Y,
      EA_C = d$A[, 1, drop = FALSE]
    ),
    "'EA_C' must have the same dimensions as 'A'",
    fixed = TRUE
  )
  expect_error(
    estimate_residualized_pair(
      Y = d$Y,
      A = d$A[-1, , drop = FALSE],
      C = d$C,
      C_fitter = target_fit_c_linear,
      verbose = FALSE
    ),
    "'A' must have the same number of rows as 'Y'",
    fixed = TRUE
  )
  expect_error(
    estimate_residualized_pair(
      Y = d$Y,
      A = d$A,
      C = d$C,
      L = nrow(d$A) + 1,
      C_fitter = target_fit_c_linear,
      verbose = FALSE
    ),
    "'L' cannot exceed the number of observations",
    fixed = TRUE
  )
})

test_that("legacy residualized-pair call path still returns aliases", {
  d <- target_sample_data(seed = 6204)
  nuis <- train_nuisance_models(
    Y = d$Y,
    X = d$A,
    C = d$C,
    fitter = target_fit_c_linear
  )

  out <- estimate_residualized_pair(
    Y = d$Y,
    X = d$A,
    C = d$C,
    C_models = nuis,
    verbose = FALSE
  )

  expect_equal(length(out$Ytilde), nrow(d$A))
  expect_equal(dim(out$Xtilde), dim(d$A))
})

test_that("csdr_target constructs individual method targets", {
  d <- target_sample_data(n = 36, seed = 6251)

  fit_ra <- csdr_target(
    Y = d$Y,
    A = d$A,
    C = d$C,
    methods = "RA",
    L = 2,
    outcome_fitter = target_fit_linear,
    folds = rep(1:2, length.out = nrow(d$A)),
    verbose = FALSE
  )
  fit_dr <- csdr_target(
    Y = d$Y,
    A = d$A,
    C = d$C,
    methods = "DR",
    L = 2,
    outcome_fitter = target_fit_linear,
    gps_fitter = mvn_fitter,
    folds = rep(1:2, length.out = nrow(d$A)),
    args_gps = list(method_gps = "linear"),
    verbose = FALSE
  )
  fit_po <- csdr_target(
    Y = d$Y,
    A = d$A,
    C = d$C,
    methods = "PO",
    L = 2,
    outcome_fitter = target_fit_linear,
    gps_fitter = mvn_fitter,
    folds = rep(1:2, length.out = nrow(d$A)),
    args_gps = list(method_gps = "linear"),
    verbose = FALSE
  )
  fit_rp <- csdr_target(
    Y = d$Y,
    A = d$A,
    C = d$C,
    methods = "RP",
    L = 2,
    C_fitter = target_fit_c_linear,
    folds = rep(1:2, length.out = nrow(d$A)),
    verbose = FALSE
  )

  expect_s3_class(fit_ra, "csdr_target")
  expect_equal(names(fit_ra$new_Y), "RA")
  expect_equal(length(fit_ra$new_Y$RA), nrow(d$A))
  expect_equal(dim(fit_ra$new_A$RA), dim(d$A))
  expect_equal(length(fit_dr$new_Y$DR), nrow(d$A))
  expect_equal(dim(fit_dr$new_A$DR), dim(d$A))
  expect_equal(length(fit_po$new_Y$PO), nrow(d$A))
  expect_equal(dim(fit_po$new_A$PO), dim(d$A))
  expect_equal(fit_po$diagnostics$po_marginalization, "crossfit")
  expect_equal(length(fit_rp$new_Y$RP), nrow(d$A))
  expect_equal(dim(fit_rp$new_A$RP), dim(d$A))
  expect_false(isTRUE(all.equal(fit_rp$new_A$RP, d$A)))
})

test_that("csdr_target constructs all requested targets with shared folds", {
  d <- target_sample_data(n = 36, seed = 6252)
  folds <- rep(1:2, length.out = nrow(d$A))

  fit <- csdr_target(
    Y = d$Y,
    A = d$A,
    C = d$C,
    methods = c("RA", "DR", "PO", "RP"),
    L = 2,
    folds = folds,
    outcome_fitter = target_fit_linear,
    gps_fitter = mvn_fitter,
    C_fitter = target_fit_c_linear,
    args_gps = list(method_gps = "linear"),
    return_nuisance = TRUE,
    verbose = FALSE
  )

  expect_equal(names(fit$new_Y), c("RA", "DR", "PO", "RP"))
  expect_equal(names(fit$new_A), c("RA", "DR", "PO", "RP"))
  expect_equal(fit$folds, folds)
  expect_equal(fit$methods, c("RA", "DR", "PO", "RP"))
  expect_equal(fit$diagnostics$L, 2L)
  expect_equal(dim(fit$new_A$RA), dim(d$A))
  expect_equal(dim(fit$new_A$DR), dim(d$A))
  expect_equal(dim(fit$new_A$PO), dim(d$A))
  expect_equal(dim(fit$new_A$RP), dim(d$A))
  expect_false(is.null(fit$nuisance))
})

test_that("compute_new_response_and_exposure works for PO, RP, and combined targets", {
  d <- target_sample_data(n = 36, seed = 6301)

  po <- compute_new_response_and_exposure(
    Y = d$Y,
    X = d$A,
    C = d$C,
    method = "PO",
    L = 2,
    outcome_fitter = target_fit_linear,
    gps_fitter = mvn_fitter,
    seed = 6301,
    args_gps = list(method_gps = "linear")
  )
  rp <- compute_new_response_and_exposure(
    Y = d$Y,
    X = d$A,
    C = d$C,
    method = "RP",
    L = 2,
    C_fitter = target_fit_c_linear,
    seed = 6301
  )
  all_targets <- compute_new_response_and_exposure(
    Y = d$Y,
    X = d$A,
    C = d$C,
    method = c("RA", "DR", "PO", "RP"),
    L = 2,
    outcome_fitter = target_fit_linear,
    gps_fitter = mvn_fitter,
    C_fitter = target_fit_c_linear,
    seed = 6301,
    args_gps = list(method_gps = "linear")
  )

  expect_equal(length(po$new_Y), nrow(d$A))
  expect_equal(dim(po$new_X), dim(d$A))
  expect_equal(length(rp$new_Y), nrow(d$A))
  expect_equal(dim(rp$new_X), dim(d$A))
  expect_equal(names(all_targets$new_Y), c("RA", "DR", "PO", "RP"))
  expect_equal(names(all_targets$new_X), c("RA", "DR", "PO", "RP"))
  expect_equal(dim(all_targets$new_X$RA), dim(d$A))
  expect_equal(dim(all_targets$new_X$DR), dim(d$A))
  expect_equal(dim(all_targets$new_X$PO), dim(d$A))
  expect_equal(dim(all_targets$new_X$RP), dim(d$A))
  expect_equal(all_targets$new_X$RA, d$A)
  expect_equal(all_targets$new_X$PO, d$A)
})
