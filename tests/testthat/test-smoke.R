sample_data <- function(n = 60, p = 2, q = 1, seed = 2026) {
  set.seed(seed)
  C <- matrix(rnorm(n * q), n, q, dimnames = list(NULL, paste0("C", seq_len(q))))
  X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("X", seq_len(p))))
  X[, 1] <- 0.4 * C[, 1] + 0.5 * X[, 1]
  X <- as.matrix(X)

  noise <- rnorm(n, sd = 0.05)
  Y <- 1 + 0.6 * X[, 1] - 0.2 * X[, 2] + 0.4 * C[, 1] + noise

  list(Y = Y, X = X, C = C)
}

fit_linear <- function(Y, XC_df, ...) {
  lm(Y ~ ., data = as.data.frame(XC_df))
}

fit_c_linear <- function(target, C_df, ...) {
  lm(target ~ ., data = as.data.frame(C_df))
}

test_that("crossfit_ERS RA smoke test is deterministic and structured", {
  d <- sample_data()
  x_eval <- d$X[1:10, , drop = FALSE]
  out1 <- csdr:::crossfit_ERS(
    Y = d$Y,
    X = d$X,
    C = d$C,
    x_eval = x_eval,
    estimator = "RA",
    L = 3,
    outcome_fitter = fit_linear,
    seed = 2026
  )
  out2 <- csdr:::crossfit_ERS(
    Y = d$Y,
    X = d$X,
    C = d$C,
    x_eval = x_eval,
    estimator = "RA",
    L = 3,
    outcome_fitter = fit_linear,
    seed = 2026
  )

  expect_type(out1, "list")
  expect_type(out2, "list")
  expect_equal(out1$results, out2$results)
  expect_equal(out1$metadata$L_folds, 3L)
  expect_equal(out1$metadata$seed, 2026L)
  expect_equal(nrow(out1$results), nrow(x_eval))
  expect_true(all(is.finite(out1$results$estimate)))
  expect_true(all(!is.na(out1$results$estimate)))
})

test_that("outcome_model and predict output match expectation", {
  d <- sample_data(seed = 2027)
  out_model <- csdr:::outcome_model(
    Y = d$Y,
    X = d$X,
    C = d$C,
    mu_fitter = fit_linear
  )

  pred <- predict(out_model, newdata = cbind(as.data.frame(d$X), as.data.frame(d$C)))

  expect_s3_class(out_model, "outcome_model")
  expect_equal(length(pred), nrow(d$X))
  expect_true(all(is.finite(pred)))
  expect_false(anyNA(pred))
})

test_that("csdr main public API returns a csdr_fit for one variant", {
  d <- sample_data(seed = 3030)

  fit <- suppressWarnings(
    csdr(
      Y = d$Y,
      A = d$X,
      C = d$C,
      variants = "DR",
      d = 1,
      L = 2L,
      seed = 3030,
      learners = csdr_learners(
        outcome = custom_regression(fit_linear, label = "linear outcome")
      ),
      verbose = FALSE
    )
  )

  expect_s3_class(fit, "csdr_fit")
  expect_s3_class(fit$target, "csdr_target")
  expect_s3_class(fit$learners, "csdr_learners")
  expect_equal(fit$variants, "DR")
  expect_equal(names(fit$fits), "DR")
  expect_true(is.data.frame(fit$learner_summary))
  expect_true("outcome" %in% fit$learner_summary$role)
  expect_equal(nrow(fit$summary), 1L)
  expect_equal(fit$summary$variant, "DR")
  expect_equal(fit$summary$d_hat, 1L)
  expect_equal(fit$fits$DR$d_hat, 1L)
  expect_true(is.matrix(coef(fit, variant = "DR")))
  expect_equal(ncol(coef(fit, variant = "DR")), 1L)
  expect_equal(nrow(scores(fit, variant = "DR")), nrow(d$X))

  target <- targets(fit, variant = "DR")
  expect_equal(length(target$target_Y), nrow(d$X))
  expect_equal(dim(target$target_A), dim(d$X))
})

test_that("csdr supports multiple variants and storage controls", {
  d <- sample_data(n = 40, seed = 3031)

  fit <- suppressWarnings(
    csdr(
      Y = d$Y,
      A = d$X,
      C = d$C,
      variants = c("RA", "DR", "PO", "RP"),
      d = 1,
      L = 2L,
      seed = 3031,
      learners = csdr_learners(
        outcome = custom_regression(fit_linear, label = "linear outcome"),
        rp_y = custom_regression(fit_c_linear, label = "linear C regression")
      ),
      keep_mave = FALSE,
      verbose = FALSE
    )
  )

  expect_s3_class(fit, "csdr_fit")
  expect_equal(names(fit$fits), c("RA", "DR", "PO", "RP"))
  expect_equal(nrow(fit$summary), 4L)
  expect_equal(fit$summary$variant, c("RA", "DR", "PO", "RP"))
  expect_true(all(fit$summary$d_hat == 1L))
  expect_true(all(vapply(coef(fit), is.matrix, logical(1))))
  expect_equal(dim(scores(fit, variant = "PO")), c(nrow(d$X), 1L))
  expect_null(fit$fits$DR$mave_fit)
  expect_null(fit$fits$DR$mave_dim_obj)
  expect_true(is.matrix(fit$fits$DR$beta))
  expect_true(is.matrix(fit$fits$DR$score))
  expect_equal(fit$summary$target_exposure[fit$summary$variant == "RP"], "residualized_A")
  expect_equal(fit$summary$target_exposure[fit$summary$variant == "DR"], "A")
  expect_true(fit$learner_summary$used[fit$learner_summary$role == "rp_a"])
  expect_match(
    fit$learner_summary$details[fit$learner_summary$role == "rp_a"],
    "not separately wired"
  )
})

test_that("csdr targets accessor errors when targets are not retained", {
  d <- sample_data(n = 40, seed = 3032)

  fit <- suppressWarnings(
    csdr(
      Y = d$Y,
      A = d$X,
      C = d$C,
      variants = "RA",
      d = 1,
      L = 2L,
      seed = 3032,
      learners = csdr_learners(
        outcome = custom_regression(fit_linear, label = "linear outcome")
      ),
      keep_targets = FALSE,
      verbose = FALSE
    )
  )

  expect_null(fit$fits$RA$target_Y)
  expect_error(targets(fit, variant = "RA"), "Targets were not retained")
})

test_that("csdr works with default and modified learner specs", {
  d <- sample_data(n = 36, seed = 3033)

  fit_default <- suppressWarnings(
    csdr(
      Y = d$Y,
      A = d$X,
      C = d$C,
      variants = "DR",
      d = 1,
      L = 2L,
      seed = 3033,
      verbose = FALSE
    )
  )
  fit_glm <- suppressWarnings(
    csdr(
      Y = d$Y,
      A = d$X,
      C = d$C,
      variants = "DR",
      d = 1,
      L = 2L,
      seed = 3033,
      learners = csdr_learners(sl_library = c("SL.glm")),
      verbose = FALSE
    )
  )

  expect_s3_class(fit_default, "csdr_fit")
  expect_s3_class(fit_glm, "csdr_fit")
  expect_equal(
    fit_default$learner_summary$engine[fit_default$learner_summary$role == "outcome"],
    "SuperLearner"
  )
  expect_equal(
    fit_default$learner_summary$engine[fit_default$learner_summary$role == "gps"],
    "MVN GPS"
  )
  expect_equal(fit_glm$learners$outcome$args$SL.library, "SL.glm")
})

test_that("csdr learner-choice messages respect verbose", {
  d <- sample_data(n = 32, seed = 3034)

  quiet_messages <- character()
  withCallingHandlers(
    suppressWarnings(
      csdr(
        Y = d$Y,
        A = d$X,
        C = d$C,
        variants = "DR",
        d = 1,
        L = 2L,
        seed = 3034,
        learners = csdr_learners(sl_library = c("SL.glm")),
        verbose = FALSE
      )
    ),
    message = function(m) {
      quiet_messages <<- c(quiet_messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  expect_false(any(grepl("CSDR learner choices", quiet_messages, fixed = TRUE)))
  expect_message(
    suppressWarnings(
      csdr(
        Y = d$Y,
        A = d$X,
        C = d$C,
        variants = "DR",
        d = 1,
        L = 2L,
        seed = 3034,
        learners = csdr_learners(sl_library = c("SL.glm")),
        verbose = TRUE
      )
    ),
    "CSDR learner choices"
  )
})

test_that("estimate_ERS RA returns finite numeric vector and is reproducible", {
  d <- sample_data(seed = 4040)
  x_eval <- d$X[1:4, , drop = FALSE]

  out_model <- outcome_model(
    Y = d$Y,
    X = d$X,
    C = d$C,
    mu_fitter = fit_linear
  )

  est1 <- suppressWarnings(
    estimate_ERS(
    Y = d$Y,
    X = d$X,
    C = d$C,
    estimator = "RA",
    x_eval = x_eval,
    return_vector = TRUE,
    out_model = out_model
    )
  )

  est2 <- suppressWarnings(
    estimate_ERS(
    Y = d$Y,
    X = d$X,
    C = d$C,
    estimator = "RA",
    x_eval = x_eval,
    return_vector = TRUE,
    out_model = out_model
    )
  )

  expect_type(est1, "double")
  expect_equal(length(est1), nrow(x_eval))
  expect_true(all(is.finite(est1)))
  expect_equal(est1, est2)
})
