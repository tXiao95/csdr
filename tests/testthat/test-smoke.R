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
