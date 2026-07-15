learner_toy_data <- function(n = 30, seed = 7201) {
  set.seed(seed)
  C <- data.frame(C1 = rnorm(n), C2 = rnorm(n))
  A <- data.frame(A1 = 0.3 * C$C1 + rnorm(n), A2 = -0.2 * C$C2 + rnorm(n))
  Y <- 1 + 0.5 * A$A1 + 0.2 * C$C1 + rnorm(n, sd = 0.1)
  list(Y = Y, A = A, C = C)
}

test_that("learner spec constructors validate and store expected fields", {
  expect_equal(
    csdr_default_sl_library(),
    c("SL.glm", "SL.glmnet", "SL.xgboost", "SL.earth")
  )

  reg <- sl_regression(SL.library = c("SL.glm"))
  gps <- mvn_gps(method = "linear")

  expect_s3_class(reg, "csdr_regression_learner")
  expect_s3_class(reg, "csdr_learner")
  expect_equal(reg$engine, "SuperLearner")
  expect_equal(reg$args$SL.library, "SL.glm")
  expect_s3_class(gps, "csdr_gps_learner")
  expect_equal(gps$engine, "MVN GPS")
  expect_equal(gps$args$method_gps, "linear")
  expect_equal(gps$args$delta_n, 1e-16)

  expect_error(custom_regression("not a function"), "'fitter' must be a function")
  expect_error(custom_gps("not a function"), "'fitter' must be a function")
  expect_error(custom_regression(function(...) NULL, args = 1), "'args' must be a list")
  expect_error(custom_gps(function(...) NULL, args = 1), "'args' must be a list")

  custom_gps_fit <- custom_gps(function(A, C, ...) {
    structure(list(n = nrow(A), q = ncol(C)), class = "toy_gps")
  })
  toy <- custom_gps_fit$fitter(X = data.frame(A1 = 1:3), C = data.frame(C1 = 4:6))
  expect_s3_class(toy, "toy_gps")
  expect_equal(toy$n, 3L)
})

test_that("csdr_learners returns and validates default learner roles", {
  learners <- csdr_learners()

  expect_s3_class(learners, "csdr_learners")
  expect_true(all(vapply(learners, is_csdr_learner, logical(1))))
  expect_true(learners$outcome$is_default)
  expect_true(learners$gps$is_default)
  expect_true(learners$rp_y$is_default)
  expect_true(learners$rp_a$is_default)
  expect_silent(validate_csdr_learners(learners))

  glm_only <- csdr_learners(sl_library = c("SL.glm"))
  expect_equal(glm_only$outcome$args$SL.library, "SL.glm")
  expect_equal(glm_only$rp_y$args$SL.library, "SL.glm")
  expect_equal(glm_only$rp_a$args$SL.library, "SL.glm")
  expect_equal(glm_only$gps$args$method_gps, "linear")
})

test_that("outcome_family configures only default outcome-response learners", {
  gaussian <- csdr_learners(sl_library = "SL.glm")
  binomial <- csdr_learners(
    sl_library = "SL.glm",
    outcome_family = "binomial"
  )

  expect_identical(gaussian$outcome$args$family$family, "gaussian")
  expect_identical(gaussian$rp_y$args$family$family, "gaussian")
  expect_identical(gaussian$rp_a$args$family$family, "gaussian")

  expect_identical(binomial$outcome$args$family$family, "binomial")
  expect_identical(binomial$rp_y$args$family$family, "binomial")
  expect_identical(binomial$rp_a$args$family$family, "gaussian")
  expect_error(csdr_learners(outcome_family = "poisson"), "should be one of")
})

test_that("explicit learner roles override outcome_family with one warning", {
  custom <- custom_regression(function(...) NULL, label = "custom")

  expect_warning(
    one_override <- csdr_learners(
      outcome = custom,
      outcome_family = "binomial",
      sl_library = "SL.glm"
    ),
    "ignored.*role: outcome"
  )
  expect_identical(one_override$outcome, custom)
  expect_identical(one_override$rp_y$args$family$family, "binomial")

  expect_warning(
    both_overrides <- csdr_learners(
      outcome = custom,
      rp_y = custom,
      outcome_family = "gaussian",
      sl_library = "SL.glm"
    ),
    "ignored.*roles: outcome, rp_y"
  )
  expect_identical(both_overrides$outcome, custom)
  expect_identical(both_overrides$rp_y, custom)

  expect_silent(csdr_learners(outcome = custom, sl_library = "SL.glm"))
})

test_that("learner reporting includes the SuperLearner family", {
  learners <- csdr_learners(
    outcome_family = "binomial",
    sl_library = "SL.glm"
  )
  printed <- capture.output(print(learners))
  expect_true(any(grepl("outcome:.*family: binomial", printed)))
  expect_true(any(grepl("rp_y:.*family: binomial", printed)))
  expect_true(any(grepl("rp_a:.*family: gaussian", printed)))
})

test_that("low-level and compatibility fitters still fit toy data", {
  d <- learner_toy_data()
  W <- cbind(d$A, d$C)

  fit1 <- sl_regression_fitter(d$Y, W, SL.library = "SL.glm")
  fit2 <- SL_outcome_fitter(d$Y, W, SL.library = "SL.glm")
  fit3 <- SL_nuisance_fitter(d$Y, d$C, SL.library = "SL.glm")
  gps1 <- mvn_gps_fitter(A = d$A, C = d$C, method_gps = "linear")
  gps2 <- mvn_fitter(X = d$A, C = d$C, method_gps = "linear")

  expect_s3_class(fit1, "SuperLearner")
  expect_s3_class(fit2, "SuperLearner")
  expect_s3_class(fit3, "SuperLearner")
  expect_s3_class(gps1, "mvn_inner")
  expect_s3_class(gps2, "mvn_inner")
  expect_equal(gps1$method, "linear")
  expect_equal(gps2$method, "linear")
})
