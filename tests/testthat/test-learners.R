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
