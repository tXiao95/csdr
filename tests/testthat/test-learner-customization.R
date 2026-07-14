learner_custom_data <- function(seed = 8101, n = 48) {
  set.seed(seed)
  C <- matrix(rnorm(n * 2), n, 2, dimnames = list(NULL, c("C1", "C2")))
  A <- matrix(rnorm(n * 2), n, 2, dimnames = list(NULL, c("A1", "A2")))
  A[, 1] <- 0.4 * C[, 1] + A[, 1]
  Y <- 1 + 0.7 * A[, 1] - 0.2 * A[, 2] + 0.3 * C[, 1] + rnorm(n, sd = 0.1)
  list(Y = Y, A = A, C = C)
}

my_lm_fitter <- function(Y, W, ...) {
  fit <- stats::lm(Y ~ ., data = data.frame(Y = Y, W))
  class(fit) <- c("lm", "my_lm_regression")
  fit
}

predict.my_lm_regression <- function(object, newdata, ...) {
  as.numeric(stats::predict(object, newdata = as.data.frame(newdata)))
}

my_gps_fitter <- function(A, C, ...) {
  mvn_gps_fitter(A = A, C = C, method_gps = "linear", ...)
}

test_that("global SL library customization is reflected in learner specs", {
  learners <- csdr_learners(sl_library = "SL.glm")

  expect_equal(learners$outcome$summary$SL.library, "SL.glm")
  expect_equal(learners$rp_y$summary$SL.library, "SL.glm")
  expect_equal(learners$rp_a$summary$SL.library, "SL.glm")
  expect_equal(learners$gps$engine, "MVN GPS")
})

test_that("role-specific SL libraries are reflected in learner specs", {
  learners <- csdr_learners(
    outcome = sl_regression(SL.library = "SL.glm"),
    rp_y = sl_regression(SL.library = "SL.glmnet"),
    rp_a = sl_regression(SL.library = "SL.glm")
  )

  expect_equal(learners$outcome$summary$SL.library, "SL.glm")
  expect_equal(learners$rp_y$summary$SL.library, "SL.glmnet")
  expect_equal(learners$rp_a$summary$SL.library, "SL.glm")
  expect_equal(learners$gps$engine, "MVN GPS")
})

test_that("custom regression learner validates and runs through csdr", {
  d <- learner_custom_data(seed = 8102)
  learners <- csdr_learners(
    outcome = custom_regression(my_lm_fitter, label = "custom lm"),
    sl_library = "SL.glm"
  )

  expect_silent(validate_csdr_learners(learners))
  fit <- suppressMessages(csdr(
    Y = d$Y,
    A = d$A,
    C = d$C,
    variants = "RA",
    d = 1,
    L = 2,
    learners = learners,
    verbose = FALSE
  ))

  expect_s3_class(fit, "csdr_fit")
  expect_equal(fit$learner_summary$label[fit$learner_summary$role == "outcome"], "custom lm")
  expect_true(is.matrix(coef(fit, variant = "RA")))
})

test_that("custom GPS learner validates and runs through csdr", {
  d <- learner_custom_data(seed = 8103)
  learners <- csdr_learners(
    gps = custom_gps(my_gps_fitter, label = "custom mvn gps"),
    sl_library = "SL.glm"
  )

  expect_silent(validate_csdr_learners(learners))
  fit <- suppressMessages(csdr(
    Y = d$Y,
    A = d$A,
    C = d$C,
    variants = "DR",
    d = 1,
    L = 2,
    learners = learners,
    verbose = FALSE
  ))

  expect_s3_class(fit, "csdr_fit")
  expect_equal(fit$learner_summary$label[fit$learner_summary$role == "gps"], "custom mvn gps")
  expect_true(is.matrix(coef(fit, variant = "DR")))
})
