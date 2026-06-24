end_to_end_data <- function(seed = 1, n = 60, p = 2, q = 3) {
  set.seed(seed)
  C <- matrix(rnorm(n * q), n, q)
  A <- matrix(0.5 * C[, 1] + rnorm(n * p), n, p)
  Y <- A[, 1] + 0.5 * A[, 2] + C[, 1] + rnorm(n)
  colnames(A) <- paste0("A", seq_len(p))
  colnames(C) <- paste0("C", seq_len(q))
  list(Y = Y, A = A, C = C)
}

test_that("csdr runs end-to-end for DR from the public API", {
  d <- end_to_end_data()

  fit <- suppressMessages(csdr(
    Y = d$Y,
    A = d$A,
    C = d$C,
    variants = "DR",
    d = 1,
    L = 2,
    learners = csdr_learners(sl_library = "SL.glm"),
    verbose = FALSE
  ))

  expect_s3_class(fit, "csdr_fit")
  expect_s3_class(fit$target, "csdr_target")
  expect_true("DR" %in% names(fit$fits))
  expect_true(is.matrix(coef(fit, variant = "DR")))
  expect_equal(nrow(scores(fit, variant = "DR")), nrow(d$A))
})

test_that("csdr runs end-to-end for all target variants", {
  d <- end_to_end_data(seed = 2)

  fit_all <- suppressMessages(csdr(
    Y = d$Y,
    A = d$A,
    C = d$C,
    variants = c("RA", "DR", "PO", "RP"),
    d = 1,
    L = 2,
    learners = csdr_learners(sl_library = "SL.glm"),
    verbose = FALSE
  ))

  expect_s3_class(fit_all, "csdr_fit")
  expect_equal(names(fit_all$fits), c("RA", "DR", "PO", "RP"))
  expect_equal(nrow(fit_all$summary), 4L)
  po_targets <- targets(fit_all, variant = "PO")
  expect_equal(length(po_targets$target_Y), nrow(d$A))
  expect_equal(dim(po_targets$target_A), dim(d$A))
  expect_equal(
    fit_all$summary$target_exposure[fit_all$summary$variant == "RP"],
    "residualized_A"
  )
})

test_that("csdr storage controls keep compact fit objects usable", {
  d <- end_to_end_data(seed = 3)

  fit_light <- suppressMessages(csdr(
    Y = d$Y,
    A = d$A,
    C = d$C,
    variants = "DR",
    d = 1,
    L = 2,
    learners = csdr_learners(sl_library = "SL.glm"),
    keep_targets = FALSE,
    keep_mave = FALSE,
    verbose = FALSE
  ))

  expect_null(fit_light$fits$DR$mave_fit)
  expect_null(fit_light$fits$DR$mave_dim_obj)
  expect_true(is.matrix(fit_light$fits$DR$beta))
  expect_true(is.matrix(fit_light$fits$DR$score))
  expect_error(targets(fit_light, variant = "DR"), "Targets were not retained")
})

test_that("csdr verbose flag controls learner-choice messages", {
  d <- end_to_end_data(seed = 4, n = 40)

  quiet_messages <- character()
  withCallingHandlers(
    suppressMessages(csdr(
      Y = d$Y,
      A = d$A,
      C = d$C,
      variants = "DR",
      d = 1,
      L = 2,
      learners = csdr_learners(sl_library = "SL.glm"),
      verbose = FALSE
    )),
    message = function(m) {
      quiet_messages <<- c(quiet_messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_false(any(grepl("CSDR learner choices", quiet_messages, fixed = TRUE)))

  expect_message(
    suppressWarnings(csdr(
      Y = d$Y,
      A = d$A,
      C = d$C,
      variants = "DR",
      d = 1,
      L = 2,
      learners = csdr_learners(sl_library = "SL.glm"),
      verbose = TRUE
    )),
    "CSDR learner choices"
  )
})

test_that("csdr verbose false does not emit legacy ERS row messages", {
  d <- end_to_end_data(seed = 5, n = 40)

  quiet_messages <- character()
  withCallingHandlers(
    suppressWarnings(csdr(
      Y = d$Y,
      A = d$A,
      C = d$C,
      variants = "DR",
      d = 1,
      L = 2,
      learners = csdr_learners(sl_library = "SL.glm"),
      verbose = FALSE
    )),
    message = function(m) {
      quiet_messages <<- c(quiet_messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  expect_false(any(grepl("Evaluating ERS row", quiet_messages, fixed = TRUE)))
})
