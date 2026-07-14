condition_test_data <- function() {
  load(testthat::test_path("..", "..", "data", "csdr_example.rda"))
  csdr_example
}

test_that("validation errors expose a structured parent condition", {
  d <- condition_test_data()
  d$A[1L, 1L] <- NA_real_

  condition <- tryCatch(
    csdr(d$Y, d$A, d$C, d = 1, verbose = FALSE),
    error = identity
  )

  expect_s3_class(condition, "csdr_validation_error")
  expect_s3_class(condition, "csdr_error")
  expect_identical(condition$stage, "validation")
  expect_s3_class(condition$parent, "error")
  expect_match(conditionMessage(condition), "A.*non-finite values")
})

test_that("nuisance failures report role and fold without losing the cause", {
  d <- condition_test_data()
  failing_learner <- function(...) stop("deliberate learner failure")

  condition <- tryCatch(
    csdr(
      d$Y, d$A, d$C,
      variants = "RA",
      d = 1,
      L = 2,
      learners = csdr_learners(
        outcome = custom_regression(failing_learner, label = "failing learner")
      ),
      verbose = FALSE
    ),
    error = identity
  )

  expect_s3_class(condition, "csdr_fit_error")
  expect_identical(condition$stage, "nuisance_fit")
  expect_identical(condition$fold, 1L)
  expect_identical(condition$role, "outcome")
  expect_null(condition$variant)
  expect_s3_class(condition$parent, "error")
  expect_match(conditionMessage(condition$parent), "deliberate learner failure")
})

test_that("MAVE failures identify the variant and stage", {
  d <- condition_test_data()
  condition <- tryCatch(
    suppressWarnings(csdr(
      d$Y, d$A, d$C,
      variants = "RA",
      d = 1,
      L = 2,
      learners = csdr_learners(sl_library = "SL.glm"),
      mave_control = list(method = "not-a-mave-method"),
      verbose = FALSE
    )),
    error = identity
  )

  expect_s3_class(condition, "csdr_fit_error")
  expect_identical(condition$stage, "mave_fit")
  expect_identical(condition$variant, "RA")
  expect_s3_class(condition$parent, "error")
})

test_that("Phase 1 first-user journey includes retained raw object access", {
  d <- condition_test_data()
  messages <- character()
  fit <- withCallingHandlers(
    suppressWarnings(csdr(
      d$Y, d$A, d$C,
      variants = "DR",
      d = 1,
      L = 2,
      learners = csdr_learners(sl_library = "SL.glm"),
      keep_nuisance = TRUE,
      verbose = FALSE
    )),
    message = function(message) {
      messages <<- c(messages, conditionMessage(message))
      invokeRestart("muffleMessage")
    }
  )

  expect_length(messages, 0L)
  expect_s3_class(summary(fit), "summary.csdr_fit")
  expect_identical(mave_fits(fit)$fit, fit$fits$DR$mave_fit)
  expect_identical(
    nuisance_fits(fit, role = "outcome")[[1L]],
    fit$target$nuisance$outcome_models[[1L]]
  )
})
