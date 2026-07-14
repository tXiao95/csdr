fit_output_data <- function() {
  data_path <- testthat::test_path("..", "..", "data", "csdr_example.rda")
  if (file.exists(data_path)) {
    data_environment <- new.env(parent = emptyenv())
    load(data_path, envir = data_environment)
    data_environment$csdr_example
  } else {
    data_environment <- new.env(parent = emptyenv())
    utils::data("csdr_example", package = "csdr", envir = data_environment)
    data_environment$csdr_example
  }
}

test_that("csdr preserves names on coefficients and scores", {
  d <- fit_output_data()
  fit <- suppressWarnings(csdr(
    d$Y, d$A, d$C,
    variants = "DR",
    d = 1,
    L = 2,
    learners = csdr_learners(sl_library = "SL.glm"),
    verbose = FALSE
  ))

  expect_identical(rownames(coef(fit)), colnames(d$A))
  expect_identical(rownames(scores(fit)), rownames(d$A))
})

test_that("csdr print and summary expose useful compact metadata", {
  d <- fit_output_data()
  fit <- suppressWarnings(csdr(
    d$Y, d$A, d$C,
    variants = c("RA", "DR"),
    d = 1,
    L = 2,
    learners = csdr_learners(sl_library = "SL.glm"),
    keep_nuisance = TRUE,
    verbose = FALSE
  ))

  printed <- capture.output(print(fit))
  expect_true(any(grepl("Exposures: 2 \\(A1, A2\\)", printed)))
  expect_true(any(grepl("Folds: 2", printed, fixed = TRUE)))
  expect_true(any(grepl("targets=yes, MAVE=yes, nuisance=yes", printed,
                        fixed = TRUE)))
  expect_true(any(grepl("RA", printed, fixed = TRUE)))
  expect_true(any(grepl("DR", printed, fixed = TRUE)))

  summarized <- capture.output(summary(fit))
  expect_true(any(grepl("Learners:", summarized, fixed = TRUE)))
  expect_true(any(grepl("Target construction:", summarized, fixed = TRUE)))
  expect_true(any(grepl("Effective folds: 2", summarized, fixed = TRUE)))
})

test_that("mave_fits returns original retained MAVE objects", {
  d <- fit_output_data()
  fit <- suppressWarnings(csdr(
    d$Y, d$A, d$C,
    variants = c("RA", "DR"),
    d = 1,
    L = 2,
    learners = csdr_learners(sl_library = "SL.glm"),
    verbose = FALSE
  ))

  dr <- mave_fits(fit, variant = "DR")
  expect_identical(dr$fit, fit$fits$DR$mave_fit)
  expect_null(dr$dimension_selection)

  all_mave <- mave_fits(fit)
  expect_named(all_mave, c("RA", "DR"))
  expect_identical(all_mave$RA$fit, fit$fits$RA$mave_fit)

  fit$control$keep_mave <- FALSE
  expect_error(mave_fits(fit), "keep_mave = TRUE")
})

test_that("mave_fits exposes raw dimension selection when it was run", {
  d <- fit_output_data()
  fit <- suppressWarnings(csdr(
    d$Y, d$A, d$C,
    variants = "RA",
    d = NULL,
    max_dim = 1,
    L = 2,
    learners = csdr_learners(sl_library = "SL.glm"),
    verbose = FALSE
  ))

  raw <- mave_fits(fit)
  expect_identical(raw$dimension_selection, fit$fits$RA$mave_dim_obj)
  expect_false(is.null(raw$dimension_selection))
})

test_that("nuisance_fits returns original named fold objects", {
  d <- fit_output_data()
  fit <- suppressWarnings(csdr(
    d$Y, d$A, d$C,
    variants = c("RA", "DR", "PO", "RP"),
    d = 1,
    L = 2,
    learners = csdr_learners(sl_library = "SL.glm"),
    keep_nuisance = TRUE,
    verbose = FALSE
  ))

  all_nuisance <- nuisance_fits(fit)
  expect_named(all_nuisance, c("outcome", "gps", "po", "rp"))
  expect_named(all_nuisance$outcome, c("fold1", "fold2"))
  expect_identical(
    all_nuisance$outcome[[1L]],
    fit$target$nuisance$outcome_models[[1L]]
  )

  gps_fold <- nuisance_fits(fit, role = "gps", fold = 1L)
  expect_named(gps_fold, "fold1")
  expect_identical(gps_fold[[1L]], fit$target$nuisance$gps_models[[1L]])
  expect_error(nuisance_fits(fit, role = "gps", fold = 3L),
               "cannot exceed.*2")
  expect_error(nuisance_fits(fit, fold = 1.5), "positive integer")

  fit$control$keep_nuisance <- FALSE
  expect_error(nuisance_fits(fit), "keep_nuisance = TRUE")
})

test_that("nuisance_fits reports roles not used by selected variants", {
  d <- fit_output_data()
  fit <- suppressWarnings(csdr(
    d$Y, d$A, d$C,
    variants = "RA",
    d = 1,
    L = 2,
    learners = csdr_learners(sl_library = "SL.glm"),
    keep_nuisance = TRUE,
    verbose = FALSE
  ))

  expect_named(nuisance_fits(fit), "outcome")
  expect_error(nuisance_fits(fit, role = "gps"), "was not used")
})
