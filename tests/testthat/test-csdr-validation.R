validation_data <- function() {
  set.seed(901)
  n <- 12L
  A <- matrix(rnorm(n * 2L), ncol = 2L,
              dimnames = list(paste0("row", seq_len(n)), c("A1", "A2")))
  C <- matrix(rnorm(n * 2L), ncol = 2L,
              dimnames = list(paste0("row", seq_len(n)), c("C1", "C2")))
  list(Y = rnorm(n), A = A, C = C)
}

expect_csdr_validation_before_fit <- function(data, pattern, ..., L = 2L) {
  learner_called <- FALSE
  stopper <- function(...) {
    learner_called <<- TRUE
    stop("learner should not run")
  }
  expect_error(
    csdr(
      data$Y, data$A, data$C,
      d = 1,
      L = L,
      learners = csdr_learners(outcome = custom_regression(stopper)),
      verbose = FALSE,
      ...
    ),
    pattern
  )
  expect_false(learner_called)
}

test_that("csdr validates dimensions and finite numeric data before fitting", {
  d <- validation_data()

  bad_rows <- d
  bad_rows$A <- bad_rows$A[-1L, , drop = FALSE]
  expect_csdr_validation_before_fit(bad_rows, "must have 12 rows")

  bad_type <- d
  bad_type$C <- as.data.frame(bad_type$C)
  bad_type$C$C2 <- letters[seq_len(nrow(bad_type$C))]
  expect_csdr_validation_before_fit(bad_type, "nonnumeric column: C2")

  bad_finite <- d
  bad_finite$A[2L, "A2"] <- Inf
  expect_csdr_validation_before_fit(bad_finite, "non-finite values.*A2")

  bad_y <- d
  bad_y$Y[1L] <- NA_real_
  expect_csdr_validation_before_fit(bad_y, "finite numeric vector")
})

test_that("csdr validates column names and variation before fitting", {
  d <- validation_data()

  unnamed <- d
  colnames(unnamed$A) <- NULL
  expect_csdr_validation_before_fit(unnamed, "nonempty column names")

  duplicated <- d
  colnames(duplicated$A) <- c("A1", "A1")
  expect_csdr_validation_before_fit(duplicated, "duplicated column names: A1")

  overlapping <- d
  colnames(overlapping$C)[1L] <- "A1"
  expect_csdr_validation_before_fit(overlapping, "duplicated across inputs: A1")

  constant <- d
  constant$C[, "C2"] <- 1
  expect_csdr_validation_before_fit(constant, "constant column: C2")
})

test_that("csdr validates fold and control options before fitting", {
  d <- validation_data()

  expect_csdr_validation_before_fit(d, "positive integer", L = 1.5)
  expect_csdr_validation_before_fit(d, "cannot exceed.*12", L = 13L)
  expect_csdr_validation_before_fit(
    d,
    "Unknown entry.*Did you mean 'po_marginalization'",
    target_control = list(po_marginalisation = "fold")
  )
  expect_csdr_validation_before_fit(
    d,
    "cannot override.*formula",
    mave_control = list(formula = stats::as.formula("Y ~ ."))
  )
  expect_csdr_validation_before_fit(
    d,
    "Unknown entry.*bandwidth",
    mave_control = list(bandwidth = 1)
  )
})

test_that("csdr normalization preserves supplied row and column names", {
  d <- validation_data()
  normalized <- csdr:::normalize_ers_inputs(d$Y, d$A, d$C, d$A)

  expect_identical(colnames(normalized$A), colnames(d$A))
  expect_identical(colnames(normalized$C), colnames(d$C))
  expect_identical(rownames(normalized$A), rownames(d$A))
  expect_identical(rownames(normalized$C), rownames(d$C))
})
