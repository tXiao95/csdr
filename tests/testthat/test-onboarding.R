test_that("bundled example supports the documented quick start", {
  data_path <- testthat::test_path("..", "..", "data", "csdr_example.rda")
  if (file.exists(data_path)) {
    load(data_path)
  } else {
    utils::data("csdr_example", package = "csdr", envir = environment())
  }

  expect_named(csdr_example, c("Y", "A", "C"))
  expect_length(csdr_example$Y, 60L)
  expect_equal(dim(csdr_example$A), c(60L, 2L))
  expect_equal(dim(csdr_example$C), c(60L, 3L))

  fit <- suppressMessages(csdr(
    Y = csdr_example$Y,
    A = csdr_example$A,
    C = csdr_example$C,
    variants = "DR",
    d = 1,
    L = 5,
    seed = 1,
    learners = csdr_learners(sl_library = "SL.glm"),
    verbose = FALSE
  ))

  expect_s3_class(fit, "csdr_fit")
  expect_equal(dim(coef(fit, variant = "DR")), c(2L, 1L))
  expect_equal(dim(scores(fit, variant = "DR")), c(60L, 1L))
  expect_s3_class(summary(fit), "summary.csdr_fit")
})
